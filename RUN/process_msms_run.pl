#!/usr/bin/perl
#-----------------------------------------------------------------------
#  Given an MS/MS run already processed through Agilent's
#  "Find Compounds by Molecular Feature", load the output files into the 
#  database to create a "mass_dataset", then process the data w/r/t a
#  "consensus" set of peaks (produced from MS-only runs).
#
#-----------------------------------------------------------------------
use strict;
use DBI;
use DBD::Pg;
#use Getopt::Long;
use POSIX;

use opts_w_file;
use ms_db;
use openms_utils;
use ms_gfdb;
use gpm;

my $dbh;

{ ###  MAIN  ###
  my %opts;
  opts_w_file::GetOpts(\%opts, @ms_db::db_opts, 'consensus_id=s', 'dataset=s',
             'inputs_folder=s', 'csv=s', 'cef=s', 'mgf=s');

  $opts{consensus_id}  or die "Missing required parameter 'consensus_id'\n";

  if ($opts{inputs_folder}) {
    -d $opts{inputs_folder} or
        die "Folder $opts{inputs_folder} does not exist.\n";
    -r $opts{inputs_folder} or
        die "Folder $opts{inputs_folder} is not readable.\n";
  }

  #-------------------------------------------------------------------
  #  The user can leave the dataset blank in the command and then be
  #  prompted for it:
  #-------------------------------------------------------------------
  if (!defined($opts{dataset})) {
    print "\n>>>>>>>>>>  ENTER DATASET >>>>\n";
    $opts{dataset} = <>;
    $opts{dataset} =~ s/^\s+//;
    $opts{dataset} =~ s/\s+$//;
  }

  #-------------------------------------------------------------------
  #  If an inputs_folder is specified, then all the input files
  #  (csv, cef, mgf) are expected to be in that folder with the dataset
  #  as base name.
  #-------------------------------------------------------------------
  if ($opts{inputs_folder}) {
    $opts{$_} ||= "$opts{inputs_folder}/$opts{dataset}.$_"
        foreach (qw(csv cef mgf));
  }
    
  #-------------------------------------------------------------------
  #  Verify that each of the input files:
  #   1) exists
  #   2) is readable
  #   3) has the correct extension
  #-------------------------------------------------------------------
  foreach (qw(csv cef mgf)) {
    $opts{$_}  or  die "Missing parameter: $_\n";
    substr($opts{$_}, -3) eq $_  or
        die "File extension (for $_ option) must be \".$_\".\n";
    -f $opts{$_} and -r $opts{$_}  or  die "Can't open/read file $opts{$_}\n";
  }

  #-------------------------------------------------------------------
  $dbh = ms_db::open_db(%opts);
  print "Database successfully opened\n\n";

  my $query = <<PARMS;
    SELECT lc_config_id, bio_context_id, cons_matching_params_key
    FROM sample_consensus
    JOIN consensus_parameters USING(consensus_parameters_key)
    WHERE consensus_id = ?
PARMS
  my $sth = $dbh->prepare($query);
  $sth->execute($opts{consensus_id})
      or die "Couldn't execute statement: " . $sth->errstr;
  my $cons_parms = $sth->fetchrow_hashref();

  $sth = $dbh->prepare('SELECT dataset FROM mass_dataset WHERE dataset = ?');
  $sth->execute($opts{dataset})
      or die "Couldn't execute statement: " . $sth->errstr;
  my $ds_row = $sth->fetchrow_hashref();
  die "Dataset $opts{dataset} already exists in the database.\n"  if $ds_row;

  #-------------------------------------------------------------------
  print "Loading MS/MS data files ($opts{dataset}).\n";
  my $sth =
      $dbh->prepare("SELECT load_molf_data(?, ?, ?, ?, ?)");
  $sth->execute($opts{dataset}, $cons_parms->{lc_config_id},
                @opts{qw(csv cef mgf)} ) 
      or die "Can't load MS/MS data: " . $sth->errstr . "\n";

  $sth = $dbh->prepare('SELECT import_molf_data(?,?)');
  $sth->execute($opts{dataset}, $cons_parms->{bio_context_id})
      or die "Couldn't execute statement: " . $sth->errstr;

  print "MS/MS data files loaded. Dataset \"$opts{dataset}\" created.\n\n";

  #-------------------------------------------------------------------
  #  Create a temp ".featureXML" file for the dataset.
  #  Use an OpenMS utility is used to adjust the retention times relative
  #  to the MS-only consensus.
  #-------------------------------------------------------------------
  print "Adjust the retention times of the MS/MS data to match consensus.\n";
  openms_utils::dataset_adjust_retention_times($dbh, $opts{dataset}, 
      "$opts{consensus_id}.featureXML");
  print "Retention values adjusted and saved to the database.\n\n";

  #-------------------------------------------------------------------
  #  Fetch matching parameters.
  #-------------------------------------------------------------------
  $query = 'SELECT proteins_fasta_file, decoy_fasta_file, min_score ' .
           'FROM cons_matching_params WHERE cons_matching_params_key = ?';
  $sth = $dbh->prepare($query);
  $sth->execute($cons_parms->{cons_matching_params_key})
      or die "Couldn't execute statement: " . $sth->errstr;
  my $cons_matching_params = $sth->fetchrow_hashref();

  #-------------------------------------------------------------------
  #  Run GPM against the dataset.
  #-------------------------------------------------------------------
  print "Run GPM.\n";
  my $gpm_match_run_id = gpm::get_gpm_matches($dbh,
      @opts{qw(gpm_parameters_key dataset mgf)}, @{$cons_matching_params}{
          qw(proteins_fasta_file decoy_fasta_file min_score)},
      $opts{preserve_gpm_output}
      );
  print "GPM completed.\n";

  print "Run MS-GFDB\n";
  my $gfdb_match_run_id = ms_gfdb::get_msgfdb_matches($dbh,
      $cons_parms->{cons_matching_params_key}, $opts{dataset}, $opts{mgf});
  print "MS-GFDB completed.\n";

  print "Combining GPM and GFDB match runs.\n";
  $sth = $dbh->prepare('SELECT combine_match_runs(?)');
  $sth->execute( [$gpm_match_run_id, $gfdb_match_run_id] )
      or die "Couldn't execute statement: " . $sth->errstr;
  my ($new_match_run_id) = $sth->fetchrow_array();

  $sth = $dbh->prepare('SELECT obs_mass_correct_from_matches(?)');
  $sth->execute($new_match_run_id)
      or die "Couldn't execute statement: " . $sth->errstr;

  $sth = $dbh->prepare('SELECT apply_matches_to_consensus(?,?)');
  $sth->execute($opts{consensus_id}, $new_match_run_id)
      or die "Couldn't execute statement: " . $sth->errstr;

  #-------------------------------------------------------------------
  #  Each time more MS/MS data is added, an exclusion list is created.
  #  The parameters guiding this are set up when creating the consensus,
  #  in order to be consistent.  It is done here to make sure it's only
  #  done once, whereas generating the report file of the exclusions that
  #  can be input to the MS instrument is done with 'exclusions_report.pl',
  #  which can be run repeatedly without side effects.
  #-------------------------------------------------------------------
  $sth = $dbh->prepare('SELECT generate_exclusion_list(?)');
  $sth->execute($opts{consensus_id})
      or die "Couldn't execute statement: " . $sth->errstr;
  my ($iter) = $sth->fetchrow_array();
  $sth->finish();

  ms_db::close_db($dbh);
}

