#!/usr/bin/perl
########################################################################
#  Setup configurations in the database.
########################################################################
use strict;
use DBI;
use DBD::Pg;

use ms_db;
use opts_w_file;
use cons_matching_proteins;

my $dbh;
my %opts;
my @opts_specs = (

  # Default input file for GPM; it will be modified to produce actual input.
  'gpm_default_input=s',
 
#---------------------------------------------------------
#  For bio_context table:
#---------------------------------------------------------
  'bio_context_id=s',
  'bio_context_description=s',
  
#---------------------------------------------------------
#  For lc_configuration table:
#---------------------------------------------------------
  # Primary key
  'lc_config_id=s',
  'lc_config_description=s',

  # end time of the LC run
  'theoretical_max_rt=f',

  # maximum mz detectable by instrument
  'theoretical_max_mz=f',

  # key to calibrant_ion table
  'calibrant_ion_list_id=s',

  # half-width of mz band for calibrant ions in the exlusion list
  'calibrant_error_ppm_limit=f',

  'charge_carrier_mass=f',
  
#-----------------------------------------------------
#  For consensus_parameters table:
#-----------------------------------------------------

  #-------------------------------------------------------------------------
  # Related to FeatureLinkerUnlabeledQT program:
  #-------------------------------------------------------------------------
  # Maximum allowed mass difference in forming consensus using
  'mass_ppm_error_limit=f',

  # Maximum allowed retention time difference in forming consensus using
  'cons_max_rt_diff=f',

  # Quality threshold for consensus groups output
  'min_featurelinker_quality=f',

  # 
  'min_featurelinker_num_items=i',

  # Quantity threshold for consensus grps output from FeatureLinkerUnlabeledQT.
  'min_quantity=i',

  # Times to attempt fragmenting compound before putting in exclusion list.
  'max_frag_attempts=i',

  # Times to attempt identifying fragmented cpd before putting in excl list.
  'max_ident_attempts=i',

  'max_excl_pref_list_length=i',

  # Used in generating mass_rt_rectangle values:
  'rect_mass_ppm_width=f',
  'rect_rt_min_height=f',
  
  'max_rt_gap=f',
  'min_num_peaks=i',
  'inputs_folder=s',
  'bio_context_id=s',

#---------------------------------------------------------------------------
# For cons_matching_params table:
#---------------------------------------------------------------------------
  'cons_matching_params_key=s',
  'min_score=f',
  'proteins_fasta_file=s',
  'use_decoy!',
  'decoy_is_random!',
  'gpm_parameters_key=s',
  'gfdb_mass_ppm_tolerance=f',
  'gfdb_use_phosphorylation!',
  'gfdb_other_options=s',

);


#-----------------------------------------------------------------------
sub print_help {
  my $fmt_str = "%25s (%7s) default=%s\n";
  printf $fmt_str, 'parameter', 'type', 'default value';
  printf $fmt_str, @{$_}{qw(p type def)}   foreach (@opts_specs);
}



{ ###  MAIN  ############################################################
  opts_w_file::GetOpts(\%opts, @ms_db::db_opts, @opts_specs,
      'calibrant_ions_file=s', 'preserve_temp_files!', 'help!'
  );

  if ($opts{help}) {
     print_help();
     exit(0);
  }

  #------------------------------------------------------------------------
  #  Default values:
  #------------------------------------------------------------------------
  $opts{max_rt_gap} = 0  unless defined($opts{max_rt_gap});
  $opts{min_num_peaks}        ||= 2;

  $dbh = ms_db::open_db(%opts);
  print "Database successfully opened\n\n";

  my $query;

 
  #------------------------------------------------------------------------
  #  configuration.gpm_default_input
  #------------------------------------------------------------------------
  my $sth = $dbh->prepare(
      "SELECT value FROM configuration WHERE tag = 'gpm_default_input'");
  $sth->execute();
  if ($sth->fetchrow_hashref()) {
    print "Parameter 'gpm_default_input' is already set in the database.\n";
  }
  else {
    die "Missing required parameter 'gpm_default_input'\n"
      unless $opts{gpm_default_input};

    $sth = $dbh->prepare('INSERT INTO configuration(tag,value) VALUES(?,?)');
    $sth->execute('gpm_default_input', $opts{gpm_default_input})
        or die "Couldn't execute statement: " . $sth->errstr;

    print "Parameter 'gpm_default_input' is set to " .
          "'$opts{gpm_default_input}' in the database\n";
  }

  #------------------------------------------------------------------------
  #  bio_context
  #------------------------------------------------------------------------
  $sth = $dbh->prepare('SELECT * FROM bio_context WHERE bio_context_id = ?');
  $sth->execute( $opts{bio_context_id} );
  if ($sth->fetchrow_hashref()) {
    print "bio_context for id=$opts{bio_context_id} is already in database.\n";
  }
  elsif (not $opts{bio_context_description}) {
    die "Error: Parameter 'bio_context_description' is required.\n";
  }
  else {
    $sth = $dbh->prepare(
        'INSERT INTO bio_context(bio_context_id, description) VALUES(?,?)');
    $sth->execute(@opts{qw(bio_context_id bio_context_description)})
        or die "Couldn't execute statement: " . $sth->errstr;
    print "bio_context($opts{bio_context_id}) saved to database.\n";
  }


  #------------------------------------------------------------------------
  #  lc_configuration
  #------------------------------------------------------------------------
  $sth = $dbh->prepare('SELECT * FROM lc_configuration WHERE lc_config_id = ?');
  $sth->execute( $opts{lc_config_id} );
  if ($sth->fetchrow_hashref()) {
    print "lc_configuration for id=$opts{lc_config_id} " . 
          "is already in database.\n";
  }
  else {
    foreach (qw(lc_config_description theoretical_max_rt calibrant_ion_list_id
                charge_carrier_mass)) {
      die "Error: Parameter '$_' is required.\n" unless $opts{$_};
    }

    my @cols = qw(lc_config_id theoretical_max_rt theoretical_max_mz
                  calibrant_ion_list_id calibrant_error_ppm_limit
                  charge_carrier_mass);

    $query = 'INSERT INTO lc_configuration(description,' .
                 join(', ', @cols) . ') ' .
                'VALUES(' . join(',', ('?') x (scalar(@cols)+1)) . ')';

    $sth = $dbh->prepare($query);
    $sth->execute($opts{lc_config_description}, @opts{@cols})
        or die "Couldn't execute statement: " . $sth->errstr;
    print "lc_configuration($opts{lc_config_id}) saved to database.\n";
  }


  #------------------------------------------------------------------------
  #  Calibrant ions
  #------------------------------------------------------------------------
  my @calibrant_mz_values = ();
  if ($opts{calibrant_ions_file}) {
    die "Calibrant file ($opts{calibrant_ions_file}) doesn't exist.\n"
        unless -f $opts{calibrant_ions_file};
    die "Can't read calibrant file ($opts{calibrant_ions_file}).\n"
        unless -r $opts{calibrant_ions_file};

    open(INF, $opts{calibrant_ions_file}) or
            die "Can't open file $opts{calibrant_ions_file}.\n";
    @calibrant_mz_values = <INF>;
    close(INF);

    s/\s//gs  foreach(@calibrant_mz_values);

    printf "%d mz values read from calibrant ions file '%s'\n",
        scalar(@calibrant_mz_values), $opts{calibrant_ions_file};
  }

  $sth = $dbh->prepare(
      'SELECT count(*) FROM calibrant_ion WHERE calibrant_ion_list_id = ?');
  $sth->execute($opts{calibrant_ion_list_id})
      or die "Couldn't execute statement: " . $sth->errstr;
  my ($calibrant_count) = $sth->fetchrow_array();

  if ($calibrant_count > 0)  {
    if (scalar(@calibrant_mz_values) > 0) {
      print "Database already had $calibrant_count calibrant ions " .
            "(id = $opts{calibrant_ion_list_id}); values provided ignored.\n"
    }
    else {
      print "Database already had $calibrant_count calibrant ions " .
           "(id = $opts{calibrant_ion_list_id})\n";
    }
  }

  else {
    die 'Error: Database has no calibrant ions for list_id = ' .
            $opts{calibrant_ion_list_id} . "\n"
        if (scalar(@calibrant_mz_values) == 0);

    $query = 'INSERT INTO calibrant_ion(calibrant_ion_list_id,mz,charge) ' .
                'VALUES(?,?,?)';
    $sth = $dbh->prepare($query);
    foreach (@calibrant_mz_values) {
      $sth->execute($opts{calibrant_ion_list_id}, $_, 1)
          or die "Couldn't execute statement: " . $sth->errstr;
    }
    print scalar(@calibrant_mz_values) .
        " calibrant ions inserted into the database\n";
  }


  #------------------------------------------------------------------------
  #  consensus_parameters
  #------------------------------------------------------------------------
  $sth = $dbh->prepare('SELECT consensus_parameters_key ' .
                       'FROM consensus_parameters ' .
                       'WHERE  consensus_parameters_key = ?');
  $sth->execute($opts{consensus_parameters_key})
        or die "Couldn't execute statement: " . $sth->errstr;
  if ($sth->fetchrow_array()) {
    print 'Database already had consensus_parameters(' .
         $opts{consensus_parameters_key} . ")\n";
  }
  else {
    my @cons_param_cols =
        qw(consensus_parameters_key lc_config_id
           mass_ppm_error_limit cons_max_rt_diff 
           min_featurelinker_quality 
           max_frag_attempts max_ident_attempts
           max_excl_pref_list_length
           rect_mass_ppm_width rect_rt_min_height
          );

    push @cons_param_cols, 'min_quantity'  if ($opts{min_quantity});

    foreach (@cons_param_cols) {
      die "Missing required parameter $_.\n"  unless exists($opts{$_});
    }

    $query = sprintf('INSERT INTO consensus_parameters(%s) VALUES(%s)',
                     join(',', @cons_param_cols),
                     join(',', ('?') x scalar(@cons_param_cols))
                    );
    $sth = $dbh->prepare($query);
    $sth->execute( @opts{@cons_param_cols} )
        or die "Couldn't execute statement: " . $sth->errstr;
  }

  #------------------------------------------------------------------------
  #  cons_matching_params
  #------------------------------------------------------------------------
  $sth = $dbh->prepare('SELECT cons_matching_params_key ' .
                       'FROM cons_matching_params ' .
                       'WHERE  cons_matching_params_key = ?');
  $sth->execute($opts{cons_matching_params_key})
        or die "Couldn't execute statement: " . $sth->errstr;
  if ($sth->fetchrow_array()) {
    print 'Database already had cons_matching_params(' .
          $opts{cons_matching_params_key} . ")\n";
  }

  else {
    #-------------------------------------------------------------------
    #  First handle the proteins data.
    #-------------------------------------------------------------------
    -f $opts{proteins_fasta_file}  or
      die "Proteins fasta file ($opts{proteins_fasta_file}) doesn't exist.\n";
    -r $opts{proteins_fasta_file}  or
      die "Proteins fasta file ($opts{proteins_fasta_file}) isn't readable.\n";

    print "Processing proteins.\n";

    my $decoys_file = cons_matching_proteins::save_proteins($dbh,
        @opts{qw(cons_matching_params_key proteins_fasta_file
                 use_decoy decoy_is_random)});

    my @cons_matching_params_cols =
        qw(cons_matching_params_key
           min_score proteins_fasta_file gpm_parameters_key
           gfdb_mass_ppm_tolerance);

    foreach (@cons_matching_params_cols) {
      die "Missing required parameter $_.\n"  unless exists($opts{$_});
    }

    push(@cons_matching_params_cols, 'gfdb_other_options')
        if $opts{gfdb_other_options};

    my $query = sprintf('INSERT INTO cons_matching_params(%s) VALUES(%s)',
        join(',', @cons_matching_params_cols,
            'decoy_fasta_file', 'decoy_is_random', 'gfdb_use_phosphorylation'),
        join(',', ('?') x (scalar(@cons_matching_params_cols) + 3))
        );
    $sth = $dbh->prepare($query);
    $sth->execute( @opts{@cons_matching_params_cols}, 
        $decoys_file,
        $opts{use_decoy} ? ($opts{decoy_is_random} ? 't' : 'f') : undef,
        $opts{gfdb_use_phosphorylation} ? 't' : 'f'
        )
            or die "Couldn't execute statement: " . $sth->errstr;

  }

}
