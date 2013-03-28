#!/usr/bin/perl
########################################################################
#  Create a consensus from multiple (at least 3) MS-only datasets.
#
#  1) Align the retention times of all dataset (MapAlignerPoseClustering).
#  2) Form a consensus out of the aligned datasets (FeatureLinkerUnlabeledQT).
#  3) Use the list of observed masses for each consensus compound to 
#     generate a new consensus with:
#     a) Retention time is an average of contributing observed mass peaks.
#     b) Retention time range is the interval that contains the retention time
#        ranges for all the observed mass peaks.
#     c) Mass is the average of the contributing observed masses.
#
########################################################################
use strict;
use DBI;
use DBD::Pg;
use Getopt::Long;
use POSIX;
use File::Spec;

use ms_db;
use opts_w_file;
use openms_utils;

#------------------------------------------------------------------
sub load_data {
  my ($dbh, $a_datasets, $inputs_folder, $lc_config_id, $bio_context_id) = @_;

  print "Loading data:\n";

  my $sth_load = $dbh->prepare('SELECT load_molf_data(?,?,?,?,?)');
  my $sth_import = $dbh->prepare('SELECT import_molf_data(?,?)');

  my $dir = File::Spec->rel2abs($inputs_folder);

  foreach my $dataset (@$a_datasets) {
    for (qw(csv cef)) {
      my $fname = "$dir/$dataset.$_";
      -f $fname  or die "File $fname doesn't exist.\n";
      -r $fname  or die "File $fname is not readable.\n";
    }

    print "Dataset = $dataset\n";

    my $f_stem = "$dir/$dataset";
    $sth_load->execute($dataset, $lc_config_id,
            $f_stem . '.csv', $f_stem . '.cef', undef) or
        die "Error (load_molf_data)\n";

    $sth_import->execute($dataset, $bio_context_id)  or
        die "Error (import_molf_data)\n";
  }

  print "Done loading data\n\n";
}


#------------------------------------------------------------------
sub gen_featurexml_from_lcms_lib {
  my ($dbh, $lcms_library_id) = @_;

  my ($row, @rows);

  my $query = <<LCMS_READ;
    SELECT
      lcms_library_compound_id AS id,
      mass,
      (rt_start + rt_end)/2::real AS rt
    FROM lcms_library_compound
    WHERE lcms_library_id = ?
LCMS_READ
  my $sth = $dbh->prepare($query);
  $sth->execute($lcms_library_id) 
      or die "Couldn't execute statement: " . $sth->errstr;
  while ($row = $sth->fetchrow_hashref()) {
    $row->{quantity} = 100;
    push @rows, $row;
  }

  my $fname = File::Spec->rel2abs(tmpnam() . '.featureXML');
  openms_utils::create_featurexml_file('lcms_ref', $fname, \@rows);

  $fname;
}


my $dbh;

{ ###  MAIN  ###
  my %opts;
  opts_w_file::GetOpts(\%opts, @ms_db::db_opts,
      'consensus_id=s', 'consensus_parameters_key=s',
      'cons_matching_params_key=s',
      'inputs_folder=s', 'bio_context_id=s',
      'lcms_library_id=s', # optional
      'preserve_temp_files!',
      );
  my @datasets = @ARGV;

  foreach (qw(consensus_id consensus_parameters_key cons_matching_params_key
              bio_context_id)) {
    $opts{$_}  or  die "Missing parameter: $_\n";
  }

  $dbh = ms_db::open_db(%opts);
  print "Database successfully opened\n\n";

  #--------------------------------------------------------------------
  #  Make sure consensus_parameters row for the given key already exists.
  #--------------------------------------------------------------------
  my $query = <<CONS_PARMS;
    SELECT
      lc_config_id,
      cons_max_rt_diff,
      mass_ppm_error_limit,
      rect_mass_ppm_width,
      rect_rt_min_height,
      min_featurelinker_quality,
      min_featurelinker_num_items,
      min_quantity
    FROM consensus_parameters
    WHERE consensus_parameters_key = ?
CONS_PARMS

  my $sth = $dbh->prepare($query);
  $sth->execute($opts{consensus_parameters_key})
          or die "Couldn't execute statement: " . $sth->errstr;
  my $cons_parms_row = $sth->fetchrow_hashref();

  die "Database does not contain consensus_parameters(" .
        "'$opts{consensus_parameters_key}')\n"
      unless ($cons_parms_row);

  #--------------------------------------------------------------------
  $sth = $dbh->prepare('SELECT cons_matching_params_key ' .
                       'FROM cons_matching_params ' .
                       'WHERE cons_matching_params_key = ?');
  $sth->execute($opts{cons_matching_params_key})
          or die "Couldn't execute statement: " . $sth->errstr;
  my $tmp_row = $sth->fetchrow_hashref();
  die "Database does not contain cons_matching_params(" .
        "'$opts{cons_matching_params_key}')\n"
      unless ($tmp_row);


  #--------------------------------------------------------------------
  load_data($dbh, \@datasets, $opts{inputs_folder},
            $cons_parms_row->{lc_config_id}, $opts{bio_context_id})
      if ($opts{inputs_folder});


  #--------------------------------------------------------------------
  #  dataset_list first
  #--------------------------------------------------------------------
  my @opt_cols = qw(consensus_id consensus_parameters_key
                    cons_matching_params_key bio_context_id);
  push @opt_cols, 'lcms_library_id'  if $opts{lcms_library_id};
  $query = 'INSERT INTO sample_consensus(dataset_list, ' .
      join(',', @opt_cols) . ') VALUES(' .
        join(',', ('?') x (scalar(@opt_cols)+1) ) . ')';

  $sth = $dbh->prepare($query);

  my $dataset_list = '{' . join(',', map("${_}_merged", @datasets)) . '}';

  $sth->execute($dataset_list, @opts{@opt_cols})
          or die "Couldn't execute statement: " . $sth->errstr;


  #-------------------------------------------------------------------
  #  Create new datasets with similar (close) peaks merged together.
  #-------------------------------------------------------------------
  print "Creating new datasets with merged peaks.\n";
  $sth = $dbh->prepare('SELECT create_dataset_from_merged_peaks(?,?,?,?)');
  foreach (@datasets) {
    $sth->execute($_, $_ . '_merged',
                  $cons_parms_row->{mass_ppm_error_limit},
                  $cons_parms_row->{max_rt_gap})
          or die "Couldn't execute statement: " . $sth->errstr;
  }

  #-------------------------------------------------------------------
  #  Create a ".featureXML" file for each dataset.
  #-------------------------------------------------------------------
  $query = <<READ_DS;
    SELECT mass, rt, quantity, o.id
    FROM observed_mass o 
    WHERE dataset = ?
READ_DS
  $sth = $dbh->prepare($query);

  my ($row, @rows);
  foreach (@datasets) {
    my $dataset = $_ . '_merged';
    $sth->execute($dataset)
          or die "Couldn't execute statement: " . $sth->errstr;
    @rows = ();
    while ($row = $sth->fetchrow_hashref()) {
      push @rows, $row;
    }
    
    openms_utils::create_featurexml_file($dataset,
        $dataset . '.featureXML', \@rows);
  }

  #-------------------------------------------------------------------
  #  Align the retention times for the datasets.
  #  If there's an LCMS library, use its compounds for rt reference.
  #  Otherwise the default is to use one of the input datasets.
  #-------------------------------------------------------------------
  my $ref_file;
  $ref_file = gen_featurexml_from_lcms_lib($dbh, $opts{lcms_library_id})
      if $opts{lcms_library_id};

  my $cmd = 'MapAlignerPoseClustering ' .
            ' -in '  . join(' ', map($_ . '_merged.featureXML', @datasets)) .
            ' -out ' . join(' ', map($_ . '_mod.featureXML', @datasets));
  $cmd .= " -reference:file $ref_file"  if $ref_file;

  my $rv = system($cmd);
  $rv==0  or  die "Error($rv):  $cmd\n";

  unlink($ref_file)  if $ref_file;

  #-------------------------------------------------------------------
  #  Combine like peaks (across datasets) into a ".consensusXML" file.
  #  First, create a custom ".ini" file with parameters.
  #-------------------------------------------------------------------
  my $ini_file = openms_utils::gen_parms_ini_file('FeatureLinkerUnlabeledQT',
      'distance_MZ:max_difference' => $cons_parms_row->{mass_ppm_error_limit},
      'distance_MZ:unit' => 'ppm',
      'distance_RT:max_difference' => $cons_parms_row->{cons_max_rt_diff},
      );

  my $cons_file = "$opts{consensus_id}.consensusXML";
  $cmd = join(' ', 'FeatureLinkerUnlabeledQT',
              '-in', map($_ . '_mod.featureXML', @datasets),
              '-out', $cons_file,
              '-ini', $ini_file);
  my $rv = system($cmd);
  $rv==0  or  die "Error($rv):  $cmd\n";

  #-------------------------------------------------------------------
  #  Parse the consensus peaks from the output file and select those
  #  that meet the specified threshold(s).
  #-------------------------------------------------------------------
  my @cons_elems = grep {
    ((not $cons_parms_row->{min_featurelinker_quality}) or
      $_->{quality} >= $cons_parms_row->{min_featurelinker_quality}) and
    ((not $cons_parms_row->{min_quantity}) or
      $_->{quantity} >= $cons_parms_row->{min_quantity})  and
    ((not $cons_parms_row->{min_featurelinker_num_items}) or
      scalar(@{$_->{obs_mass_list}}) >=
          $cons_parms_row->{min_featurelinker_num_items})
  } openms_utils::parse_consensus_file($cons_file);

  unless ($opts{preserve_temp_files}) {
    unlink($ini_file);
    unlink($cons_file);
    foreach (@datasets) {
      unlink($_ . '_merged.featureXML');
      unlink($_ . '_mod.featureXML');
    }
  }

  #-------------------------------------------------------------------
  #  Use the consensus results as clusters of observed_mass ids.
  #  For each cluster, make a consensus_compound row.
  #  However, ignore the other data from the consensus output; compute the
  #  values for retention time, etc. by combining (usually averaging) the
  #  values from the original observed_mass rows.
  #
  #  Use the COPY command (from a file) for an efficient bulk upload.
  #-------------------------------------------------------------------
  my $tmpfile = File::Spec->rel2abs(tmpnam());
  open(OUTF, ">$tmpfile")  or  die "Can't create $tmpfile\n";
  my $counter = 1;
  foreach my $e (@cons_elems) {
    print OUTF join("\t", $counter, $_, $e->{quality}) . "\n" 
        foreach (@{$e->{obs_mass_list}});
    ++$counter;
  }
  close(OUTF);

  my $query = <<CREATE;
    CREATE TEMP TABLE tmp_cluster(
      cluster_id  integer,
      id          integer,
      quality     real
    );
CREATE
  $dbh->do($query);

  $dbh->do('CREATE INDEX tmp_clust_ind ON tmp_cluster(cluster_id)');

  $dbh->do("COPY tmp_cluster FROM '$tmpfile';");
  unlink($tmpfile)  unless ($opts{preserve_temp_files});

  $query = <<INSERT_COMPOUND;
    INSERT INTO consensus_compound(consensus_id,
        mass, rt, rt_width_at_half_ht, rt_start, rt_end,
        quantity, rel_quantity, mass_rt_rectangle, obs_mass_list,
        quality, dataset_list, msms_only, 
        min_z, max_z, dominant_z)
      WITH tmp_merged AS (SELECT * FROM generate_merged_peaks(NULL))
      SELECT
        ? AS consensus_id,
        mass, rt, rt_width_at_half_ht, rt_start, rt_end,
        quantity, rel_quantity,
        get_compound_mass_rt_rectangle(mass, rt_start, rt_end, NULL, ?, ?)
            AS mass_rt_rectangle,
        merged_ids AS obs_mass_list,

        (SELECT quality FROM tmp_cluster c
         WHERE c.cluster_id = tmp_merged.cluster_id
         LIMIT 1)  AS quality,

        (SELECT array_unique_vals (array_agg(dataset))
         FROM tmp_cluster c
         JOIN observed_mass o ON (o.id = c.id)
         WHERE c.cluster_id = tmp_merged.cluster_id)
            AS dataset_list,

         False,  
         min_z,
         max_z,
         dominant_z

      FROM tmp_merged
INSERT_COMPOUND
  $sth = $dbh->prepare($query);
  $sth->execute($opts{consensus_id},
                $cons_parms_row->{rect_mass_ppm_width},
                $cons_parms_row->{rect_rt_min_height})
          or die "Couldn't execute statement: " . $sth->errstr;

  #---------------------------------------------------------------------
  #  Get observed_mass rows that aren't included in consensus_compound rows.
  #---------------------------------------------------------------------
  $query = <<TMP_ORPHANS;
    CREATE TEMP TABLE tmp_orphan AS
      WITH pass_1 AS (
        SELECT id
        FROM observed_mass
        WHERE dataset IN
          (SELECT unnest(dataset_list)
           FROM sample_consensus WHERE consensus_id = ?)
        EXCEPT
        SELECT unnest(obs_mass_list) AS id
        FROM consensus_compound
        WHERE consensus_id = ?
        )
      SELECT
        id AS obs_mass_id,
        mass,
        rt_start,
        rt_end,
        get_compound_mass_rt_rectangle(mass, rt_start, rt_end,
                rt_adjustment_to_consensus, ?, ?)
            AS mass_rt_rectangle,
        min_z,
        max_z,
        dominant_z
      FROM pass_1 JOIN observed_mass USING(id)
TMP_ORPHANS
  $sth = $dbh->prepare($query);
  $sth->execute( $opts{consensus_id}, $opts{consensus_id}, 
                $cons_parms_row->{rect_mass_ppm_width},
                $cons_parms_row->{rect_rt_min_height})
      or die "Couldn't execute statement: " . $sth->errstr;

  #-------------------------------------------------------------------
  #  Special handling if this sample has a LCMS library for reference.
  #-------------------------------------------------------------------
  if ($opts{lcms_library_id}) {

    $query = <<INSERT_LCMS_MAP;
      INSERT INTO consensus_cpd_lcms_map(consensus_id,
          consensus_compound_id, lcms_library_compound_id)
        SELECT ?, cc.consensus_compound_id, lc.lcms_library_compound_id
        FROM consensus_compound cc, lcms_library_compound  lc
        WHERE cc.consensus_id = ? AND 
              cc.mass_rt_rectangle && lc.mass_rt_rectangle
INSERT_LCMS_MAP
    $sth = $dbh->prepare($query);
    $sth->execute($opts{consensus_id}, $opts{consensus_id})
        or die "Couldn't execute statement: " . $sth->errstr;


    #-------------------------------------------------------------------
    #  Find unclustered peaks that didn't cluster but that match a compound
    #  in LCMS library that isn't already matched with a consensus compound.
    #
    #  Assume there are few enough of those that there's no efficiency
    #  issue of searching for overlaps again after inserting to
    #  consensus_compound.
    #-------------------------------------------------------------------
    $query = <<LCMS_ORPHAN;
      CREATE TEMP TABLE tmp_orphan_to_promote AS
        --  LCMS compounds not already matched to consensus compounds:
        WITH lcms_unmatched AS (
          SELECT lcms_library_compound_id
          FROM lcms_library_compound
          WHERE lcms_library_id = ?
            EXCEPT
          SELECT lcms_library_compound_id
          FROM consensus_cpd_lcms_map
          WHERE consensus_id = ?
          )
        SELECT DISTINCT o.obs_mass_id, lcms_library_id
        FROM lcms_unmatched
        JOIN lcms_library_compound lc USING(lcms_library_compound_id),
             tmp_orphan o
        WHERE o.mass_rt_rectangle && lc.mass_rt_rectangle
LCMS_ORPHAN
    $sth = $dbh->prepare($query);
    $sth->execute($opts{lcms_library_id}, $opts{consensus_id})
        or die "Couldn't execute statement: " . $sth->errstr;
    
    $query = <<PROMOTE;
      INSERT INTO consensus_compound(consensus_id,
          mass, rt, rt_width_at_half_ht, rt_start, rt_end,
          quantity, rel_quantity, mass_rt_rectangle, obs_mass_list,
          dataset_list, msms_only, min_z, max_z, dominant_z)
        SELECT
          ? AS consensus_id,
          mass, rt, rt_width_at_half_ht, rt_start, rt_end,
          quantity, rel_quantity,
          get_compound_mass_rt_rectangle(mass, rt_start, rt_end, NULL, ?, ?)
              AS mass_rt_rectangle,
          ARRAY[obs_mass_id] AS obs_mass_list,
          ARRAY[dataset] AS dataset_list,
          False,  
          min_z,
          max_z,
          dominant_z
        FROM tmp_orphan_to_promote p
        JOIN observed_mass o ON (o.id = p.obs_mass_id)
PROMOTE
    $sth = $dbh->prepare($query);
    $sth->execute($opts{consensus_id},
                  $cons_parms_row->{rect_mass_ppm_width},
                  $cons_parms_row->{rect_rt_min_height})
        or die "Couldn't execute statement: " . $sth->errstr;

    $query = <<DELETE;
      DELETE FROM tmp_orphan o
      USING tmp_orphan_to_promote p
      WHERE o.obs_mass_id = p.obs_mass_id
DELETE
    $dbh->do($query)  or  die "Error: $query\n";

    $dbh->do('DROP TABLE tmp_orphan_to_promote') or die "Error: $query\n";

    $query = <<UPDATE_LCMS;
      UPDATE consensus_compound cc
      SET exclude_compound = True
      FROM consensus_cpd_lcms_map m
      WHERE m.consensus_id = ? AND
            cc.consensus_compound_id = m.consensus_compound_id
UPDATE_LCMS
    $sth = $dbh->prepare($query);
    $sth->execute($opts{consensus_id})
        or die "Couldn't execute statement: " . $sth->errstr;
  }

  #-------------------------------------------------------------------
  #  Create a ".featureXML" file from the consensus, to be used as a 
  #  reference file later in MapAlignerPoseClustering to align the
  #  retention times of subsequent MS/MS runs with the consensus.
  #-------------------------------------------------------------------
  my $cons_file = "$opts{consensus_id}.featureXML";
  openms_utils::create_featurexml_file($opts{consensus_id}, $cons_file,
                                       \@cons_elems);
 
  openms_utils::dataset_adjust_retention_times($dbh, "${_}_merged", $cons_file)
      foreach (@datasets);

  $query = <<ORPHANS;
    INSERT INTO obs_mass_no_cons(obs_mass_id, consensus_id, mass,
        rt_start, rt_end, mass_rt_rectangle, min_z, max_z, dominant_z)
      SELECT obs_mass_id, ?, mass,
            rt_start, rt_end, mass_rt_rectangle, min_z, max_z, dominant_z
      FROM tmp_orphan;
ORPHANS
  $sth = $dbh->prepare($query);
  $sth->execute($opts{consensus_id})
      or die "Couldn't execute statement: " . $sth->errstr;

  ms_db::close_db($dbh);
}

