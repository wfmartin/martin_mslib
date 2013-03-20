package ms_gfdb;
use strict;
use File::Spec;
use POSIX;

my $msgfdb_fmt_str =
  'java -Xmx2000M -jar MSGFDB.jar -s %s -d %s -t %fppm -o %s %s';

my %spec_map;
my ($sth_peptide_pos, $sth_cpd);

#--------------------------------------------------------------------------
#  Since the output from ms-gfdb refers to the spectra it matches by an index
#  offset into the mgf file, create a mapping between index and obs_mass_id.
#--------------------------------------------------------------------------
sub create_spec_map {
  my ($dbh, $dataset) = @_;

  my $query = <<SPEC_MAP;
    SELECT DISTINCT mgf_index, spectrum_id, om.id AS obs_mass_id
    FROM molf_spectrum mf
    JOIN observed_mass om
      ON (mf.dataset = om.dataset AND mf.cpd = om.peak_id)
    WHERE mf.dataset = '$dataset'
SPEC_MAP
  my $sth = $dbh->prepare($query);
  $sth->execute();
  my $row;
  while ($row = $sth->fetchrow_hashref()) {
    $spec_map{$row->{mgf_index}} = [$row->{obs_mass_id}, $row->{spectrum_id}];
  }
}


#--------------------------------------------------------------------------
#  Returns:  peptide_seq, position, ng_mods_str
#--------------------------------------------------------------------------
sub parse_peptide {
   my ($cons_matching_params_key, $protein_name, $seq_col, $is_decoy) = @_;

   my ($aa_left, $seq, $aa_right) = ($seq_col =~ /\A([A-Z_])\.(.*)\.([A-Z_])/);
   $aa_left = ''   if ($aa_left eq '_');
   $aa_right = ''  if ($aa_right eq '_');

   my @seq_parts = split /\+(\d+\.\d+)/, $seq;
   my $num_mods = (scalar(@seq_parts)-1)/2;

   my $peptide_seq = $seq_parts[0];

   my @mods = ();
   local($_);
   foreach my $i (0 .. $num_mods-1) {
     $peptide_seq .= $seq_parts[2*$i]  if ($i>0);
     $_ = $seq_parts[2*$i+1];
     my $mod_label = /^79\.96/ ? 'phosphate' :
                      (/^15\.9/ ? 'oxidation' :
                        (/^42\.01/ ? 'acetate' : die "unknown mod: $_")
                      );
     push @mods, "$mod_label\@" . length($peptide_seq);
   }
  
   my $seq_with_flanking = $aa_left . $peptide_seq . $aa_right;

   $sth_peptide_pos->execute($cons_matching_params_key, $protein_name,
       $is_decoy, $seq_with_flanking);

   my ($seq_pos) = $sth_peptide_pos->fetchrow_array();
   $seq_pos += length($aa_left);

  my $ng_mods_str = join(',', @mods);

  ($peptide_seq, $seq_pos, $ng_mods_str);
}


#--------------------------------------------------------------------------
sub parse_output {
  my ($cons_matching_params_key, $is_decoy, $out_file) = @_;

  open(INF, $out_file)  or  die "Can't open $out_file for read.\n";
  do { $_ = <INF> } until /^#SpecFile/;
  my @cols = split/\t/;
  my @lines = <INF>;
  chomp(@lines);
  close(INF);

  map {
    my $in_rec = {};
    @{$in_rec}{@cols} = split /\t/;

    my $out_rec = { };
    @{$out_rec}{'obs_mass_id', 'spectrum_id'} =
        @{$spec_map{$in_rec->{SpecIndex}}};

    $out_rec->{score} = -log( $in_rec->{'P-value'} )/log(10);

    my ($protein_name) = ($in_rec->{Protein} =~ /\A(?:sp|tr)\|\w+\|(\w+)/);
    ($protein_name) = ($in_rec->{Protein} =~ /\A(\w+)/)  unless $protein_name;

    my ($peptide_seq, $position, $ng_mods_str) =
        parse_peptide($cons_matching_params_key, $protein_name,
                      $in_rec->{Peptide}, $is_decoy);

    $sth_cpd->execute($protein_name, $peptide_seq, $position, $ng_mods_str);
    $out_rec->{compound_w_origin_id} = $sth_cpd->fetchrow_array();

    $out_rec;
  } splice(@lines, 0, scalar(@lines)-2);
}


#--------------------------------------------------------------------------
#  Parse the hits and create a file in a format for the db COPY command.
#--------------------------------------------------------------------------
sub process_hits {
  my ($dbh, $dataset, $cons_matching_params_key, $min_score,
      $gfdb_mass_ppm_tolerance,
      $normal_outf, $decoys_outf) = @_;

  my $outf_name = File::Spec->rel2abs(tmpnam());

  my @normal_hits = grep { $_->{score} >= $min_score }
      parse_output($cons_matching_params_key, 0, $normal_outf);
  $_->{is_decoy} = 'f'  foreach (@normal_hits);

  my @decoy_hits  = grep { $_->{score} >= $min_score }
      parse_output($cons_matching_params_key, 1, $decoys_outf);
  $_->{is_decoy} = 't'  foreach (@decoy_hits);


  my $query = <<INSERT_MATCH_RUN;
    INSERT INTO match_run(dataset, match_method, match_params,
        score_method, error_ppm_threshold,
        num_matches, num_decoy_matches)
    VALUES(?,?,?,?,?,?,?)
    RETURNING match_run_id
INSERT_MATCH_RUN
  my $sth = $dbh->prepare($query);
  $sth->execute($dataset, 'MS-GFDB', $cons_matching_params_key,
                '-log(expect)', $gfdb_mass_ppm_tolerance, 
                scalar(@normal_hits), scalar(@decoy_hits) )
      or die "Database error: " . $sth->errstr . "\n";
  my ($match_run_id) = $sth->fetchrow_array();

  open(HITS_OUTF, ">$outf_name")  or  die "Can't create $outf_name\n";
  my @cols = qw(obs_mass_id compound_w_origin_id score spectrum_id is_decoy);

  print HITS_OUTF join("\t", @{$_}{@cols}), "\n"  foreach (@normal_hits);
  print HITS_OUTF join("\t", @{$_}{@cols}), "\n"  foreach (@decoy_hits);
  close(HITS_OUTF);

  $query = <<CREATE_TMP_MATCH;
    CREATE TEMP TABLE tmp_match(
      obs_mass_id            integer,
      compound_w_origin_id   integer,
      score                  real,
      spectrum_id            integer,
      is_decoy               boolean
    );
CREATE_TMP_MATCH
  $dbh->do($query);

  $dbh->do('COPY tmp_match(' . join(',', @cols) . ") FROM '$outf_name'");

  $query = 'INSERT INTO match(match_run_id,' . join(',', @cols) .
    ', misc) ' .
    "SELECT $match_run_id, " . join(',', @cols) .
    ", 'GFDB' FROM tmp_match";
  $dbh->do($query)  or  die "Database error, creating tmp_match:\n$query\n"; 

  $dbh->do('DROP TABLE tmp_match');

  $match_run_id;
}


#--------------------------------------------------------------------------
#  Run ms-gfdb against the proteins database and against the decoy database
#  if one is provided.
#
#  Return the match_run_id.
#--------------------------------------------------------------------------
sub get_msgfdb_matches {
  my ($dbh, $cons_matching_params_key, $dataset, $mgf_in_file) = @_;

  create_spec_map($dbh, $dataset);
  $sth_peptide_pos = $dbh->prepare('SELECT get_peptide_position(?,?,?,?)');

  $sth_cpd = $dbh->prepare('SELECT get_compound_for_match(?,?,?,?)');

  my $normal_outf = File::Spec->rel2abs(tmpnam());
  my $decoys_outf;

  my $query = <<CONS_MATCH;
    SELECT
      min_score,
      proteins_fasta_file,
      decoy_fasta_file,
      gfdb_mass_ppm_tolerance,
      gfdb_use_phosphorylation,
      gfdb_other_options
    FROM cons_matching_params
    WHERE cons_matching_params_key = ?
CONS_MATCH

  my $sth = $dbh->prepare($query);
  $sth->execute($cons_matching_params_key)
      or die "Database error: " . $sth->errstr . "\n";
  my $cons_match_row = $sth->fetchrow_hashref();

  my $msgfdb_fmt_str =
    'java -Xmx2000M -jar MSGFDB.jar -s %s -d %s -t %fppm -o %s %s';

  my $opts_at_end = join(' ', 
      ($cons_match_row->{gfdb_use_phosphorylation} eq 't') ? '-protocol 1':'',
      $cons_match_row->{gfdb_other_options});

  #--------------------------------------------------------------------
  #  Run with real proteins.
  #--------------------------------------------------------------------
  my $cmd = sprintf($msgfdb_fmt_str, $mgf_in_file,
      $cons_match_row->{proteins_fasta_file},
      $cons_match_row->{gfdb_mass_ppm_tolerance},
      $normal_outf, $opts_at_end);

  print "Running msgfdb (proteins):\n$cmd\n\n";
  system($cmd) == 0  or die "Error running msgfdb: $cmd\n";

  #--------------------------------------------------------------------
  #  Run with decoy proteins if decoys file provided.
  #  Don't use the built-in decoy mechanism.
  #--------------------------------------------------------------------
  if ($cons_match_row->{decoy_fasta_file}) {

    $decoys_outf = File::Spec->rel2abs(tmpnam());

    $cmd = sprintf($msgfdb_fmt_str, $mgf_in_file,
        $cons_match_row->{decoy_fasta_file},
        $cons_match_row->{gfdb_mass_ppm_tolerance},
        $decoys_outf, $opts_at_end);

    print "Running msgfdb (decoys):\n$cmd\n\n";
    system($cmd) == 0  or die "Error running msgfdb: $cmd\n";
  }

  my $match_run_id = process_hits($dbh, $dataset, $cons_matching_params_key,
          $cons_match_row->{min_score},
          $cons_match_row->{gfdb_mass_ppm_tolerance},
          $normal_outf, $decoys_outf);

  $sth_peptide_pos->finish();
  $sth_cpd->finish();

  $match_run_id;
}

1;
