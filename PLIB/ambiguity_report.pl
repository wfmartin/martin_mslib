#!/usr/bin/perl
#----------------------------------------------------------------------------
#  For a given sample consensus, generate an exclusion list. This includes:
#   1) Already identified compounds.
#   2) Compounds where the maximum number of attempts at either fragmentation
#      or identification have been reached.
#   3) All peaks (cumulative) that don't match the consensus.
#   4) Calibrant ions.
#
#  For categories #1 and #2, various charge states are iterated through.
#
#----------------------------------------------------------------------------
use strict;
use DBI;
use DBD::Pg;
use Getopt::Long;
use POSIX;
use File::Basename;

use ms_db;
use opts_w_file;

my $dbh;

{ ###  MAIN  ###
  my ($progname) = fileparse($0);
  my $usage = <<USAGE;

Usage: $progname
    --db_name       <database name>
  [ --db_user       <database user name> ]
    --consensus_id  <consensus_id>

  The $progname program generates a report (with tab-separated
  columns) about ambiguous consensus compounds -- those that
  have multiple peptides identified for them.

USAGE
  die $usage unless (scalar(@ARGV) > 0);

  my %opts;
  opts_w_file::GetOpts(\%opts, @ms_db::db_opts, 'consensus_id=s');
  $opts{consensus_id} or  die "Missing option consensus_id\n";

  $dbh = ms_db::open_db(%opts);

  my @cols = qw(match_id consensus_compound_id cw.compound_w_origin_id
      protein_name start_pos_on_protein aa_seq ng_mods_str match_run_id
      m.score cw.mass error_ppm spectrum_id);
  my $cols_str = join(',', @cols);

  my $query = <<REPORT;
    WITH ambiguous_cpd AS (
      SELECT consensus_compound_id
      FROM consensus_cpd_map
      WHERE consensus_id = ? AND NOT is_decoy
      GROUP BY consensus_compound_id
      HAVING COUNT(DISTINCT compound_w_origin_id) > 1
      )
      SELECT $cols_str ,
        (SELECT COALESCE(mgf_rt, component_rt) FROM molf_spectrum ms
         WHERE ms.dataset = mr.dataset AND ms.spectrum_id = m.spectrum_id
        ) AS component_rt,
      get_peptide_mol_formula_str(aa_seq, false) AS mol_formula
    FROM ambiguous_cpd a
    JOIN consensus_cpd_map cm USING(consensus_compound_id)
    JOIN compound_w_origin cw USING(compound_w_origin_id)
    JOIN compound co USING(compound_id)
    JOIN peptide p USING(peptide_id)
    JOIN match m ON (m.id = cm.match_id)
    JOIN match_run mr USING(match_run_id)
    WHERE NOT cm.is_decoy
    ORDER BY consensus_compound_id, compound_w_origin_id
REPORT
  my $sth = $dbh->prepare($query);
  $sth->execute($opts{consensus_id});

  s/\A\w+\.//  foreach @cols;  # remove table name prefix

  push @cols, 'component_rt', 'mol_formula';

  print join("\t", @cols), "\n";
  my $row;
  while ($row = $sth->fetchrow_hashref()) {
    print join("\t", @{$row}{@cols}), "\n";
  }

}
