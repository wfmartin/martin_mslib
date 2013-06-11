#!/usr/bin/perl
#----------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------
use strict;
use DBI;
use DBD::Pg;
use Getopt::Long;

use ms_db;

my $dbh;

{ ###  MAIN  ###
  my ($progname) = fileparse($0);
  my $usage = <<USAGE;

Usage: $progname
    --db_name       <database name>
  [ --db_user       <database user name> ]
    --consensus_id  <consensus_id>
    --max_fdr       <false discovery rate limit>

  The $progname program sorts the consensus compound identifications in order
  of descending p-values.  The list is traversed, counting real compound
  identifications and decoy compound identifications until the proportion of
  decoy hits encountered exceeds the allowable false discovery rate.
  Real hits are put into the table 'final_consensus_compound'.

USAGE
  die $usage unless (scalar(@ARGV) > 0);

  my %opts;
  GetOptions(\%opts, @ms_db::db_opts, 'consensus_id=s', 'min_rt_width=f');
  $opts{consensus_id} or  die "Missing option consensus_id\n";
  $opts{min_rt_width} or  die "Missing option min_rt_width\n";

  $dbh = ms_db::open_db(%opts);

  my $sth = $dbh->prepare('SELECT finalize_consensus_cpds(?,?)');
  $sth->execute($opts{consensus_id}, $opts{max_fdr})
      or die "Couldn't execute statement: " . $sth->errstr;
}
