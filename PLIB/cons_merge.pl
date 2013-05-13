#!/usr/bin/perl
#----------------------------------------------------------------------------
#  For a given sample consensus, merge consensus_compound values with 
#
#----------------------------------------------------------------------------
use strict;
use DBI;
use DBD::Pg;
use Getopt::Long;
use POSIX;
use File::Basename;

use ms_db;

my $dbh;

{ ###  MAIN  ###
  my ($progname) = fileparse($0);
  my $usage = <<USAGE;

Usage: $progname
    --db_name       <database name>
  [ --db_user       <database user name> ]
    --consensus_id  <consensus_id>

  The $progname program merges consensus compounds 
  columns) about ambiguous consensus compounds -- those that
  have multiple peptides identified for them.

USAGE
  die $usage unless (scalar(@ARGV) > 0);

  my %opts;
  GetOptions(\%opts, @ms_db::db_opts, 'consensus_id=s, min_rt_width=f');
  $opts{consensus_id} or  die "Missing option consensus_id\n";
  $opts{min_rt_width} or  die "Missing option min_rt_width\n";

  $dbh = ms_db::open_db(%opts);

  my $sth = $dbh->prepare('SELECT consensus_merge_overlapping(?,?)');
  $sth->execute($opts{consensus_id}, $opts{min_rt_width})
      or die "Couldn't execute statement: " . $sth->errstr;
  my ($num_merged) = $sth->fetchrow_array();

  print "Consensus $opts{consensus_id}: merged $num_merged items\n";
}
