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

use ms_db;

my $dbh;

{ ###  MAIN  ###
  my %opts;
  GetOptions(\%opts, @ms_db::db_opts, 'consensus_id=s');
  $opts{consensus_id} or  die "Missing option consensus_id\n";

  $dbh = ms_db::open_db(%opts);

  my $sth = $dbh->prepare('SELECT consensus_merge_overlapping(?)');
  $sth->execute($opts{consensus_id})
      or die "Couldn't execute statement: " . $sth->errstr;
  my ($num_merged) = $sth->fetchrow_array();

  print "Consensus $opts{consensus_id}: merged $num_merged items\n";
}
