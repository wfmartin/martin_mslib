#!/usr/bin/perl
#----------------------------------------------------------------------------
#  For a given sample consensus, remove the given matches.
#
#  STDIN should contain two columns with consensus_compound_id, match_id
#  (separated by a comma, with no header)
#----------------------------------------------------------------------------
use strict;
use DBI;
use DBD::Pg;
use Getopt::Long;
use POSIX;

use ms_db;
use opts_w_file;

my $dbh;

{ ###  MAIN  ###
  my %opts;
  opts_w_file::GetOpts(\%opts, @ms_db::db_opts, 'consensus_id=s');
  $opts{consensus_id} or  die "Missing option consensus_id\n";

  my @data;
  while ($_ = <>) {
    s/\A\s+//g;
    s/\s+\Z//g;
    push @data, [ split /,/ ];
  }
  scalar(@data) > 0  or  die "Removed items expected in STDIN\n";

  $dbh = ms_db::open_db(%opts);

  my $sth = $dbh->prepare('SELECT consensus_remove_match(?,?,?)');
  $sth->execute($_->[0], $opts{consensus_id}, $_->[1])  foreach (@data);
}
