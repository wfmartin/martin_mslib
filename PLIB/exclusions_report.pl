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

use ms_db;
use opts_w_file;

my $dbh;

{ ###  MAIN  ###
  my %opts;
  opts_w_file::GetOpts(\%opts, @ms_db::db_opts,
    'consensus_id=s', 'mz_width_string=s');
  $opts{consensus_id} or  die "Missing option consensus_id\n";
  $opts{mz_width_string} ||= 'Narrow (~1.3 m/z)';

  $dbh = ms_db::open_db(%opts);

  my $sth = $dbh->prepare('SELECT * FROM ms_exclusions_report(?)');
  $sth->execute($opts{consensus_id});

  print 'AutoPreferredExcludeMSMSTable,,,,,,,,', "\r\n";
  print 'On,Prec. m/z,Delta m/z (ppm),Z,Prec. Type,Ret. Time (min),' .
        'Delta Ret. Time (min),Iso. Width,Collision Energy', "\r\n";

  my $row;
  while ($row = $sth->fetchrow_hashref()) {
    print join(',', 'True', @{$row}{qw(mz delta_mz_ppm charge)}, 'Exclude',
               @{$row}{qw(rt delta_rt)}, $opts{mz_width_string}, ''), "\r\n";
  }

}

