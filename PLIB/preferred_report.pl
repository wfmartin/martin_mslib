#!/usr/bin/perl
#----------------------------------------------------------------------------
#  For a given sample consensus, generate a preferred list. This includes:
#   1) Compounds where the maximum number of attempts at either fragmentation
#      or identification have not been reached.
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
    'consensus_id=s', 'mz_width_string=s', 'max_num_entries=i');
  $opts{consensus_id} or  die "Missing option consensus_id\n";
  $opts{mz_width_string} ||= 'Narrow (~1.3 m/z)';
  $opts{max_num_entries} ||= 9999;

  $dbh = ms_db::open_db(%opts);

  my $query = <<LC;
    SELECT lc_config_id 
    FROM consensus_parameters cp
    WHERE consensus_parameters_key = (
      SELECT consensus_parameters_key 
      FROM sample_consensus
      WHERE consensus_id = ?
      )
LC
  my $sth = $dbh->prepare($query);
  $sth->execute($opts{consensus_id});
  my ($lc_config_id) = $sth->fetchrow_array();

  #----------------------------------------------------------------------
  #  Exclude calibrant ions.
  #----------------------------------------------------------------------
  $sth = $dbh->prepare('SELECT * FROM calibrant_ions_report(?)');
  $sth->execute($lc_config_id);
  
  print 'AutoPreferredExcludeMSMSTable,,,,,,,,', "\r\n";
  print 'On,Prec. m/z,Delta m/z (ppm),Z,Prec. Type,Ret. Time (min),' .
        'Delta Ret. Time (min),Iso. Width,Collision Energy', "\r\n";

  my $row;
  while ($row = $sth->fetchrow_hashref()) {
    print join(',', 'True', @{$row}{qw(mz delta_mz_ppm charge)}, 'Exclude',
               @{$row}{qw(rt delta_rt)}, $opts{mz_width_string}, ''), "\r\n";
  }

  $sth = $dbh->prepare('SELECT * FROM ms_inclusions_report(?)');
  $sth->execute($opts{consensus_id});

  my $count = 0;
  my $row;
  while (($count++) <= $opts{max_num_entries} and 
         $row = $sth->fetchrow_hashref()) {
    print join(',', 'True', @{$row}{qw(mz delta_mz_ppm charge)}, 'Preferred',
               @{$row}{qw(rt delta_rt)}, $opts{mz_width_string}, ''), "\r\n";
  }

}

#   print join(',', 'True', @{$row}{qw(mz charge rt delta_rt)}, 
#       @{$row}{qw(mz delta_rt)}, $opts{mz_width_string}, '', ''), "\r\n";
