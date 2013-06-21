#!/usr/bin/perl
##########################################################################
# Copyright (C) 2013 William F. Martin
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation;
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
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

  my $sth = $dbh->prepare('SELECT * FROM ms_inclusions_report(?)');
  $sth->execute($opts{consensus_id});

  print 'TargetedMSMSTable,,,,,,,', "\r\n";
  print 'On,Prec. m/z,Z,Ret. Time (min),Delta Ret. Time (min),' .
        'Iso. Width,Collision Energy,Acquisition Time (ms/spec)', "\r\n";

  my $count = 0;
  my $row;
  while (($count++) <= $opts{max_num_entries} and 
         $row = $sth->fetchrow_hashref()) {
    print join(',', 'True', @{$row}{qw(mz charge rt delta_rt)}, 
        @{$row}{qw(mz delta_rt)}, $opts{mz_width_string}, '', ''), "\r\n";
  }

}

