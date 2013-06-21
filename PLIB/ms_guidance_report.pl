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
#  For a given sample consensus, generate a excluded/preferred list.
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

  $dbh = ms_db::open_db(%opts);

  my $sth = $dbh->prepare('SELECT * FROM msms_guidance_report(?)');
  $sth->execute($opts{consensus_id}) 
      or die "Couldn't execute statement: " . $sth->errstr;

  print 'AutoPreferredExcludeMSMSTable,,,,,,,,', "\r\n";
  print 'On,Prec. m/z,Delta m/z (ppm),Z,Prec. Type,Ret. Time (min),' .
        'Delta Ret. Time (min),Iso. Width,Collision Energy', "\r\n";

  my @cols = qw(mz delta_mz_ppm charge excl_or_pref rt delta_rt);

  my $row;
  while ($row = $sth->fetchrow_hashref()) {
    $row->{excl_or_pref} = $row->{exclude} ? 'Exclude' : 'Preferred';

    print join(',', 'True', @{$row}{@cols}, $opts{mz_width_string}, ''),
               "\r\n";
  }

}
