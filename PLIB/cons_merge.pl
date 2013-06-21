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
    --min_rt_width  <minimum retention time width>

  The $progname program merges consensus compounds that have the same
  compound identification with overlapping retention times.
  Consensus compounds whose retention time width is less than 'min_rt_width'
  are widened for the purpose of testing for overlaps.

USAGE
  die $usage unless (scalar(@ARGV) > 0);

  my %opts;
  GetOptions(\%opts, @ms_db::db_opts, 'consensus_id=s', 'min_rt_width=f');
  $opts{consensus_id} or  die "Missing option consensus_id\n";
  $opts{min_rt_width} or  die "Missing option min_rt_width\n";

  $dbh = ms_db::open_db(%opts);

  my $sth = $dbh->prepare('SELECT consensus_merge_overlapping(?,?)');
  $sth->execute($opts{consensus_id}, $opts{min_rt_width})
      or die "Couldn't execute statement: " . $sth->errstr;
  my ($num_merged) = $sth->fetchrow_array();

  print "Consensus $opts{consensus_id}: merged $num_merged items\n";
}
