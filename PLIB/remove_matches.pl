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
