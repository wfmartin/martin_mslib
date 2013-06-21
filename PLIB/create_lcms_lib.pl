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
########################################################################
#  Create an LCMS library from a sample consensus.
########################################################################
use strict;
use DBI;
use DBD::Pg;

use ms_db;
use opts_w_file;

my $dbh;

{ ###  MAIN  ###
  my %opts;
  opts_w_file::GetOpts(\%opts, @ms_db::db_opts,
      'lcms_library_id=s', 'consensus_id=s', 'max_fdr=f');

  foreach (qw(lcms_library_id consensus_id)) {
    $opts{$_}  or  die "Missing parameter: $_\n";
  }

  $dbh = ms_db::open_db(%opts);
  print "Database successfully opened\n\n";

  #--------------------------------------------------------------------
  #  Make sure consensus_parameters row for the given key already exists.
  #--------------------------------------------------------------------
  my $sth =
      $dbh->prepare('SELECT * FROM lcms_lib_create_from_consensus(?,?,?)');
  $sth->execute(@opts{qw(lcms_library_id consensus_id max_fdr)} )
      or die "Couldn't execute statement: " . $sth->errstr;
  my $res = $sth->fetchrow_hashref();
  $sth->finish();

  print "$res->{num_cpds} compounds saved in LCMS library " .
        $opts{lcms_library_id} . "\n";
  print "FDR = $res->{fdr}\n";

  ms_db::close_db($dbh);
}

