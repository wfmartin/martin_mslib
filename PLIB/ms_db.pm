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
package ms_db;
use strict;
use DBI;
use DBD::Pg;

our @db_opts = ('db_name=s', 'db_user=s', 'trace=s');

my $trace;

#--------------------------------------------------------------------------
sub open_db {
  my %opts = @_;

  $trace = $opts{trace};

  my $db_user = $opts{db_user} || $ENV{USER};
  my $db_name = $opts{db_name} || $ENV{MS_DB_NAME};
  $db_name   or  die "Missing option db_name\n";

  my $dbh =
      DBI->connect("DBI:Pg:database=$db_name;host=localhost", $db_user, '')
        or die "Unable to connect: $DBI::errstr\n";

  if ($trace) {
    $dbh->do("SET track_functions='all'");
    $dbh->do("SELECT pg_stat_reset()");
  }

  $dbh;
}


#--------------------------------------------------------------------------
sub functions_report {
  my $dbh = shift;

  open(TRACE_OUT, ">$trace")  or  die "Can't create trace file $trace.\n";

  my @cols = qw(funcname calls total_time self_time);
  print TRACE_OUT join("\t", @cols), "\n";

  my $sth = $dbh->prepare(
      'SELECT * FROM pg_stat_user_functions ORDER BY self_time DESC');
  $sth->execute();

  my $row;
  while ($row = $sth->fetchrow_hashref()) {
    print TRACE_OUT join("\t", @{$row}{@cols}), "\n";
  }
  close(TRACE_OUT);
}


#--------------------------------------------------------------------------
sub close_db {
  my $dbh = shift;

  functions_report($dbh)  if ($trace);

  $dbh->disconnect();
}

1;
