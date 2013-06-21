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
package opts_w_file;
use strict;

use Getopt::Long;

sub GetOpts {
  my ($h_opts, @opts) = @_;

  my %tmp_opts;
  GetOptions(\%tmp_opts, @opts, 'opts_file=s');
  
  if ($tmp_opts{opts_file}) {
     open(OPTS_FILE, $tmp_opts{opts_file})  or
         die "Can't open $tmp_opts{opts_file}\n";
     my @lines = <OPTS_FILE>;
     close(OPTS_FILE);
     s/^\s+//  foreach @lines;
     s/\s+$//s foreach @lines;
     @lines = grep { length() > 0 and !/^#/ }  @lines;
     my @no_value = grep /^\w+!$/, @lines;
     s/!$//  foreach @no_value;
 
     my @with_value = grep /^\w+\s*=/, @lines;
     %$h_opts = map {
         my ($key, $val) = (/^(\w+)\s*=\s*(.*)$/);
         $val =~ s/'([^'])+]/\1/;
         ($key, $val);
       } grep /^\w+\s*=/, @lines;

     $h_opts->{$_} = 1  foreach @no_value;
   }

   %$h_opts = (%$h_opts, %tmp_opts);
}

1;
