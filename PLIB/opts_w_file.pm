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
