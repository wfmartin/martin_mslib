package cons_matching_proteins;

use strict;
use POSIX;
use Cwd;
use DBI;
use DBD::Pg;
use List::Util('shuffle');
use File::Spec;


#-------------------------------------------------------------------------
sub gen_decoy_file_name {
  my ($cons_matching_params_key) = @_;

  $cons_matching_params_key =~ s/ /_/g;
  $cons_matching_params_key =~ s/\W//g;

  my $i;
  while (-f "$cons_matching_params_key$i.fasta") {
    ++$i;
  }

  File::Spec->rel2abs("$cons_matching_params_key$i.fasta");
}



#-------------------------------------------------------------------------
sub get_decoy_seq {
  my ($seq, $decoy_is_random) = @_;

  my $mod_seq;

  if ($decoy_is_random eq 't') {
    my @chars = split //, $seq;
    my @shuffled = @chars[ shuffle(0..$#chars) ];
    $mod_seq = join('', @shuffled);
  }
  
  else {
    $mod_seq = reverse($seq);
  }

  $mod_seq;
}


#-------------------------------------------------------------------------
#
#-------------------------------------------------------------------------
sub save_proteins {
  my ($dbh, $cons_matching_params_key, $proteins_fasta_file, 
      $use_decoy, $decoy_is_random) = @_;

  my $decoy_file;

  if ($use_decoy) {
    $decoy_file = gen_decoy_file_name($cons_matching_params_key);
    open(DECOY, ">$decoy_file")  or die "Can't create $decoy_file.\n";
  }

  my $sth = $dbh->prepare('INSERT INTO cons_matching_protein(' .
    'cons_matching_params_key, protein_name, decoy, aa_seq) VALUES(?,?,?,?)');

  local($/);
  open(INF, $proteins_fasta_file)  or
      die("Can't open file $proteins_fasta_file");
  my $proteins_text = <INF>;
  close(INF);
  $proteins_text =~ s/\r//g;

  my ($ignore, @sections) = split /^>/m, $proteins_text;
  foreach my $sect (@sections) {
    my $nl_pos = index($sect, "\n");
    my $header = substr($sect, 0, $nl_pos);
    my $seq =    substr($sect, $nl_pos+1);
    $seq =~ s/\s//g;

    my ($protein_name) = ($header =~ /\A(?:sp|tr)\|\w+\|(\w+) /);
    ($protein_name) = ($header =~ /\A(\w+)/)  unless ($protein_name);

    $sth->execute($cons_matching_params_key, $protein_name, 0, $seq)
       or die "Couldn't execute statement: " . $sth->errstr;

    if ($use_decoy) {
      my $decoy_aa_seq = get_decoy_seq($seq, $decoy_is_random);

      my $decoy_header = $header;
      $decoy_header =~ s/^(\S+)/\1_DECOY/;

      $sth->execute($cons_matching_params_key, $protein_name . '_DECOY',
                    1, $decoy_aa_seq)
         or die "Couldn't execute statement: " . $sth->errstr;

      print DECOY '>', $decoy_header, "\n";
      print DECOY join("\n", grep $_, split /(.{60})/, $decoy_aa_seq), "\n";
    }
  }

  close(DECOY)  if ($use_decoy);
  
  $decoy_file;
}

1;
