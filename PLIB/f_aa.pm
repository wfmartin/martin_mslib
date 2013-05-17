package aa;

use strict;

my %aa_map = (
  A =>  { N => 1, O => 1, S => 0, C => 3, H =>  5 },
  R =>  { N => 4, O => 1, S => 0, C => 6, H => 12 },
  N =>  { N => 2, O => 2, S => 0, C => 4, H =>  6 },
  D =>  { N => 1, O => 3, S => 0, C => 4, H =>  5 },
  C =>  { N => 1, O => 1, S => 1, C => 3, H =>  5 },
  E =>  { N => 1, O => 3, S => 0, C => 5, H =>  7 },
  Q =>  { N => 2, O => 2, S => 0, C => 5, H =>  8 },
  G =>  { N => 1, O => 1, S => 0, C => 2, H =>  3 },
  H =>  { N => 3, O => 1, S => 0, C => 6, H =>  7 },
  I =>  { N => 1, O => 1, S => 0, C => 6, H => 11 },
  L =>  { N => 1, O => 1, S => 0, C => 6, H => 11 },
  K =>  { N => 2, O => 1, S => 0, C => 6, H => 12 },
  M =>  { N => 1, O => 1, S => 1, C => 5, H =>  9 },
  F =>  { N => 1, O => 1, S => 0, C => 9, H =>  9 },
  P =>  { N => 1, O => 1, S => 0, C => 5, H =>  7 },
  S =>  { N => 1, O => 2, S => 0, C => 3, H =>  5 },
  T =>  { N => 1, O => 2, S => 0, C => 4, H =>  7 },
  W =>  { N => 2, O => 1, S => 0, C =>11, H => 10 },
  Y =>  { N => 1, O => 2, S => 0, C => 9, H =>  9 },
  V =>  { N => 1, O => 1, S => 0, C => 5, H =>  9 },
);

my @aa_list = keys %aa_map;
my @elem_list = keys %{$aa_map{A}};


#-----------------------------------------------------------------------
sub get_peptide_mol_formula {
  my ($aa_seq, $without_water) = @_;

  my @aa = split //, $aa_seq;
  my %mf = map { ($_ => 0) } @elem_list;
  
  foreach my $aa (@aa) {
    $mf{$_} += $aa_map{$aa}{$_}  foreach (@elem_list);
  }

  unless ($without_water) {
    $mf{H} += 2;
    $mf{O} += 1;
  }

  \%mf;
}


#-----------------------------------------------------------------------
sub get_peptide_mol_formula_str {

  my $mol_formula = get_peptide_mol_formula(@_);

  my $str = join('', map "$_$mol_formula->{$_}", qw(C H O N));

  $str .= "S$_$mol_formula->{S}"  if $mol_formula->{S} > 0;

  $str;
}

1;
