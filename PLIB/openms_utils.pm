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
package openms_utils;
use strict;
use POSIX;
use DBI;
use DBD::Pg;
use XML::Writer;
use XML::Parser;
use IO::File;

#---------------------------------------------------------------------------
#  Creates a modified parameters (.ini) file.
#
#  Given parameter values passed in a hash with key/values as follows:
#   "node:item"          (as hash key) 
#      where "node" is the "name" attribute of a NODE element
#            "item" is the "name" attribute of an ITEM element under the node.
#   "value" attribute    (as hash value)
#
#---------------------------------------------------------------------------
sub gen_parms_ini_file {
  my ($prog_name, %values) = @_;

  my $ini_file = tmpnam();
  my $cmd = "$prog_name -write_ini $ini_file";
  system($cmd)==0  or  die "Error: $cmd\n";

  local($/) = undef;
  open(INF, $ini_file)  or die "Can't open $ini_file\n";
  my $text = <INF>;
  close(INF);
  unlink($ini_file);

  my ($key, $val);
  while (($key, $val) = each %values) {
    my ($node, $item) = ($key =~ /(\w+):(\w+)/);
    my $pattern = "<NODE name=\"$node\".*?<ITEM name=\"$item\" value=";
    $text =~ s/($pattern)"[^"]+"/\1"$val"/s;
  }
  my $out_file = tmpnam();
  open(OUTF, ">$out_file")  or  die "Can't create $out_file\n";
  print OUTF $text;
  close(OUTF);
  
  $out_file;
}


#---------------------------------------------------------------------------
sub handle_consensusElement {
  my $elem_node = shift;

  my $cons_id = $elem_node->[0]{id};
  my $centroid = $elem_node->[4][0];

  my $e_list = $elem_node->[8];
  my $num_elems = (scalar(@$e_list) - 2)/4;

  my @elem_ids = map { $e_list->[4*$_ + 4][0]{id} } (0..$num_elems-1);

  my $elem = {
    consensus_compound_id  => $elem_node->[0]{id},
    quality                => $elem_node->[0]{quality},
    mass                   => $centroid->{mz},
    rt                     => $centroid->{rt},
    quantity               => $centroid->{it},
    obs_mass_list          => [ @elem_ids ],
  };
  
  $elem;
}


#---------------------------------------------------------------------------
sub get_parser_root {
  my ($in_file) = @_;

  my $p = XML::Parser->new(Style => 'Tree');
  
  my $root;
  if ($in_file) {
    $root = $p->parsefile($in_file);
  }
  else {
    local($/);
    my $text = <>;
    $root = $p->parse($text);
  }

  $root;
}


#---------------------------------------------------------------------------
sub parse_consensus_file {
  my ($cons_file) = @_;

  my $root = get_parser_root($cons_file);

  die "XML is not as expected"
      unless $root->[0] eq 'consensusXML' and
      $root->[1][11] eq 'consensusElementList';

  my $cons_elems_list = $root->[1][12];
  my $num_cons_elems = (scalar(@{$cons_elems_list}) - 3) / 4;
  
  map {
    die "XML is not as expected"
        unless $cons_elems_list->[4*$_ + 3] eq 'consensusElement';
    handle_consensusElement($cons_elems_list->[4*$_ + 4]);
  }  (0..$num_cons_elems-1);

}


#---------------------------------------------------------------------------
#  Return hash with feature_id as key and adjusted retention time as value.
#---------------------------------------------------------------------------
sub parse_rt_from_featurexml_file {
  my ($feat_file) = @_;

  my $root = get_parser_root($feat_file);
  die "XML is not as expected"
      unless $root->[0] eq 'featureMap' and
             $root->[1][7] eq 'featureList';

  my $feature_root = $root->[1][8];
  my $num_features = $feature_root->[0]{count};

  my %h = map {
    my $f = $feature_root->[4*$_ + 4];
    ($f->[0]{id}, $f->[4][2]);

  } (0..$num_features-1);

}


#---------------------------------------------------------------------------
sub output_peptide_hit {
  my ($w, $row, $match_run_id) = shift;

  die "Peptide hits not implemented yet\n";

  $w->startTag('PeptideIdentification', 
               identification_run_ref=>$match_run_id,
               higher_score_better=>'true',
               significance_threshold=>1.3
               );
  $w->endTag();
}


#---------------------------------------------------------------------------
#  Given an array of rows with fields (mass, rt, quantity, id), 
#  create a ".featureXML" file.
#---------------------------------------------------------------------------
sub create_featurexml_file {
  my ($featuremap_name, $outfile_name, $a_rows) = @_;

  my $outf = new IO::File();
  $outf->open(">$outfile_name");
  my $w = XML::Writer->new(OUTPUT=>$outf, DATA_INDENT=>2);

  print $outf '<?xml version="1.0" encoding="ISO-8859-1"?>', "\n";

  $w->startTag('featureMap', id=>$featuremap_name);
    print $outf "\n";
    $w->startTag('featureList', count=>scalar(@$a_rows) );
    print $outf "\n  ";
    my $cnt = 1;
    foreach my $row (@$a_rows) {
      $w->startTag('feature', id=>$row->{id});
        print $outf "\n    ";
        $w->dataElement('position', $row->{rt}, dim=>0); 
  
        print $outf "\n    ";
        $w->dataElement('position', $row->{mass}, dim=>1); 
  
        print $outf "\n    ";
        $w->dataElement('intensity', sprintf("%d", $row->{quantity}) );
  
        print $outf "\n    ";
        $w->dataElement('charge', '0');
        print $outf "\n  ";

      $w->endTag(); # feature
      print $outf "\n  ";
    }
    $w->endTag();   # featureList
  print $outf "\n";
  $w->endTag();     # featureMap

  $w->end();

  $outf->close();
}


#------------------------------------------------------
sub dataset_adjust_retention_times {
  my ($dbh, $dataset, $mso_consensus_file) = @_;

  my $sth = $dbh->prepare(
      'SELECT mass, rt, quantity, id FROM observed_mass WHERE dataset = ?');
  $sth->execute($dataset);

  my ($row, @rows);
  while ($row = $sth->fetchrow_hashref()) {
    push @rows, $row;
  }

  my $tmp_features_file = tmpnam() . '.featureXML';
  #-------------------------------------------------------------------
  #  Create a ".featureXML" file for the ms/ms dataset.
  #  Determine adjusted retention times to align with the MS-only consensus.
  #-------------------------------------------------------------------
  openms_utils::create_featurexml_file($dataset, $tmp_features_file, \@rows);

  my $feat_file = tmpnam() . '.featureXML';
  my $cmd = "MapAlignerPoseClustering -in $tmp_features_file " . 
            "-out $feat_file " .
            "-reference:file  $mso_consensus_file";
  system($cmd)==0  or  die "Error: $cmd\n";
  unlink($tmp_features_file);

  #-----------------------------------------------------------------------
  #  Parse the modified features file.
  #  A hash of obs_mass_id => rt is returned.
  #  Since the obs_mass_id has "f_" prepended, strip that off.
  #-----------------------------------------------------------------------
  my @feature_data = openms_utils::parse_rt_from_featurexml_file($feat_file);
  s/^f_//  foreach @feature_data;
  my %feat_rt = @feature_data;

  my $tmp_copy_file = tmpnam();
  open(OUTF, ">$tmp_copy_file")  or die "Can't create file $tmp_copy_file\n";
  my ($id, $rt);
  while (($id, $rt) = each(%feat_rt)) {
    print OUTF "$id\t$rt\n";
  }
  close(OUTF);

  #-------------------------------------------------------------------
  #  Load modified retention times into a temporary table in the database.
  #-------------------------------------------------------------------
  $dbh->do('CREATE TEMP TABLE tmp_mod_rt(id integer, rt real);') or
      die "Can't create temp table tmp_mod_rt\n";

  $dbh->do("COPY tmp_mod_rt(id,rt) FROM '$tmp_copy_file';")  or
      die "Error copying data to tmp_mod_rt\n";

  unlink($tmp_copy_file);

  my $query = <<UPDATE_QUERY;
    UPDATE observed_mass o
    SET rt_adjustment_to_consensus = m.rt - o.rt
    FROM tmp_mod_rt m
    WHERE o.id = m.id
UPDATE_QUERY
  $dbh->do($query) or die "Error: $query\n";

  $dbh->do('DROP TABLE tmp_mod_rt')  or die "Error dropping tmp_mod_rt\n";
}

1;
