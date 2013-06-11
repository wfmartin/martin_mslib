package gpm;

use strict;
#use POSIX;
use Cwd;
use DBI;
use DBD::Pg;
use File::Spec;


#-------------------------------------------------------------------------
sub create_taxonomy_file {
  my ($fasta_file) = @_;

  my $taxonomy_text = <<TAXONOMY;
<?xml version="1.0"?>
<bioml label="x! taxon-to-file matching list">
        <taxon label="protein_lib">
                <file format="peptide" URL="$fasta_file"/>
        </taxon>
</bioml>
TAXONOMY
  open(TAX_F, '>taxonomy.xml')  or die "Can't create taxonomy.xml\n";
  print TAX_F $taxonomy_text;
  close(TAX_F);
}


#-------------------------------------------------------------------------
#
#-------------------------------------------------------------------------
sub create_gpm_input_file {
  my ($dbh, $mgf_input_file, $gpm_parameters_key, $dataset) = @_;


  my $sth = $dbh->prepare('SELECT gpm_gen_input_xml(?,?,?)');
  $sth->execute($gpm_parameters_key, $dataset, $mgf_input_file);
  my ($input_text) = $sth->fetchrow_array();
  open(GPM_IN, '>input.xml')  or  die "Can't create input.xml\n";
  print GPM_IN $input_text;
  close(GPM_IN);
}


#-------------------------------------------------------------------------
sub find_tandem_exe {
  my $loc = `which tandem.exe`;
  return 'tandem.exe'  if $loc;

  $loc = `which tandem`;
  return 'tandem'  if $loc;

  die "Can't find tandem.exe tandem\n";
}


#-------------------------------------------------------------------------
#  Run GPM (global peptide machine) and return the name of the output file.
#-------------------------------------------------------------------------
sub run_gpm {
  my $ex_file = find_tandem_exe();
  my $cmd = "$ex_file input.xml > gpm_stdout";
  system($cmd)==0  or  die "Error: $cmd\n";

  my $dir = cwd();
  opendir(INDIR, $dir) or die "Can't read directory: $dir";

  my ($out_file, @ignore) =
      sort { -M "$dir/$a" <=> -M "$dir/$b" }
      grep { /^output\..*\.xml$/ } readdir(INDIR);
  closedir(INDIR);

  File::Spec->rel2abs($out_file);
}


#-------------------------------------------------------------------------
#
#-------------------------------------------------------------------------
sub get_gpm_matches {
  my ($dbh, $gpm_parameters_key, $dataset, $mgf_input_file,
      $proteins_fasta, $decoys_fasta, $min_score,
      $preserve_gpm_output) = @_;

  create_taxonomy_file($proteins_fasta);
  create_gpm_input_file($dbh, $mgf_input_file, $gpm_parameters_key, $dataset);

  my $gpm_outf = run_gpm();

  my $gpm_decoys_outf;
  if ($decoys_fasta) {
    create_taxonomy_file($decoys_fasta);
    $gpm_decoys_outf = run_gpm();
  }

  my $sth = $dbh->prepare('SELECT process_gpm_output(?,?,?,?,?)');
  $sth->execute($gpm_parameters_key, $dataset, $gpm_outf, $gpm_decoys_outf,
                $min_score)
      or die "Couldn't execute statement: " . $sth->errstr;

  my ($run_id) = $sth->fetchrow_array();

  unless ($preserve_gpm_output) {
    unlink($gpm_outf);
    unlink($gpm_decoys_outf)  if $decoys_fasta;
  }

  $run_id;
}

1;
