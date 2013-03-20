--#########################################################################
-- Copyright (C) 2012 William F. Martin
--
-- This program is free software; you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by the
-- Free Software Foundation;
--
-- This program is distributed in the hope that it will be useful, but
-- WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
-- See the GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License along
-- with this program; if not, write to the Free Software Foundation, Inc.,
-- 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
--#########################################################################

-----------------------------------------------------------------------
--  Alterations:
--   If the CEF file showed neutrons, then subtract their mass from mz,
--   which is output to PEPMASS line and in the TITLE line.
--
--   Replace TITLE line with some content that is useful for us.
-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION agilent_export_mgf_for_gpm(
    p_dataset                varchar)
    RETURNS varchar AS $$

  my ($dataset) = @_;

  my $query = <<CHG_CARRIER;
    SELECT frag_data_is_raw, charge_carrier_mass
    FROM molf_dataset
    JOIN lc_configuration USING(lc_config_id)
    WHERE dataset = '$dataset'
CHG_CARRIER

  my $rv = spi_exec_query($query);
  my $charge_carrier_mass = $rv->{rows}[0]{charge_carrier_mass};
  my $frag_data_is_raw = $rv->{rows}[0]{frag_data_is_raw};
  my $charge_mass_mod = $frag_data_is_raw ? 0 : $charge_carrier_mass;

  my $query = <<QUERY;
    SELECT cpd, mgf_rt, mz, component_charge, spectrum_id,
           num_neutrons
    FROM molf_spectrum AS s
    WHERE dataset = '$dataset' AND peaks IS NOT NULL
    ORDER BY mgf_rt
QUERY
  my $rv = spi_exec_query($query);

  $query = <<FRAG_QUERY;
SELECT *
FROM molf_fragment_peak
WHERE dataset='$dataset' AND spectrum_id=\$1
FRAG_QUERY

  my $peak_sth = spi_prepare($query, 'INTEGER');

  my $fmt_str = <<FMT;
BEGIN IONS
PEPMASS=%.7f
CHARGE=%s
TITLE=dataset=%s spectrum_id=%d  cpd=%d  mz=%.6f  z=%d
RTINSECONDS=%.3f
%s
END IONS
FMT

  my $out_str = '';

  foreach my $spec_row (@{$rv->{rows}}) {
    my $backwards_charge = abs($spec_row->{component_charge}) .
        ($spec_row->{component_charge} > 0 ? '+' : '-');
    
    my $peak_rv = spi_exec_prepared($peak_sth, $spec_row->{spectrum_id});

    my $peak_str = join("\n", map {
        my $mz = $_->{mass} + $charge_mass_mod;
        "$mz\t$_->{ion_count}";
      } @{$peak_rv->{rows}} );

    my $adjusted_mz = $spec_row->{mz} -
        ($spec_row->{num_neutrons} * 1.008664916 /
          $spec_row->{component_charge});

    $out_str .= sprintf($fmt_str, 
        $adjusted_mz,
        $backwards_charge,
        $dataset, 
        $spec_row->{spectrum_id},
        $spec_row->{cpd},
        $adjusted_mz,
        $spec_row->{component_charge}, 
        60.0 * $spec_row->{mgf_rt}, 
        $peak_str);
  }

  $out_str;

$$ LANGUAGE plperl;
