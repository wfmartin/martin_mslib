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
CREATE TABLE molf_dataset(
  dataset                    varchar  PRIMARY KEY,
  lc_config_id               varchar  NOT NULL,
  csv_filename               varchar,
  cef_filename               varchar,
  mgf_filename               varchar,
  frag_data_is_raw           boolean  DEFAULT False
);


-------------------------------------------------------------------------
CREATE TYPE molf_peak AS (
  mass                       double precision,
  ion_count                  real
);


-----------------------------------------------------------------------
-- One row per spectrum in the CEF file
-- (<p> elements within <MsPeaks>
--
--  Note:  mz  comes from the MGF file (because of better precision)
--         when there is tandem data, otherwise from cef file.

--  The columns that come from the CSV file are a denormalized part of this
--  table and duplicated for each component spectrum.  Most notably, the
--  volume and volume_pct columns would be multiply counted if summed across
--  rows.  They are distinct by dataset+cpd, however.
-----------------------------------------------------------------------
CREATE TABLE molf_spectrum(
  dataset                    varchar,
  spectrum_id                integer, -- automatically generated (per spectrum)
  cpd                        integer,
  csv_mass                   double precision,
  component_charge           integer,   -- from CEF
  component_recipe           varchar,   -- from CEF
  min_z                      integer,   -- from CSV 
  max_z                      integer,   -- from CSV
  mz                         double precision,
  mgf_rt                     real,
  csv_rt                     real,
  rt_width_at_half_ht        real,
  rt_start                   real,
  rt_end                     real,
  num_neutrons               integer,
  volume                     real,
  volume_pct                 real,
  height                     real,
  peaks                      molf_peak[],
  collision_energy           real,
  mgf_index                  integer,  -- one based index into the MGF file
  mgf_title                  varchar,

  PRIMARY KEY(dataset, spectrum_id)
);
CREATE INDEX molf_spectrum_cpd_key ON molf_spectrum(dataset, cpd);
CREATE INDEX molf_title_key ON molf_spectrum(mgf_title);


-----------------------------------------------------------------------
--  Fragment peaks (integrated together) from Find Compound by Mol Feature.
-----------------------------------------------------------------------
CREATE TABLE molf_fragment_peak(
  dataset                    varchar,
  spectrum_id                integer,
  peak_index                 integer,
  mass                       double precision,
  ion_count                  real,

  PRIMARY KEY(dataset, spectrum_id, peak_index)
);
