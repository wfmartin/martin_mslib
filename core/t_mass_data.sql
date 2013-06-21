--#########################################################################
-- Copyright (C) 2013 William F. Martin
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

--*************************************************************************
--**  Observed MASSES data  (tables)
--*************************************************************************

-------------------------------------------------------------------------
--  Note that these values are not masses, but rather mz values.
-------------------------------------------------------------------------
CREATE TABLE calibrant_ion(
  calibrant_ion_list_id      varchar,
  mz                         double precision,
  charge                     integer,
  mass                       double precision
);


-------------------------------------------------------------------------
--  An entity to encapsulate parameters in the liquid chromotography that
--  would affect retention times (for a compound library).
--  The obvious example is the solvent gradient.
--
--  This table could be filled out with multiple columns for structured
--  information, but now it's just free-form text in the description column.
-------------------------------------------------------------------------
CREATE TABLE lc_configuration(
  lc_config_id               varchar PRIMARY KEY,
  description                varchar,
  theoretical_max_mz         real,  -- limits of instrument's sensitivity
  theoretical_max_rt         real,  -- time at end of rt run
  calibrant_ion_list_id      varchar,
  calibrant_error_ppm_limit  real,
  charge_carrier_mass        double precision
);


-------------------------------------------------------------------------
--  This table indicates the kind of sample (i.e., human milk) was input
--  to the LCMS instrument.
--  
--  A compound library makes most sense relative to the kind of sample
--  analyzed.
-------------------------------------------------------------------------
CREATE TABLE bio_context(
  bio_context_id             varchar PRIMARY KEY,
  description                varchar
  --  proteins_list_id  (add later)
);


-------------------------------------------------------------------------
--  A summary record of a dataset.
-------------------------------------------------------------------------
CREATE TABLE mass_dataset(
  dataset                    varchar PRIMARY KEY,

  bio_context_id             varchar,
  lc_config_id               varchar,
  data_source_type           varchar,

  num_peaks                  integer,
  frag_data_is_raw           boolean  DEFAULT False,
  rt_adjustment_ts           timestamp,

  consensus_id               varchar,
  merge_source_dataset       varchar,

  create_ts                  timestamp  DEFAULT now()
);


-------------------------------------------------------------------------
--  MS (not fragment) peak.
-------------------------------------------------------------------------
CREATE TABLE observed_mass (
  id                         SERIAL PRIMARY KEY,
  dataset                    varchar,
  peak_id                    integer,

  mass                       double precision,
  rt                         real,
  aligned_rt                 real,
  rt_width_at_half_ht        real,
  rt_start                   real,
  rt_end                     real,
  quantity                   real,
  rel_quantity               real,
  min_z                      integer,
  max_z                      integer,
  dominant_z                 integer,

  merged_ids                 integer[],
  rt_adjustment_to_consensus real,

  UNIQUE(dataset, peak_id)
);
CREATE INDEX observed_by_mass ON observed_mass(dataset, mass);

CREATE SEQUENCE obs_mass_clone_seq INCREMENT -1  MAXVALUE -1 OWNED BY observed_mass.id;


-------------------------------------------------------------------------
CREATE TABLE obs_mass_correction (
  new_obs_mass_id            integer PRIMARY KEY,
  orig_obs_mass_id           integer
);


-------------------------------------------------------------------------
--  MS peaks summed together to make 'observed_mass',
--  possibly with corresponding fragment spectra.
------------------------------------------------------------------------
CREATE TABLE component_spectrum(
  id                         SERIAL PRIMARY KEY,
  dataset                    varchar,
  spectrum_id                integer,

  parent_peak_id             integer, 
  charge                     integer,
  num_peaks                  integer,
  num_neutrons               integer  DEFAULT 0,
  collision_energy           real,

  UNIQUE(dataset, spectrum_id),
  FOREIGN KEY(dataset, parent_peak_id)
      REFERENCES observed_mass(dataset, peak_id)
);
CREATE INDEX fs_by_peak_id ON component_spectrum(dataset, parent_peak_id);


-------------------------------------------------------------------------
--  A peak in a fragment spectrum.
-------------------------------------------------------------------------
CREATE TABLE fragment_mass (
  id                         SERIAL PRIMARY KEY,
  dataset                    varchar,
  spectrum_id                integer,
  fragment_peak_id           integer,
  mass                       double precision,

  quantity                   real,

  UNIQUE(dataset, fragment_peak_id),
  FOREIGN KEY(dataset, spectrum_id) 
      REFERENCES component_spectrum(dataset, spectrum_id)
);


-------------------------------------------------------------------------
CREATE VIEW fragment_data AS
  SELECT
    num_peaks AS num_fragments,
    collision_energy,
    o.mass
  FROM component_spectrum s 
  JOIN observed_mass o  ON 
       (o.dataset = s.dataset AND o.peak_id = s.parent_peak_id)
  WHERE collision_energy > 0 AND num_peaks > 0;

