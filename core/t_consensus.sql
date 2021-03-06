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
-------------------------------------------------------------------------
-- Parameters for managing a sample_consensus that are likely to be shared
-- for multiple samples.
-------------------------------------------------------------------------
CREATE TABLE consensus_parameters(
  consensus_parameters_key   varchar PRIMARY KEY,

  lc_config_id               varchar,

  -- Parameters for FeatureLinkerUnlabeledQT (for mso consensus):
  mass_ppm_error_limit       real,
  cons_max_rt_diff           real,
  min_featurelinker_quality  real,
  min_featurelinker_num_items integer,
  min_quantity               real,

  -- Times to attempt identifying fragmented cpd before putting in excl list
  max_ident_attempts         integer,

  -- Times on the preferred list before excluding.
  max_times_preferred        integer,

  -- Maximum number of items in the exclusion/preferred list.
  max_excl_pref_list_length  integer,
  max_exclusion_iterations   integer,
  
  -- Minimum width (mass) and height (rt) for mass_rt_rectangle
  rect_mass_ppm_width        real,
  rect_rt_min_width          real,

  min_excl_pref_rt_width     real
);


-------------------------------------------------------------------------
CREATE TABLE cons_matching_params(
  cons_matching_params_key   varchar PRIMARY KEY,

  min_score                  real,  --   = -log(pvalue)

  proteins_fasta_file        varchar,  -- name of input file
  decoy_fasta_file           varchar,  -- NULL if no decoys
  decoy_is_random            boolean DEFAULT True,

  gpm_parameters_key         varchar,
  
  gfdb_mass_ppm_tolerance    real,
  gfdb_use_phosphorylation   boolean,
  gfdb_other_options         varchar
);


-------------------------------------------------------------------------
CREATE TABLE cons_matching_protein(
  cons_matching_params_key   varchar,
  protein_name               varchar,
  decoy                      boolean,
  aa_seq                     varchar,

  UNIQUE(cons_matching_params_key, protein_name, decoy)
);


-------------------------------------------------------------------------
--
-------------------------------------------------------------------------
CREATE TABLE sample_consensus(
  consensus_id               varchar PRIMARY KEY,

  consensus_parameters_key   varchar,
  cons_matching_params_key   varchar,

  bio_context_id             varchar,

  dataset_list               varchar[],
  create_ts                  timestamp  DEFAULT now(),

  list_produced              boolean    DEFAULT False,
  list_iteration             integer    DEFAULT 0,

  lcms_library_id            varchar
);
  

-------------------------------------------------------------------------
--
-------------------------------------------------------------------------
CREATE TABLE consensus_compound (
  consensus_compound_id      SERIAL PRIMARY KEY,
  consensus_id               varchar,
  quality                    real,
  mass                       double precision,
  rt                         real,
  rt_width_at_half_ht        real,  -- max of orig obs_mass
  rt_start                   real,  -- min of orig obs_mass
  rt_end                     real,  -- max of orig obs_mass
  quantity                   real,
  rel_quantity               real,

  ---------------------------------------------------------------------------
  -- The box type encapsulates a range of both mass (x-axis) and
  -- retention time (y-axis).  There are convenient and efficient operators
  -- for checking overlaps.
  ---------------------------------------------------------------------------
  mass_rt_rectangle          box,

  ---------------------------------------------------------------------------
  -- Columns for which MS-only observed masses comprise/create the consensus.
  -- (or in a few cases, MS/MS-only peaks when a GPM match peak didn't
  --  correspond to any existing compound consensus)
  --
  -- This doesn't include the MS/MS peaks which identify a compound but are
  -- not used to create the consensus.
  ---------------------------------------------------------------------------
  obs_mass_list              integer[],
  dataset_list               varchar[],
  msms_only                  boolean,

  min_z                      integer,
  max_z                      integer,
  dominant_z                 integer,

  ----------------------------------------------------------------------
  -- For optimization, these status values are kept here.
  -- When an identification has already been made or when the maximum
  -- number of attempts at fragmentation or identification have been made,
  -- then this compound is kept on the exclusion list.
  -- Once the 'exclude_compound' is set to true, the 3 columns above it
  -- are no longer updated.
  ----------------------------------------------------------------------
  has_non_decoy_ident        boolean  DEFAULT False,
  num_ident_attempts         integer  DEFAULT 0,
  num_times_preferred        integer  DEFAULT 0,
  exclude_compound           boolean  DEFAULT False,

  corrected_mass             double precision
);


-------------------------------------------------------------------------
CREATE TABLE consensus_cpd_merge_pair(
  consensus_id               varchar,
  receiver_cons_cpd_id       integer,
  donor_cons_cpd_id          integer
);

CREATE INDEX cc_merge_pair_ind ON consensus_cpd_merge_pair(consensus_id);


-------------------------------------------------------------------------
--  The output from FDR filtering.
-------------------------------------------------------------------------
CREATE TABLE final_consensus_compound (
  f_cons_cpd_id              SERIAL PRIMARY KEY,
  consensus_id               varchar,
  consensus_compound_id      integer,
  compound_w_origin_id       integer,
  rt_start                   real,
  rt_end                     real,
  min_z                      integer,
  max_z                      integer,
  quantity                   real,
  rel_quantity               real,
  score                      real
);


-------------------------------------------------------------------------
CREATE TABLE final_cons_cpd_cluster (
  cluster_id                 integer PRIMARY KEY,
  consensus_id               varchar,
  protein_name               varchar,
  start_pos                  integer,
  end_pos                    integer,
  seq_len                    integer,

  cpd_ids                    integer[]
);


-------------------------------------------------------------------------
CREATE TABLE cons_cpd_mods_map(
  consensus_id               varchar,
  base_cons_cpd_id           integer,
  mods_cons_cpd_id           integer,
  base_cpd_w_origin_id       integer,
  mods_cpd_w_origin_id       integer,
  protein_name               varchar,
  start_pos_on_protein       integer,
  seq_len                    integer,
  aa_seq                     varchar,
  ng_mods_str                varchar
);
 

-------------------------------------------------------------------------
CREATE TABLE consensus_cpd_lcms_map(
  consensus_compound_id      integer,
  consensus_id               varchar,

  lcms_library_compound_id   integer
);


-------------------------------------------------------------------------
--  A map of the MS/MS observed_mass rows to consensus_compound_rows,
--  along with matches (if any).
-------------------------------------------------------------------------
CREATE TABLE consensus_cpd_map(
  consensus_compound_id      integer,
  consensus_id               varchar,
  obs_mass_id                integer,
  min_z                      integer,
  max_z                      integer,

  is_fragmented              boolean,

  match_id                   integer,
  is_decoy                   boolean,
  compound_w_origin_id       integer,
  score                      real
);

CREATE INDEX cpd_cons_map_index ON consensus_cpd_map(consensus_compound_id);


-------------------------------------------------------------------------
CREATE TABLE consensus_removed_match(
  consensus_compound_id      integer,
  consensus_id               varchar,
  match_id                   integer
);


-------------------------------------------------------------------------
--  When identified compounds associated with the consensus_compound are
--  added, set the corrected_mass to the theoretical mass
--  (if not already set, and if the compound_w_origin mass clearly 
--   represents the same exact compound, including mods).
-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION corrected_mass_updater()
    RETURNS TRIGGER
    LANGUAGE plpgsql AS $$
DECLARE
  v_cons_cpd                 consensus_compound;
  v_cpd_mass                 numeric(10,5);
  v_max_ident_attempts       integer;
BEGIN
  SELECT *
  INTO v_cons_cpd
  FROM consensus_compound
  WHERE consensus_compound_id = NEW.consensus_compound_id;

  -------------------------------------------------------------------------
  --  If there's a new GPM identification:
  -------------------------------------------------------------------------
  IF NEW.match_id IS NOT NULL AND (NOT NEW.is_decoy) THEN
    v_cons_cpd.has_non_decoy_ident := True;
    v_cons_cpd.exclude_compound  := True;

    IF v_cons_cpd.corrected_mass IS NULL THEN
      SELECT mass INTO v_cpd_mass
      FROM compound_w_origin
      WHERE compound_w_origin_id = NEW.compound_w_origin_id;

      IF masses_within_error_limit(v_cpd_mass::double precision,
                                   v_cons_cpd.mass, 20::real)  THEN
        v_cons_cpd.corrected_mass := v_cpd_mass::double precision;
      END IF;
    END IF;

  ELSIF NEW.is_fragmented THEN
    v_cons_cpd.num_ident_attempts := v_cons_cpd.num_ident_attempts + 1;
  END IF;

  IF v_cons_cpd.min_z IS NULL OR NEW.min_z < v_cons_cpd.min_z
  THEN
    v_cons_cpd.min_z := NEW.min_z;
  END IF;

  IF v_cons_cpd.max_z IS NULL OR NEW.max_z > v_cons_cpd.max_z
  THEN
    v_cons_cpd.max_z := NEW.max_z;
  END IF;

  IF NOT v_cons_cpd.exclude_compound THEN
    SELECT max_ident_attempts INTO v_max_ident_attempts
    FROM sample_consensus
    JOIN consensus_parameters USING(consensus_parameters_key)
    WHERE consensus_id = NEW.consensus_id;

    IF v_cons_cpd.num_ident_attempts >= v_max_ident_attempts  THEN
      v_cons_cpd.exclude_compound := True;
    END IF;
  END IF;

  UPDATE consensus_compound
  SET min_z                =  v_cons_cpd.min_z,
      max_z                =  v_cons_cpd.max_z,
      corrected_mass       =  v_cons_cpd.corrected_mass,
      has_non_decoy_ident  =  v_cons_cpd.has_non_decoy_ident,
      num_ident_attempts   =  v_cons_cpd.num_ident_attempts,
      exclude_compound     =  v_cons_cpd.exclude_compound
  WHERE consensus_compound_id = NEW.consensus_compound_id;

  RETURN NEW;
END;
$$;

CREATE TRIGGER corrected_mass_updater_trigger
  AFTER INSERT on consensus_cpd_map
  FOR EACH ROW EXECUTE PROCEDURE corrected_mass_updater();


-- -------------------------------------------------------------------------
-- CREATE TABLE obs_mass_to_consensus_cpd(
--   obs_mass_id                integer,
--   consensus_id               varchar,
--   consensus_compound_id      integer
-- );


-------------------------------------------------------------------------
CREATE TABLE obs_mass_no_cons(
  obs_mass_id                integer,
  consensus_id               varchar,
  mass                       double precision,
  rt_start                   real,
  rt_end                     real,
  mass_rt_rectangle          box,
  min_z                      integer,
  max_z                      integer,
  dominant_z                 integer
);

CREATE INDEX obs_mass_no_cons_by_mass ON obs_mass_no_cons(mass);
--  To improve efficiency of retrieving in ascending mass order.


-------------------------------------------------------------------------
--  For observed_mass rows that didn't correspond to consensus_compound rows
--  when first encountered, but overlap later observed_mass rows that have
--  identified compound matches (and thus are added as consensus_compounds).
--
--  This corresponds to rows that are moved from 'obs_mass_no_cons' to
--  'consensus_cpd_map'.
-------------------------------------------------------------------------
CREATE TABLE obs_mass_promoted(
  obs_mass_id                integer,
  match_run_id               integer,   -- run causing promotion
  consensus_compound_id      integer
);

