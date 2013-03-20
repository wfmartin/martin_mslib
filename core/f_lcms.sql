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


----------------------------------------------------------------------------
--  CREATE a temp table  tmp_cpd  with compounds, keeping FDR under limit.
--
--  Iterate through the compounds in the consensus in order of decreasing score 
--  returning the putative (non-decoy) compounds until the portion of decoys
--  exceeds the given maximum for false discovery rate (p_max_fdr).
----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION get_cpds_from_consensus(
    p_consensus_id           varchar,
    p_max_fdr                real)
    RETURNS real 
    LANGUAGE plpgsql AS $$
DECLARE
  v_cpd                      record;
  v_num_cpds                 integer;
  v_num_decoy_cpds           integer;
  v_new_fdr                  real;
  v_fdr                      real;
  v_lcms_pseudo_score        real;
BEGIN
  v_lcms_pseudo_score := 1e30;
  v_num_cpds := 0;
  v_num_decoy_cpds := 0;

  CREATE TEMP TABLE tmp_cpd(
    consensus_compound_id    integer,
    compound_w_origin_id     integer,
    rt_start                 real,
    rt_end                   real,
    min_z                    integer,
    max_z                    integer,
    mass                     double precision,
    mass_rt_rectangle        box
  );

  FOR v_cpd IN
    WITH tmp_map_1 AS (
      SELECT
        consensus_compound_id,
        min_z,
        max_z,
        False AS is_decoy,
        unnest(compound_w_origin_ids) AS compound_w_origin_id,
        v_lcms_pseudo_score  AS score
      FROM consensus_cpd_lcms_map m
      JOIN lcms_library_compound c USING(lcms_library_compound_id)
      WHERE consensus_id = p_consensus_id
        UNION
      SELECT
        consensus_compound_id,
        min_z,
        max_z,
        is_decoy,
        compound_w_origin_id,
        score
      FROM consensus_cpd_map
      WHERE consensus_id = p_consensus_id  AND
            compound_w_origin_id IS NOT NULL
      ),
    tmp_map_2 AS (
      SELECT
        consensus_compound_id,
        min(min_z) AS min_z,
        max(max_z) AS max_z,
        is_decoy,
        compound_w_origin_id,
        max(score) AS score
      FROM tmp_map_1
      GROUP BY consensus_compound_id, is_decoy, compound_w_origin_id
      )
    SELECT
      consensus_compound_id,
      compound_w_origin_id,
      mass,
      score,
      is_decoy,
      rt_start,
      rt_end,
      m.min_z,
      m.max_z,
      is_decoy,
      mass_rt_rectangle
    FROM tmp_map_2 m
    JOIN consensus_compound USING(consensus_compound_id)
    WHERE consensus_id = p_consensus_id
    ORDER BY score DESC
  LOOP
    IF v_cpd.is_decoy THEN
      v_num_decoy_cpds := v_num_decoy_cpds + 1;
      EXIT WHEN p_max_fdr IS NOT NULL AND v_num_decoy_cpds >= 2 AND
        v_num_decoy_cpds::real/(v_num_cpds + v_num_decoy_cpds)::real >
            p_max_fdr;
    ELSE
      --------------------------------------------------------------------
      --  Exclude matches to the LCMS library from the FDR computation.
      --------------------------------------------------------------------
      IF (v_cpd.score < v_lcms_pseudo_score) THEN
        v_num_cpds := v_num_cpds + 1;
      END IF;

      INSERT INTO tmp_cpd(consensus_compound_id, compound_w_origin_id,
              rt_start, rt_end, min_z, max_z, mass, mass_rt_rectangle)
      VALUES(v_cpd.consensus_compound_id, v_cpd.compound_w_origin_id,
          v_cpd.rt_start, v_cpd.rt_end, v_cpd.min_z, v_cpd.max_z, v_cpd.mass,
          v_cpd.mass_rt_rectangle);
    END IF;

    CONTINUE WHEN v_num_cpds = 0;

    v_new_fdr := v_num_decoy_cpds::real/(v_num_cpds + v_num_decoy_cpds)::real;
    EXIT WHEN p_max_fdr IS NOT NULL AND v_num_decoy_cpds >= 2 AND
        v_new_fdr > p_max_fdr;

    v_fdr := v_new_fdr;
  END LOOP;

  RETURN v_fdr;
END;
$$;


----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION lcms_lib_create_from_consensus(
    IN p_lcms_library_id     varchar,
    IN p_consensus_id        varchar,
    IN p_max_fdr             real,
    OUT num_cpds             integer,
    OUT fdr                  real)
    LANGUAGE plpgsql AS $$
BEGIN
  INSERT INTO lcms_library(lcms_library_id, bio_context_id, lc_config_id,
      consensus_ids)
    SELECT
      p_lcms_library_id,
      bio_context_id,
      lc_config_id,
      ARRAY[consensus_id]::varchar[] AS consensus_ids
    FROM sample_consensus
    JOIN consensus_parameters USING(consensus_parameters_key)
    WHERE consensus_id = p_consensus_id;

  fdr := get_cpds_from_consensus(p_consensus_id, p_max_fdr);

  INSERT INTO lcms_library_compound(lcms_library_id, compound_w_origin_ids,
      rt_start, rt_end, min_z, max_z, mass, mass_rt_rectangle)
    SELECT
      p_lcms_library_id,
      array_agg(compound_w_origin_id) AS compound_w_origin_ids,
      rt_start, rt_end, min_z, max_z, mass,
      (array_agg(mass_rt_rectangle))[1]
    FROM tmp_cpd
    GROUP BY consensus_compound_id, rt_start, rt_end, min_z, max_z, mass;

  GET DIAGNOSTICS num_cpds = ROW_COUNT;

  DROP TABLE tmp_cpd;
END;
$$;


----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION validate_lib_consensus_params(
    p_lcms_library_id        varchar,
    p_consensus_id           varchar)
    RETURNS void 
    LANGUAGE plpgsql AS $$
DECLARE
  v_lcms_lib                 lcms_library;
  v_bio_context_id           varchar;
  v_lc_config_id             varchar;
BEGIN
  SELECT * INTO v_lcms_lib
  FROM lcms_library
  WHERE lcms_library_id = p_lcms_library_id;

  -------------------------------------------------------------------------
  -- Validate parameters:
  --  1) Specified LCMS library exists.
  --  2) Specified sample_consensus exists.
  --  3) Bio-context and lc_config_id match between LCMS library and
  --     sample_consensus.
  -------------------------------------------------------------------------
  IF NOT FOUND THEN
    RAISE EXCEPTION 'LCMS library (%) not found.', p_lcms_library_id;
  END IF;

  IF p_consensus_id = ANY (v_lcms_lib.consensus_ids) THEN
    RAISE EXCEPTION 'LCMS library (%) already contains consensus (%).',
        p_lcms_library_id, p_consensus_id;
  END IF;

  SELECT bio_context_id, lc_config_id
  INTO v_bio_context_id, v_lc_config_id
  FROM sample_consensus
  JOIN consensus_parameters USING(consensus_parameters_key)
  WHERE consensus_id = p_consensus_id;

  IF NOT FOUND THEN
    RAISE EXCEPTION 'Cum matches group (%) not found.', p_consensus_id;
  END IF;

  IF v_lcms_lib.bio_context_id <> v_bio_context_id  THEN
    RAISE EXCEPTION 'The values for bio_context_id do not match between '
                'lcms_library (%) and sample_consensus (%).',
                v_lcms_lib.bio_context_id, v_bio_context_id;
  END IF;

  IF v_lcms_lib.lc_config_id <> v_lc_config_id  THEN
    RAISE EXCEPTION 'The values for lc_config_id do not match between '
                'lcms_library (%) and sample_consensus (%).',
                v_lcms_lib.lc_config_id, v_lc_config_id;
  END IF;
END;
$$;
 

----------------------------------------------------------------------------
--  Add cumulative matches to an MS library of compounds
----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION lcms_lib_add_consensus(
    p_lcms_library_id        varchar,
    p_consensus_id           varchar,
    p_max_fdr                real)
    RETURNS integer
    LANGUAGE plpgsql AS $$
DECLARE
  v_num_cpds                 integer;
BEGIN
  PERFORM validate_lib_consensus_params(p_lcms_library_id, p_consensus_id);

  --------------------------------------------------------------------------
  UPDATE lcms_library
  SET consensus_ids = consensus_ids || p_consensus_id
  WHERE lcms_library_id = p_lcms_library_id;

  --  Create tmp_cpd table :
  PERFORM get_cpds_from_consensus(p_consensus_id, p_max_fdr);

  --------------------------------------------------------------------------
  --  Determine overlaps of new cumulative match group compounds with
  --  those already in the library.
  --------------------------------------------------------------------------
  CREATE TEMP TABLE tmp_compound_map AS
    WITH tmp_consensus_cpd AS (
      SELECT
        array_agg(compound_w_origin_id) AS compound_w_origin_ids,
        rt_start, rt_end, min_z, max_z, mass,
        (array_agg(mass_rt_rectangle))[1] AS mass_rt_rectangle
      FROM tmp_cpd
      GROUP BY consensus_compound_id, rt_start, rt_end, min_z, max_z, mass
      ),
    tmp_join AS (
      SELECT cc.consensus_compound_id, lc.lcms_library_compound_id
      FROM tmp_consensus_cpd cc, lcms_library_compound lc
      WHERE lcms_library_id = p_lcms_library_id  AND
            cc.mass_rt_rectangle && lc.mass_rt_rectangle
      )
    SELECT
      lc.id AS lcms_compound_id,
      cmc.consensus_compound_id,
      array_unique_vals(lc.compound_w_origin_ids || cmc.compound_w_origin_ids)
          AS compound_w_origin_ids,
      LEAST(rt_start, rt_start)      AS rt_start,
      GREATEST(rt_end, rt_end)       AS rt_end,
      LEAST(cmc.min_z, lc.min_z)     AS min_z,
      GREATEST(cmc.max_z, lc.max_z)  AS max_z,
      cmc.mass,
      cmc.mass_rt_rectangle,
      CASE
        WHEN lc.id IS NULL THEN 1   -- new compound
        WHEN rt_start BETWEEN rt_start and rt_end  AND
             rt_end   BETWEEN rt_start and rt_end
            THEN 0                  -- no change needed
        ELSE 2                      -- update
      END
          AS update_mode
    FROM tmp_consensus_cpd cmc
    LEFT OUTER JOIN tmp_join USING(consensus_compound_id)
    LEFT OUTER JOIN lcms_library_compound USING(lcms_library_id);

--  LEFT OUTER JOIN lcms_library_compound lc ON 
--   (lc.lcms_library_id = p_lcms_library_id  AND
--    lc.compound_w_origin_ids && cmc.compound_w_origin_ids  AND
--    ranges_overlap(rt_start, rt_end, rt_start, rt_end));

  DROP TABLE tmp_cmd;

  --------------------------------------------------------------------------
  --  Compounds not already in the library:
  --------------------------------------------------------------------------
  INSERT INTO lcms_library_compound(lcms_library_id, compound_w_origin_ids,
                           rt_start, rt_end, min_z, max_z, mass)
    SELECT
      p_lcms_library_id,
      compound_w_origin_ids,
      rt_start,
      rt_end,
      min_z,
      max_z,
      mass
    FROM tmp_compound_map
    WHERE update_mode = 1;

  GET DIAGNOSTICS v_num_cpds = ROW_COUNT;

  --------------------------------------------------------------------------
  --  If a consensus_compound overlaps multiple lcms_library_compounds, then they
  --  need to be merged.
  --  Keep the lowest id number and delete the others.
  --------------------------------------------------------------------------
  CREATE TEMP TABLE tmp_merged  AS
    SELECT
      min(lcms_compound_id)        AS lcms_compound_id,
      min(rt_start)                AS rt_start,
      max(rt_end)                  AS rt_end,
      min(min_z)                   AS min_z,
      max(max_z)                   AS max_z,
      array_agg(lcms_compound_id)  AS all_ids
    FROM tmp_compound_map
    WHERE update_mode = 2
    GROUP BY consensus_compound_id;

  DROP TABLE tmp_compound_map;

  UPDATE lcms_library_compound lc
  SET
    rt_start = m.rt_start,
    rt_end   = m.rt_end,
    min_z    = m.min_z,
    max_z    = m.max_z
  FROM tmp_merged m
  WHERE lc.id = m.lcms_compound_id;

  WITH pass_1 AS (
    SELECT
      lcms_compound_id,
      unnest(all_ids) AS delete_id
    FROM tmp_merged
    ),
  delete_list AS (
    SELECT delete_id
    FROM pass_1
    WHERE delete_id <> lcms_compound_id
  )
  DELETE FROM lcms_library_compound
  USING delete_list
  WHERE id = delete_id;

  RETURN v_num_cpds;
END;
$$;


-----------------------------------------------------------------------
--
-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION match_against_lcms_lib(
    p_lcms_library_id        varchar,
    p_dataset                varchar,
    p_error_ppm_limit        real)
    RETURNS integer
    LANGUAGE plpgsql AS $$
DECLARE
  v_match_run_id             integer;
  v_num_matches              integer;
BEGIN
  INSERT INTO match_run(dataset, match_method, match_source,
                        error_ppm_threshold)
  VALUES(p_dataset, 'lcms_library', p_lcms_library_id, p_error_ppm_limit)
  RETURNING match_run_id INTO v_match_run_id;

  INSERT INTO match(match_run_id, obs_mass_id, compound_w_origin_id, error_ppm)
    SELECT
      v_match_run_id  AS match_run_id,
      o.id AS obs_mass_id,
      unnest(compound_w_origin_ids) AS compound_w_origin_id,
      ppm_error(c.mass, o.mass) AS error_ppm
    FROM lcms_library_compound c
    JOIN observed_mass o ON (
          ranges_overlap(o.rt_start, o.rt_end, c.rt_start, c.rt_end)  AND
          masses_within_error_limit(c.mass, o.mass, p_error_ppm_limit) )
    WHERE c.lcms_library_id = p_lcms_library_id  AND
          o.dataset = p_dataset;

  GET DIAGNOSTICS v_num_matches = ROW_COUNT;

  UPDATE match_run
  SET num_matches = v_num_matches
  WHERE match_run_id = v_match_run_id;

  RETURN v_match_run_id;
END;
$$;
