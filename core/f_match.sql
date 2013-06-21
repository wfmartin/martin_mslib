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

CREATE OR REPLACE FUNCTION combine_match_runs(
    p_match_run_ids          integer[])
    RETURNS integer
    LANGUAGE plpgsql AS $$
DECLARE
  v_new_match_run_id         integer;
  v_dataset                  varchar;
  v_num_matches              integer;
  v_num_decoy_matches        integer;
BEGIN
  INSERT INTO match_run(dataset, match_method, match_source,
                        error_ppm_threshold)
    SELECT
      dataset,
      'combine match_run rows',
      string_agg(match_run_id::varchar, ','  ORDER BY match_run_id)
          AS match_source,
      max(error_ppm_threshold)  AS error_ppm_threshold
    FROM match_run
    WHERE match_run_id = ANY($1)
    GROUP BY dataset
  RETURNING match_run_id, dataset INTO v_new_match_run_id, v_dataset;

  --------------------------------------------------------------------
  --  For any peak and compound, chose the higher score between the
  --  two methods.
  --------------------------------------------------------------------
  INSERT INTO match(match_run_id, obs_mass_id, compound_w_origin_id,
                    error_ppm, score, spectrum_id, is_decoy, misc)
    SELECT
      v_new_match_run_id AS match_run_id,
      obs_mass_id,
      compound_w_origin_id,
      (array_agg(error_ppm ORDER BY score DESC))[1]  AS error_ppm,
      max(score) AS score,
      spectrum_id,
      is_decoy,
      string_agg(misc, '+')  AS misc
    FROM match
    WHERE match_run_id = ANY($1)
    GROUP BY obs_mass_id, spectrum_id, compound_w_origin_id, is_decoy;

  SELECT
    sum( case when is_decoy then 0 else 1 end )::integer AS num_matches,
    sum( case when is_decoy then 1 else 0 end )::integer AS num_decoy_matches
  INTO v_num_matches, v_num_decoy_matches
  FROM match
  WHERE match_run_id = v_new_match_run_id;

  UPDATE match_run
  SET num_matches = v_num_matches,
      num_decoy_matches = v_num_decoy_matches
  WHERE match_run_id = v_new_match_run_id;

  RETURN v_new_match_run_id;
END;
$$;



--------------------------------------------------------------------------
--  When there's a significant mass error in a match, that probably means
--  that the computed mass of the compound found by the matching program
--  is probably correct rather than the mass determined by
--  "Find Compounds by Molecular Feature" (which could have errors in 
--  determining charge, number of neutrons, or modifications).
--
--  This function will attempt to correct for that scenario by splitting 
--  off the spectrum from the 'observed_mass' row, creating a new row
--  using the computed compound mass instead of an observed mass.
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION obs_mass_correct_from_matches(
    p_match_run_id           integer)
    RETURNS void
    LANGUAGE plpgsql AS $$
DECLARE
  v_dataset                  varchar;
  v_min_peak_id              integer;
  v_seq_start                integer;
BEGIN
  SELECT dataset INTO v_dataset
  FROM match_run
  WHERE match_run_id = p_match_run_id;

  CREATE TABLE tmp_bad_mass AS 
    SELECT
      m.id AS  match_id,
      m.obs_mass_id AS orig_obs_mass_id,
      nextval('obs_mass_clone_seq')::integer AS obs_mass_id,
      c.mass AS mass,
      s.component_charge AS min_z,
      s.component_charge AS max_z,
      s.rt_start,
      s.rt_end
    FROM match m
    JOIN compound_w_origin c USING(compound_w_origin_id)
    JOIN molf_spectrum s
         ON (s.dataset = v_dataset AND s.spectrum_id = m.spectrum_id)
    WHERE match_run_id = p_match_run_id AND ABS(error_ppm) > 30;

  INSERT INTO obs_mass_correction(new_obs_mass_id, orig_obs_mass_id)
    SELECT obs_mass_id, orig_obs_mass_id
    FROM tmp_bad_mass;

  UPDATE match m
  SET obs_mass_id = t.obs_mass_id,
      error_ppm = 0,
      misc = COALESCE(misc,'') || ' (corrected observed_mass)'
  FROM tmp_bad_mass t
  WHERE m.id = t.match_id;

  SELECT least(min(peak_id),0) - 1  INTO v_seq_start
  FROM observed_mass
  WHERE dataset = v_dataset;

  EXECUTE 'CREATE TEMP SEQUENCE tmp_seq INCREMENT -1  START ' || v_seq_start;

  INSERT INTO observed_mass(id, dataset, peak_id, mass, rt, rt_start, rt_end,
      min_z, max_z, rt_adjustment_to_consensus)
    SELECT 
      m.obs_mass_id AS id,
      v_dataset,
      nextval('tmp_seq')::integer AS peak_id,
      m.mass,
      (m.rt_start + m.rt_end)/2::real  AS rt,
      m.rt_start,
      m.rt_end,
      m.min_z,
      m.max_z,
      o.rt_adjustment_to_consensus
    FROM tmp_bad_mass m
    JOIN observed_mass o ON (m.orig_obs_mass_id = o.id);

  DROP SEQUENCE tmp_seq;
  DROP TABLE tmp_bad_mass;
END;
$$;

