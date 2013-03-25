CREATE OR REPLACE FUNCTION create_global_map(
    p_consensus_id           varchar,
    p_match_run_ids          integer[])
    RETURNS TABLE(
      match_run_id           integer,
      obs_mass_id            integer,
      consensus_compound_id  integer,
      has_match              boolean,
      is_fragmented          boolean
    )
    LANGUAGE plpgsql AS $$
DECLARE
  v_match_run_id             integer;
BEGIN
  
  FOREACH v_match_run_id IN ARRAY p_match_run_ids  LOOP

    RETURN QUERY
    WITH tmp_map AS (
      -- Cols: obs_mass_id, mass_rt_rectangle, match_ids[], consensus_cpd_id
      SELECT *
      FROM map_msms_dataset(p_consensus_id, v_match_run_id, 1.3::real)
      )
    SELECT
      v_match_run_id AS match_run_id,
      t.obs_mass_id,
      t.consensus_compound_id,

      (SELECT bool_or(NOT is_decoy) 
       FROM match WHERE id = ANY(match_ids)
      )  AS has_match,

      (SELECT bool_or(num_peaks > 100)
       FROM observed_mass o
       JOIN component_spectrum cs
         ON (cs.dataset = o.dataset AND
             cs.parent_peak_id = o.peak_id)
       WHERE o.id = t.obs_mass_id
       )  AS is_fragmented
    FROM tmp_map t
    WHERE t.consensus_compound_id IS NOT NULL;

  END LOOP;

END;
$$;


CREATE TABLE cpd_ident_stats AS
  WITH pass_1 AS (
    SELECT * FROM create_global_map('mso',
        ARRAY[3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57])
    ),
  pass_2 AS (
    SELECT
      consensus_compound_id,
      count(*) AS num_peaks,
      (sum( CASE WHEN is_fragmented or has_match THEN 1 ELSE 0 END ))::integer
          AS num_times_fragmented,
      (sum( CASE WHEN has_match THEN 1 ELSE 0 END ))::integer
          AS num_times_matched
    FROM pass_1
    GROUP BY consensus_compound_id
    )
  SELECT * FROM pass_2
  WHERE num_times_matched > 0;
