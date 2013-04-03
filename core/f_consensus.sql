-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION closest_consensus_cpd(
    p_obs_mass_id            integer,
    p_matched_cpd_ids        integer[],
    p_consensus_compound_ids integer[])
    RETURNS integer
    LANGUAGE SQL AS $$

    WITH om AS (  -- Fetch the observed_mass row for comparison
      SELECT mass, rt_start, rt_end
      FROM observed_mass
      WHERE id = $1
      ),

    tmp_cc_ids(consensus_compound_id) AS (SELECT unnest($3)),

    detailed AS (
      SELECT
        cc.consensus_compound_id,

        (EXISTS (SELECT *
         FROM consensus_cpd_map ccm
	 WHERE ccm.consensus_compound_id = cc.consensus_compound_id AND
               compound_w_origin_id = ANY($2) )
        )  AS has_cpds_in_common,

        ABS( ppm_error(COALESCE(corrected_mass,cc.mass), om.mass) )
            AS abs_ppm_error,

        gap_between_ranges(cc.rt_start, cc.rt_end, om.rt_start, om.rt_end)
            AS rt_gap

      FROM tmp_cc_ids JOIN consensus_compound cc USING(consensus_compound_id),
         om 
      )

    SELECT consensus_compound_id
    FROM detailed
    ORDER BY
      has_cpds_in_common DESC,
      (CASE WHEN abs_ppm_error < 4 THEN 0::real ELSE abs_ppm_error END),
      rt_gap
    LIMIT 1;
$$;


-----------------------------------------------------------------------
--  Generate a rectangle (box) to represent the range of mass / retention time
--  for a compound.
--  This rectangle can be simply tested for overlap with other compounds.
--  Include in consensus_compound and obs_mass_no_cons rows.
-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION get_compound_mass_rt_rectangle(
    p_mass                   double precision,
    p_rt_start               real,
    p_rt_end                 real,
    p_rt_adjustment_to_consensus  real,
    p_rect_mass_ppm_width    real,
    p_rect_rt_min_width      real
    )
    RETURNS box
    LANGUAGE SQL AS $$

  WITH pass_1 AS (
    SELECT
      -----------------------------------------------------------------------
      -- How much should the width of the retention time window be expanded?
      -- If (end-start) + adjustment is less than the minimum rt width, then
      -- add enough to expand it to the minimum width.
      -----------------------------------------------------------------------
      GREATEST($6 - ($3 - $2 + ABS(COALESCE($4,0::real))), 0::real)
          AS rt_expansion,

      -----------------------------------------------------------------------
      -- Use the rt adjustment to widen the window rather than shift it,
      -- erring in the direction of being more inclusive.
      -- If the adjustment is negative, shift rt_start to the left,
      -- otherwise shift rt_end to the right.
      -----------------------------------------------------------------------
      $2 + LEAST   ($4, 0::real)  AS rt_start,
      $3 + GREATEST($4, 0::real)  AS rt_end
    )
  SELECT
    -----------------------------------------------------------------------
    --  Add half of the rt expansion to each side of the retention window.
    -----------------------------------------------------------------------
    box(
      point(pass_1.rt_start - pass_1.rt_expansion/2::real,
            $1 * (1 - 1e-6* $5):: double precision),
      point(pass_1.rt_end + pass_1.rt_expansion/2::real,
            $1 * (1 + 1e-6* $5):: double precision)
      )
  FROM pass_1;
$$;


-----------------------------------------------------------------------
--  Map observed masses to :
--    1) consensus compounds.  In most cases, there are zero or one
--       consensus_compound rows for each observed_mass row.
--       However, occasionally there are multiples, and (for simplicity's sake)
--       they must be resolved to one.
--    2) matches (compound identification).  There could be zero, 1, or more.
-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION map_msms_dataset(
    p_consensus_id           varchar,
    p_match_run_id           integer,
    p_min_score              real)
    RETURNS TABLE (
      obs_mass_id            integer,
      mass_rt_rectangle      box,
      match_ids              integer[],
      consensus_compound_id  integer
      )
    LANGUAGE SQL AS $$

  WITH
    tmp_parms AS (
      SELECT rect_mass_ppm_width, rect_rt_min_width
      FROM sample_consensus
      JOIN consensus_parameters USING(consensus_parameters_key)
      WHERE consensus_id = $1
      ),
    tmp_om AS (
      SELECT
        id,
        get_compound_mass_rt_rectangle(mass,
                rt_start, rt_end, rt_adjustment_to_consensus,
                t.rect_mass_ppm_width, t.rect_rt_min_width)
            AS mass_rt_rectangle
      FROM observed_mass o, tmp_parms t
      WHERE o.dataset = (SELECT dataset FROM match_run WHERE match_run_id = $2)
    ),
    tmp_match AS (
      SELECT
        m.obs_mass_id,
        array_unique_vals(array_agg(m.id)) AS match_ids,
        array_unique_vals(array_agg(m.compound_w_origin_id))
            AS compound_w_origin_ids
      FROM match m 
      WHERE match_run_id = $2 AND score >= $3
      GROUP BY m.obs_mass_id
      )
  SELECT
    o.id AS obs_mass_id,
    (array_agg(o.mass_rt_rectangle))[1]
        AS mass_rt_rectangle, -- can't put in GROUP BY since no equal operator
    m.match_ids,
    CASE WHEN count(*) <= 1  THEN (array_agg(c.consensus_compound_id))[1]
         ELSE closest_consensus_cpd(
                  o.id,
                  m.compound_w_origin_ids,
                  array_agg(c.consensus_compound_id)
                  )
    END
        AS consensus_compound_id
  FROM tmp_om o
  LEFT OUTER JOIN tmp_match m ON (o.id = m.obs_mass_id)
  LEFT OUTER JOIN consensus_compound c 
       ON (o.mass_rt_rectangle && c.mass_rt_rectangle)
  GROUP BY o.id, m.match_ids, m.compound_w_origin_ids;
$$;


------------------------------------------------------------------------
--  Handle GPM matches that don't correspond to compound consensus:
--
--   1) Create standard temp input tables for the ids clustering function.
--   2) Cluster (nearly) overlapping peaks (with the 'cluster_ids' function)
--      into the 'tmp_clust' table inside the query.
--   3) Link the ids from the clusters to information about the peaks
--      (contained in 'tmp_map').
--   4) Combine all the observed masses for each cluster into one consensus.
--
--  Input table is tmp_map.
------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION get_clustered_orphan_matches()
    RETURNS TABLE (
      consensus_compound_id  integer,
      obs_mass_list          integer[],
      match_ids              integer[],
      mass                   double precision,
      rt                     real,
      rt_width_at_half_ht    real,
      rt_start               real,
      rt_end                 real,
      quantity               real,
      rel_quantity           real,
      mass_rt_rectangle      box
    )
    LANGUAGE plpgsql AS $$
BEGIN
  CREATE TEMP TABLE tmp_match_no_cons AS
    SELECT m.obs_mass_id, m.mass_rt_rectangle
    FROM tmp_map m 
    WHERE array_length(m.match_ids,1) >= 1 AND m.consensus_compound_id IS NULL;

  CREATE TEMP TABLE tmp_all_ids AS
    SELECT obs_mass_id AS id FROM tmp_match_no_cons;

  CREATE TEMP TABLE tmp_pair AS
    WITH pass_1 AS (
      SELECT
        t1.obs_mass_id AS id_1, t2.obs_mass_id AS id_2
      FROM tmp_match_no_cons t1, tmp_match_no_cons t2
      WHERE t1.obs_mass_id < t2.obs_mass_id AND
            t1.mass_rt_rectangle && t2.mass_rt_rectangle
      )
    SELECT id_1, id_2  FROM pass_1
    UNION ALL
    SELECT id_2, id_1  FROM pass_1;
      
  --------------------------------------------------------------------
  --  Combine all observed masses in each cluster together to create
  --  new consensus compounds.
  --------------------------------------------------------------------
  RETURN QUERY 
    WITH
      tmp_clust AS ( SELECT * FROM cluster_ids() ),
      tmp_grouped AS (
        SELECT
          array_unique_vals(array_agg(obs_mass_id)) AS obs_mass_list,
          combine_arrays(m.match_ids)  AS match_ids,
          avg(o.mass) AS mass,
          avg(o.rt)::real AS rt,
          max(o.rt_width_at_half_ht)  AS  rt_width_at_half_ht,
          encompassing_box(array_agg(m.mass_rt_rectangle))
              AS mass_rt_rectangle,
          sum(o.quantity)::real AS quantity,
          sum(o.rel_quantity)::real AS rel_quantity
        FROM tmp_clust c
        JOIN tmp_map m ON (c.id = m.obs_mass_id)
        JOIN observed_mass o ON (o.id = m.obs_mass_id)
        GROUP BY cluster_id
        )
    SELECT
      nextval('consensus_compound_consensus_compound_id_seq')::integer
          AS consensus_compound_id,
      t.obs_mass_list,
      t.match_ids,
      t.mass,
      t.rt,
      t.rt_width_at_half_ht,
      LEAST( (t.mass_rt_rectangle[0])[0], (t.mass_rt_rectangle[1])[0])::real
          AS rt_start,
      GREATEST((t.mass_rt_rectangle[0])[0], (t.mass_rt_rectangle[1])[0])::real
          AS rt_end,

      t.quantity,
      t.rel_quantity,
      t.mass_rt_rectangle
    FROM tmp_grouped t;

  DROP TABLE tmp_match_no_cons, tmp_all_ids, tmp_pair;
END;
$$;


-----------------------------------------------------------------------
--  For an MS/MS dataset, update the compound consensus rows to represent
--  GPM matches and failed attempts to fragment or identify compounds.
--
--  For matches that don't correspond to existing consensus, create new
--  consensus rows.
--
--  Keep a list of other mz/rt peaks (without GPM identifications) that
--  don't correspond to the consensus; they will be put on exclusion lists.
-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION apply_matches_to_consensus(
    p_consensus_id           varchar,
    p_match_run_id           integer)
    RETURNS void
    LANGUAGE plpgsql AS $$
DECLARE
  v_dataset                  varchar;
  v_min_score                real;
BEGIN
  SELECT dataset INTO v_dataset
  FROM match_run 
  WHERE match_run_id = p_match_run_id;

  SELECT min_score INTO v_min_score
  FROM sample_consensus
  JOIN cons_matching_params USING(cons_matching_params_key)
  WHERE consensus_id = p_consensus_id;

  CREATE TEMP TABLE tmp_map AS 
    SELECT *
    FROM map_msms_dataset(p_consensus_id, p_match_run_id, v_min_score);

  ------------------------------------------------------------------------
  --  Get GPM matches that don't correspond to existing consensus compounds.
  ------------------------------------------------------------------------
  CREATE TEMP TABLE tmp_match_2_cons AS
    SELECT * FROM get_clustered_orphan_matches();

  ------------------------------------------------------------------------
  --  Find old obs_mass_no_cons rows that overlap new consensus compounds
  --  that are to be created from the matches.
  --  Promote those observed masses to be associate with consensus compounds.
  ------------------------------------------------------------------------
  CREATE TEMP TABLE tmp_need_promoting AS
    SELECT
      consensus_compound_id,
      obs_mass_id,
      nc.rt_start,
      nc.rt_end,
      nc.mass_rt_rectangle,
      nc.min_z,
      nc.max_z
    FROM tmp_match_2_cons t
    JOIN obs_mass_no_cons nc ON (t.mass_rt_rectangle && nc.mass_rt_rectangle);

  INSERT INTO obs_mass_promoted(obs_mass_id, match_run_id,
                                consensus_compound_id)
    SELECT obs_mass_id, p_match_run_id, consensus_compound_id
    FROM tmp_need_promoting;

  DELETE FROM obs_mass_no_cons nc
  USING tmp_need_promoting t
  WHERE nc.obs_mass_id = t.obs_mass_id;

  INSERT INTO consensus_cpd_map(consensus_compound_id, consensus_id,
                                obs_mass_id, min_z, max_z, is_fragmented)
    SELECT
      consensus_compound_id, p_consensus_id, obs_mass_id, min_z, max_z, 
      (obs_mass_times_fragmented(obs_mass_id) > 0)  AS is_fragmented
    FROM tmp_need_promoting;

  DROP TABLE tmp_need_promoting;

  ------------------------------------------------------------------------
  --  Observed masses that don't correspond to consensus compound don't
  --  correspond to an existing compound AND there are no GPM matches
  --  (because new consensus compounds will be created for matches).
  ------------------------------------------------------------------------
  INSERT INTO obs_mass_no_cons(obs_mass_id, consensus_id)
    SELECT obs_mass_id, p_consensus_id
    FROM tmp_map
    WHERE consensus_compound_id IS NULL AND
          COALESCE(array_length(match_ids,1),0) = 0;
    

  INSERT INTO consensus_compound(
      consensus_compound_id,
      consensus_id,
      mass,
      rt,
      rt_width_at_half_ht,
      rt_start,
      rt_end,
      quantity,
      rel_quantity,
      mass_rt_rectangle,
      obs_mass_list,
      dataset_list,
      msms_only)
    SELECT
      consensus_compound_id,
      p_consensus_id,
      mass,
      rt,
      rt_width_at_half_ht,
      rt_start,
      rt_end,
      quantity,
      rel_quantity,
      mass_rt_rectangle,
      obs_mass_list,
      (WITH x(id) AS (SELECT unnest(obs_mass_list))
       SELECT array_unique_vals(array_agg(dataset))
       FROM x JOIN observed_mass o USING(id)
      )  AS dataset_list,
      True
    FROM tmp_match_2_cons;

  INSERT INTO consensus_cpd_map(consensus_compound_id, consensus_id,
          obs_mass_id, min_z, max_z, is_fragmented,
          match_id, is_decoy, compound_w_origin_id, score)
    WITH tmp_spectrum AS (
      SELECT
        o.id AS obs_mass_id,
        min(charge) AS min_z,
        max(charge) AS max_z,
        (max(num_peaks) > 0) AS is_fragmented
      FROM component_spectrum cs
      JOIN observed_mass o
          ON (o.dataset = v_dataset AND o.peak_id = cs.parent_peak_id)
      WHERE cs.dataset = v_dataset
      GROUP BY o.id
      ),
    pass_1 AS (
      -------------------------------------------------------------------
      -- GPM matches with newly created consensus compounds
      -------------------------------------------------------------------
      SELECT
        consensus_compound_id,
        unnest(match_ids) AS match_id
      FROM tmp_match_2_cons
      UNION ALL
      -------------------------------------------------------------------
      --  Existing consensus compounds with GPM matches
      -------------------------------------------------------------------
      SELECT
        consensus_compound_id,
        unnest(match_ids) AS match_id
      FROM tmp_map
      WHERE consensus_compound_id IS NOT NULL
      ),
    pass_2 AS (
      SELECT
        consensus_compound_id,
        obs_mass_id,
        match_id,
        is_decoy,
        compound_w_origin_id,
        score
      FROM pass_1 JOIN match m ON (m.id = match_id)
      UNION ALL
      SELECT
        consensus_compound_id,
        obs_mass_id,
        NULL AS match_id,
        NULL AS is_decoy,
        NULL AS compound_w_origin_id,
        NULL AS score
      FROM tmp_map
      WHERE consensus_compound_id IS NOT NULL AND
            COALESCE(array_length(match_ids,1),0) = 0
      )
    SELECT
      consensus_compound_id,
      p_consensus_id,
      obs_mass_id,
      min_z,
      max_z,
      is_fragmented,
      match_id,
      is_decoy,
      compound_w_origin_id,
      score
    FROM pass_2 JOIN tmp_spectrum USING(obs_mass_id);

  DROP TABLE tmp_map, tmp_match_2_cons;

  UPDATE sample_consensus
  SET list_produced = False,
      list_iteration = list_iteration + 1
  WHERE consensus_id = p_consensus_id;

END;
$$;


--------------------------------------------------------------------------
--  Returns the position (one-based index) of the peptide in the protein.
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION get_peptide_position(
    p_cons_matching_params_key varchar,
    p_protein_name           varchar,
    p_decoy                  boolean,
    p_peptide_seq            varchar)
    RETURNS integer
    LANGUAGE SQL AS $$

  SELECT strpos(aa_seq, $4)
  FROM cons_matching_protein
  WHERE cons_matching_params_key = $1 AND protein_name = $2 AND decoy = $3;

$$;


--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION consensus_remove_match(
    p_consensus_compound_id  integer,
    p_consensus_id           varchar,
    p_match_id               integer)
    RETURNS void
    LANGUAGE SQL AS $$

  DELETE FROM consensus_cpd_map
  WHERE consensus_compound_id = $1 AND consensus_id = $2 AND match_id = $3;

  INSERT INTO consensus_removed_match(consensus_compound_id, consensus_id,
                                      match_id)
    VALUES($1,$2,$3);
$$;


--------------------------------------------------------------------------
--
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION msms_guidance_report(
    p_consensus_id           varchar)
    RETURNS TABLE(
      mz                     double precision,
      delta_mz_ppm           real,
      charge                 integer,
      rt                     real,
      delta_rt               real,
      exclude                boolean
      )
    LANGUAGE plpgsql AS $$
DECLARE
  v_cons_parms               consensus_parameters;
  v_lc_conf                  lc_configuration;
  v_max_excl_pref_list_length integer;
  v_num_calibrants           integer;
  v_num_exclude_peaks        integer;
  v_list_produced            boolean;
  v_list_iteration           integer;
  v_excl_list_len            integer;
  v_pref_list_len            integer;
BEGIN
  SELECT list_iteration INTO v_list_iteration
  FROM sample_consensus
  WHERE consensus_id = p_consensus_id;

  SELECT * INTO v_cons_parms
  FROM consensus_parameters
  WHERE consensus_parameters_key =
        (SELECT consensus_parameters_key FROM sample_consensus
         WHERE consensus_id = p_consensus_id);

  SELECT * INTO v_lc_conf
  FROM lc_configuration
  WHERE lc_config_id = v_cons_parms.lc_config_id;

  -----------------------------------------------------------------------
  --  Always exclude the calibrant ions.
  -----------------------------------------------------------------------
  RETURN QUERY 
    SELECT
      ci.mz,
      2::real*v_lc_conf.calibrant_error_ppm_limit  AS delta_mz_ppm,
      ci.charge,
      (v_lc_conf.theoretical_max_rt/2)::real AS rt,
      v_lc_conf.theoretical_max_rt::real AS delta_rt,
      True AS exclude
    FROM calibrant_ion ci
    WHERE ci.calibrant_ion_list_id = v_lc_conf.calibrant_ion_list_id;

  SELECT count(*) INTO v_num_calibrants
  FROM calibrant_ion
  WHERE calibrant_ion_list_id = v_lc_conf.calibrant_ion_list_id;

  -----------------------------------------------------------------------
  --  Determine the number of items that would be in the exclusion or
  --  preferred list.
  -----------------------------------------------------------------------
  SELECT
    SUM(CASE WHEN exclude_compound THEN max_z-min_z+1 ELSE 0 END)::integer,
    SUM(CASE WHEN exclude_compound THEN 0             ELSE 1 END)::integer
  INTO v_excl_list_len, v_pref_list_len
  FROM consensus_compound
  WHERE consensus_id = p_consensus_id;

  -----------------------------------------------------------------------
  --  If the exclusion list is not too long, use it.
  -----------------------------------------------------------------------
  IF v_list_iteration < v_cons_parms.max_exclusion_iterations  AND
     v_excl_list_len < v_pref_list_len   THEN

    RETURN QUERY
      SELECT
        ((cc.mass/z) + v_lc_conf.charge_carrier_mass)::double precision AS mz,
        (1e6 * height(mass_rt_rectangle)/mass)::real AS delta_mz_ppm,
        z AS charge,
        cc.rt,
        GREATEST(v_cons_parms.min_excl_pref_rt_width, 
                 (width(mass_rt_rectangle))::real)
            AS delta_rt,
        True AS exclude
      FROM consensus_compound cc, generate_series(1, 9) AS z
      WHERE consensus_id = p_consensus_id  AND exclude_compound  AND
            z BETWEEN min_z and max_z  AND
            (cc.mass/z) + v_lc_conf.charge_carrier_mass
              < v_lc_conf.theoretical_max_mz
      ORDER BY quantity DESC
      LIMIT (v_max_excl_pref_list_length - v_num_calibrants);

  -----------------------------------------------------------------------
  --  Otherwise use a preferred list.
  -----------------------------------------------------------------------
  ELSE
    CREATE TEMP TABLE tmp_pref_list AS 
      SELECT
        consensus_compound_id,
        ((cc.mass/cc.dominant_z) + v_lc_conf.charge_carrier_mass)
            ::double precision AS mz,
        (1e6 * height(mass_rt_rectangle)/mass)::real AS delta_mz_ppm,
        cc.dominant_z  AS charge,
        cc.rt,
        GREATEST(v_cons_parms.min_excl_pref_rt_width, 
                 (width(mass_rt_rectangle))::real)
            AS delta_rt
      FROM consensus_compound cc
      WHERE consensus_id = p_consensus_id  AND NOT exclude_compound
      ORDER BY quantity DESC
      LIMIT (v_max_excl_pref_list_length - v_num_calibrants);

    -----------------------------------------------------------------------
    --  Increment the num_times_preferred counter for each compound put in
    --  the preferred list.  However, make sure not to repeat this if the
    --  same report is generated multiple times.
    --
    --  For each compound, mark it as excluded if the limit of the number
    --  of times preferred is reached.
    -----------------------------------------------------------------------
    SELECT list_produced INTO v_list_produced
    FROM sample_consensus
    WHERE consensus_id = p_consensus_id;

    IF NOT v_list_produced THEN
      UPDATE consensus_compound cc
      SET num_times_preferred = num_times_preferred + 1
      FROM tmp_pref_list pf
      WHERE cc.consensus_compound_id = pf.consensus_compound_id;

      UPDATE consensus_compound
      SET exclude_compound = True
      WHERE consensus_id = p_consensus_id  AND
            num_times_preferred + 1 > v_cons_parms.max_times_preferred;
    END IF;

    RETURN QUERY
      SELECT t.mz, t.delta_mz_ppm, t.charge, t.rt, t.delta_rt, False AS exclude
      FROM tmp_pref_list t;

    DROP TABLE tmp_pref_list;

  END IF;

  UPDATE sample_consensus
  SET list_produced = True
  WHERE consensus_id = p_consensus_id;

END;
$$;
