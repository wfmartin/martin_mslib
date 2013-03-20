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
    p_rect_rt_min_height     real
    )
    RETURNS box
    LANGUAGE SQL AS $$

  WITH pass_1 AS (
    SELECT
      CASE
        WHEN (ABS($4) + $3 - $2) < $6  THEN ($6 - ABS($4) + $3 - $2)
        ELSE 0::real
      END  AS rt_expansion,

      CASE WHEN $4 < 0 THEN $2 + $4 ELSE $2 END  AS rt_start,

      CASE WHEN $4 > 0 THEN $3 + $4 ELSE $3 END  AS rt_end
    )
  SELECT
    box(
      point($1 * (1 - 1e-6* $5):: double precision,
            pass_1.rt_start - pass_1.rt_expansion),
      point($1 * (1 + 1e-6* $5):: double precision,
            pass_1.rt_end - pass_1.rt_expansion)
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
      SELECT rect_mass_ppm_width, rect_rt_min_height
      FROM sample_consensus
      JOIN consensus_parameters USING(consensus_parameters_key)
      WHERE consensus_id = $1
      ),
    tmp_om AS (
      SELECT
        id,
        get_compound_mass_rt_rectangle(mass,
                rt_start, rt_end, rt_adjustment_to_consensus,
                t.rect_mass_ppm_width, t.rect_rt_min_height)
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
      LEAST( (t.mass_rt_rectangle[0])[1], (t.mass_rt_rectangle[1])[1])::real
          AS rt_start,
      GREATEST((t.mass_rt_rectangle[0])[1], (t.mass_rt_rectangle[1])[1])::real
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

END;
$$;


--------------------------------------------------------------------------
--  Return large rectangles in the mass x rt plane (constrained to the ranges
--  of interesting masses and retention times) that should be excluded
--  from tandem fragmentation.
--
--  Look for relatively large mass ranges (relative to measurement precision)
--  that are mostly to be excluded, i.e., they contain at most 2 compounds
--  which would break up the rectangle.
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION get_wide_mass_bands(
    p_consensus_id           varchar)
    RETURNS SETOF box
    LANGUAGE plpgsql AS $$
DECLARE
  v_cons_parms               record;
  v_cpds                     cons_cpd_info[];
  v_total_num_cpds           integer;
  v_num_cpds                 integer;
  v_inds_used                integer[];
  v_inds                     integer[];
  v_i                        integer;
  v_lo_mass                  double precision;
  v_hi_mass                  double precision;
  v_lo_rt                    double precision;
  v_hi_rt                    double precision;
BEGIN
  SELECT * INTO v_cons_parms
  FROM sample_consensus
  JOIN consensus_parameters USING(consensus_parameters_key)
  WHERE consensus_id = p_consensus_id;

  SELECT
    array_agg((consensus_compound_id,mass,rt_start,rt_end)::cons_cpd_info
        ORDER BY mass),
    count(*)
  INTO v_cpds, v_total_num_cpds
  FROM consensus_compound
  WHERE consensus_id = p_consensus_id AND
        mass BETWEEN v_cons_parms.min_mass AND v_cons_parms.max_mass;

  FOR v_i, v_num_cpds IN
    WITH 
      -------------------------------------------------------------------
      --  num_cpds         : Number of compounds disrupting band
      --  min_mass_ppm_gap : Minimum mass width of band, depends upon num_cpds
      -------------------------------------------------------------------
      tmp_lim(num_cpds, min_mass_ppm_gap) AS (VALUES
        (0,             v_cons_parms.band_mass_width_ppm),
        (1, 1.3::real * v_cons_parms.band_mass_width_ppm),
        (2, 1.5::real * v_cons_parms.band_mass_width_ppm)
        ),
      -------------------------------------------------------------------
      --  Generate all combinations of starting index (i) and number of 
      --  compound masses in the band, plus for each of those:
      --   1) mass ppm difference between lowest and highest masses
      --   2) max ppm mass difference for number of compounds
      -------------------------------------------------------------------
      tmp_pass_1 AS (
        SELECT
          i,
          num_cpds,
          min_mass_ppm_gap,
          (1e6*(v_cpds[i+num_cpds+1].mass - v_cpds[i].mass)/v_cpds[i].mass)
              AS mass_ppm_diff
        FROM tmp_lim,
             generate_series(0, v_total_num_cpds-2) AS i
        )
      SELECT i, num_cpds
      FROM tmp_pass_1
      WHERE mass_ppm_diff >= min_mass_ppm_gap
      ORDER BY (mass_ppm_diff/min_mass_ppm_gap) DESC
  LOOP
    --------------------------------------------------------------------
    --  To efficiently check for band overlaps, the range of indices
    --  (inclusive of the starting index, exclusive of the end index) is
    --  compared to the accumulation of such indices from previous iterations.
    --------------------------------------------------------------------
    v_inds := array[]::integer[];
    FOR i IN 0 .. v_num_cpds  LOOP
      v_inds[i] := v_i+i;
    END LOOP;

    --  Skip if the some of the given compounds have already been handled.
    CONTINUE WHEN COALESCE(v_inds && v_inds_used, False);

    --------------------------------------------------------------------------
    --  The bounds of the mass bands.
    --------------------------------------------------------------------------
    v_lo_mass := v_cpds[v_i].mass
                 * (1+1e-6*v_cons_parms.band_mass_margin_ppm);
    v_hi_mass := v_cpds[v_i+v_num_cpds+1].mass
                 * (1-1e-6*v_cons_parms.band_mass_margin_ppm);

    IF v_num_cpds = 2 THEN
      ------------------------------------------------------------------------
      -- Combine retention time ranges if they overlap; treat as num_cpds = 1
      ------------------------------------------------------------------------
      IF ranges_overlap(v_cpds[v_i+1].rt_start, v_cpds[v_i+1].rt_end,
          v_cpds[v_i+2].rt_start - v_cons_parms.band_rt_margin,
          v_cpds[v_i+2].rt_end   + v_cons_parms.band_rt_margin)   THEN

        -- upper bound for rt starting at min_rt
        v_lo_rt := LEAST(v_cpds[v_i+1].rt_start, v_cpds[v_i+2].rt_start)
            - v_cons_parms.band_rt_margin;
        -- lower bound for rt ending at max_rt
        v_hi_rt := GREATEST(v_cpds[v_i+1].rt_end, v_cpds[v_i+2].rt_end)
            + v_cons_parms.band_rt_margin;
        v_num_cpds := 1;

      ELSE  -- trifurcate region
        v_lo_rt := v_cons_parms.min_rt;
        v_hi_rt := LEAST(v_cpds[v_i+1].rt_start, v_cpds[v_i+2].rt_start)
            - v_cons_parms.band_rt_margin;
        IF (v_hi_rt - v_lo_rt >= v_cons_parms.band_min_rt_height) THEN
          RETURN NEXT box(
            point(v_lo_mass, v_lo_rt),
            point(v_hi_mass, v_hi_rt)
            );
        END IF;

        v_lo_rt := LEAST(v_cpds[v_i+1].rt_end, v_cpds[v_i+2].rt_end)
                + v_cons_parms.band_rt_margin;
        v_hi_rt := GREATEST(v_cpds[v_i+1].rt_start, v_cpds[v_i+2].rt_start)
                - v_cons_parms.band_rt_margin;

        IF (v_hi_rt - v_lo_rt >= v_cons_parms.band_min_rt_height) THEN
          RETURN NEXT box(
            point(v_lo_mass, v_lo_rt),
            point(v_hi_mass, v_hi_rt)
            );
        END IF;

        v_lo_rt := GREATEST(v_cpds[v_i+1].rt_end, v_cpds[v_i+2].rt_end)
            + v_cons_parms.band_rt_margin;
        v_hi_rt := v_cons_parms.max_rt;

        IF (v_hi_rt - v_lo_rt >= v_cons_parms.band_min_rt_height) THEN
          RETURN NEXT box(
            point(v_lo_mass, v_lo_rt),
            point(v_hi_mass, v_hi_rt)
            );
        END IF;
      END IF;
    END IF;

    IF v_num_cpds = 1 THEN  -- bifurcate region
      -- Low and high rt for the compound that breaks up the band:
      v_lo_rt := COALESCE(v_lo_rt,
              v_cpds[v_i+1].rt_start + v_cons_parms.band_rt_margin);
      v_hi_rt := COALESCE(v_hi_rt,
              v_cpds[v_i+1].rt_end   - v_cons_parms.band_rt_margin);

      IF (v_lo_rt - v_cons_parms.min_rt >= v_cons_parms.band_min_rt_height)
      THEN
        RETURN NEXT box(
          point(v_lo_mass, v_cons_parms.min_rt::double precision),
          point(v_hi_mass, v_lo_rt)
          );
      END IF;

      IF (v_cons_parms.max_rt - v_hi_rt >= v_cons_parms.band_min_rt_height)
      THEN
        RETURN NEXT box(
          point(v_lo_mass, v_hi_rt),
          point(v_hi_mass, v_cons_parms.max_rt::double precision)
          );
      END IF;

    ELSIF v_num_cpds = 0 THEN
      RETURN NEXT box(
        point(v_lo_mass, v_cons_parms.min_rt::double precision),
        point(v_hi_mass, v_cons_parms.max_rt::double precision)
        );
    END IF;
   
    v_inds_used := v_inds_used || v_inds;

  END LOOP;

END;
$$;


--------------------------------------------------------------------------
--  Determine which mass/charge combinations are to be excluded in the
--  next mass run.
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION determine_ms_exclusions(
    p_consensus_id           varchar)
    RETURNS TABLE(
      excluded_area          box,
      min_z                  integer,
      max_z                  integer,
      excl_type              varchar
    )
    LANGUAGE plpgsql AS $$
DECLARE
  v_default_min_z            integer;
  v_default_max_z            integer;
  v_redundant_overlap        real;
  v_origin                   point;
BEGIN
  SELECT redundant_overlap, default_min_z, default_max_z
  INTO v_redundant_overlap, v_default_min_z, v_default_max_z
  FROM sample_consensus
  JOIN consensus_parameters cp USING(consensus_parameters_key)
  JOIN lc_configuration USING(lc_config_id)
  WHERE consensus_id = p_consensus_id;

  ---------------------------------------------------------------------------
  --  Put regions for excluded compounds in the tmp_comp table.
  --  Make wide mass bands of excluded rectangles in mass X rt plane.
  ---------------------------------------------------------------------------
  CREATE TEMP SEQUENCE tmp_seq;  -- temp ids involved in clustering

  CREATE TEMP TABLE tmp_comp AS
    SELECT
      nextval('tmp_seq')::integer AS comp_id,
      mass_rt_rectangle,
      o.min_z,
      o.max_z
    FROM obs_mass_no_cons o
    WHERE consensus_id = p_consensus_id
    UNION ALL
    SELECT
      nextval('tmp_seq')::integer AS comp_id,
      mass_rt_rectangle,
      cc.min_z,
      cc.max_z
    FROM consensus_compound cc
    WHERE consensus_id = p_consensus_id AND exclude_compound;

  --  Make wide mass bands of excluded rectangles in mass X rt plane.
  CREATE TEMP TABLE tmp_bands AS 
    SELECT
      nextval('tmp_seq')::integer AS band_id,
      get_wide_mass_bands(p_consensus_id) AS excluded_area;

  DROP SEQUENCE tmp_seq;

  -----------------------------------------------------------------------
  --  Join mass bands to individual excluded compounds to get the charges,
  --  and also to determine which compounds are separate from the mass bands.
  -----------------------------------------------------------------------
  RETURN QUERY 
    SELECT
      (array_agg(b.excluded_area))[1] AS excluded_area,
      COALESCE(min(c.min_z), v_default_min_z)  AS min_z,
      COALESCE(max(c.max_z), v_default_max_z)  AS max_z,
      'mass_band'::varchar  AS excl_type
    FROM tmp_bands b LEFT OUTER JOIN tmp_comp c
         ON( b.excluded_area && mass_rt_rectangle)
    GROUP BY band_id;

  -----------------------------------------------------------------------
  --  Join to the mass bands to figure out which compounds are redundant;
  --  delete the redundant compounds.
  -----------------------------------------------------------------------
  DELETE FROM tmp_comp c
  USING tmp_bands b
  WHERE box_overlap_proportion(mass_rt_rectangle, b.excluded_area)
        >= v_redundant_overlap;

  ---------------------------------------------------------------------------
  --  Cluster compounds by overlaps and combine them.
  ---------------------------------------------------------------------------
  CREATE TEMP TABLE tmp_pair AS
    WITH pass_1 AS (
      SELECT t1.comp_id AS id_1, t2.comp_id AS id_2
      FROM tmp_comp t1, tmp_comp t2
      WHERE t1.comp_id < t2.comp_id  AND
            t1.mass_rt_rectangle && t2.mass_rt_rectangle
      )
    SELECT id_1, id_2
    FROM pass_1
    UNION
    SELECT id_2, id_1
    FROM pass_1;

  CREATE TEMP TABLE tmp_all_ids AS SELECT comp_id AS id FROM tmp_comp;
  CREATE TEMP TABLE tmp_cluster AS SELECT * FROM cluster_ids();

  RETURN QUERY
    SELECT
      encompassing_box(array_agg(mass_rt_rectangle))  AS excluded_area,
      min(c.min_z) AS min_z,
      max(c.max_z) AS max_z,
      'obs_mass_no_cons'::varchar  AS excl_type
    FROM tmp_cluster JOIN tmp_comp c ON(id = comp_id)
    GROUP BY cluster_id;

  DROP TABLE tmp_comp, tmp_bands, tmp_pair, tmp_all_ids, tmp_cluster;
END;
$$;


--------------------------------------------------------------------------
--
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION generate_exclusion_list(
    p_consensus_id           varchar)
    RETURNS integer
    LANGUAGE plpgsql AS $$
DECLARE
  v_iteration                integer;
  v_bounding_area            real;
BEGIN
  CREATE TEMP TABLE tmp_excl AS
    SELECT * FROM determine_ms_exclusions(p_consensus_id);

  SELECT COALESCE(max(iteration), 0) + 1, (array_agg(bounding_area))[1]
  INTO v_iteration, v_bounding_area
  FROM ms_exclusion_iteration
  WHERE consensus_id = p_consensus_id;

  IF v_bounding_area IS NULL THEN
    SELECT (max_mass-min_mass)*(max_rt-min_rt)
    INTO v_bounding_area
    FROM sample_consensus 
    JOIN consensus_parameters USING(consensus_parameters_key)
    WHERE consensus_id = p_consensus_id;
  END IF;

  INSERT INTO ms_exclusion(consensus_id, iteration, excluded_area,
                           min_z, max_z, excl_type)
    SELECT
      p_consensus_id,
      v_iteration,
      excluded_area,
      min_z,
      max_z,
      excl_type
    FROM tmp_excl;

  INSERT INTO ms_exclusion_iteration(consensus_id, iteration, area_count,
                                     total_area, bounding_area)
    SELECT
      p_consensus_id,
      v_iteration,
      count(*) AS area_count,
      sum(area(excluded_area))::real AS total_area,
      v_bounding_area
    FROM tmp_excl;

  RETURN v_iteration;
END;
$$;


--------------------------------------------------------------------------
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION get_perimeter_exclusions(
    p_lc_config_id           varchar,
    p_min_mass               real,   -- range of interesting masses
    p_max_mass               real,
    p_min_rt                 real,   -- range of interesting rt
    p_max_rt                 real,
    p_min_z                  integer,
    p_max_z                  integer)
    RETURNS TABLE(
      mz                     real,
      delta_mz_ppm           real,
      charge                 integer,
      rt                     real,
      delta_rt               real
      )
    LANGUAGE plpgsql AS $$
DECLARE
  v_lc_conf                  lc_configuration;
  v_mz_lower_bound           real;
BEGIN
  SELECT * INTO v_lc_conf
  FROM lc_configuration
  WHERE lc_config_id = p_lc_config_id;

  ------------------------------------------------------------------------
  -- Scenarios of full mz range:
  --  1) rt < interesting rt range
  --  2) rt > interesting rt range
  ------------------------------------------------------------------------
  RETURN QUERY
  WITH tmp_rt(rt, delta_rt) AS ( VALUES
     ( p_min_rt/2::real, p_min_rt/2::real),
     ( (v_lc_conf.theoretical_max_rt + p_max_rt)/2::real,
       (v_lc_conf.theoretical_max_rt - p_max_rt)/2::real )
     )
  SELECT 
    v_lc_conf.theoretical_max_mz::real/2::real  AS mz,
    5e5::real  AS delta_mz_ppm,  -- (50%)
    z AS charge,
    tmp_rt.rt  AS rt,
    tmp_rt.delta_rt  AS delta_rt
  FROM tmp_rt
  CROSS JOIN generate_series(p_min_z, p_max_z) AS z;

  ------------------------------------------------------------------------
  --  Scenarios of full rt range:
  --  mz values must be computed for all charges
  ------------------------------------------------------------------------
  rt := v_lc_conf.theoretical_max_rt/2::real;
  delta_rt := rt;

  FOR z IN p_min_z .. p_max_z  LOOP
    charge := z;
    ------------------------------------------------------------------------
    --  mass < interesting masses
    --  (mz is based upon midpoint between zero and minimum interesting mass)
    ------------------------------------------------------------------------
    mz := 0.5 * p_min_mass/charge + v_lc_conf.charge_carrier_mass;
    delta_mz_ppm := 5e5;  -- ( 50% )
    RETURN NEXT;

    ------------------------------------------------------------------------
    --  mass > interesting masses
    --  (mz is based upon midpoint maximum interesting mass and 
    --   instrument's max mz sensitivity)
    --  rt is full range
    ------------------------------------------------------------------------
    v_mz_lower_bound := p_max_mass/charge + v_lc_conf.charge_carrier_mass;
    IF v_mz_lower_bound < v_lc_conf.theoretical_max_mz  THEN
      mz := (v_mz_lower_bound + v_lc_conf.theoretical_max_mz)/2;
      delta_mz_ppm :=
          1e6 * (v_lc_conf.theoretical_max_mz - v_mz_lower_bound )/(2*mz);
      RETURN NEXT;
    END IF;
  END LOOP;

END;
$$;


--------------------------------------------------------------------------
--  Get the  mz / z / rt  values to exclude from the next ms/ms run
--  (not including the calibration ions).
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION ms_exclusions_report(
    p_consensus_id           varchar,
    p_iteration              integer)
    RETURNS TABLE(
      mz                     real,
      delta_mz_ppm           real,
      charge                 integer,
      rt                     real,
      delta_rt               real
      )
    LANGUAGE plpgsql AS $$
DECLARE
  v_iteration                integer;
  v_cons_parms               consensus_parameters;
  v_lc_config                lc_configuration;
BEGIN
  v_iteration := p_iteration;

  SELECT * INTO v_cons_parms
  FROM consensus_parameters
  WHERE consensus_parameters_key =
        (SELECT consensus_parameters_key FROM sample_consensus
         WHERE consensus_id = p_consensus_id);

  SELECT * INTO v_lc_config
  FROM lc_configuration 
  WHERE lc_config_id = v_cons_parms.lc_config_id;

  IF v_iteration IS NULL THEN
    SELECT max(iteration) INTO v_iteration
    FROM ms_exclusion_iteration
    WHERE consensus_id = p_consensus_id;
  END IF;

  -----------------------------------------------------------------------
  --  Calibrant ions
  -----------------------------------------------------------------------
  RETURN QUERY
    SELECT
      ci.mz::real,
      v_lc_config.calibrant_error_ppm_limit,
      ci.charge,
      (v_lc_config.theoretical_max_rt/2)::real AS rt,
      (v_lc_config.theoretical_max_rt/2)::real AS delta_rt
    FROM calibrant_ion ci
    WHERE calibrant_ion_list_id = v_lc_config.calibrant_ion_list_id;

  -----------------------------------------------------------------------
  --  Perimeter (outside of interesting masses or retention times).
  -----------------------------------------------------------------------
  RETURN QUERY
    SELECT * FROM get_perimeter_exclusions(
        v_cons_parms.lc_config_id,
        v_cons_parms.min_mass,
        v_cons_parms.max_mass,
        v_cons_parms.min_rt,
        v_cons_parms.max_rt,
        v_cons_parms.default_min_z,
        v_cons_parms.default_max_z);
   
  -----------------------------------------------------------------------
  --  Regions informed by MS/MS runs and GPM identifications, as stored in
  --  the ms_exclusion table.
  -----------------------------------------------------------------------
  IF v_iteration IS NOT NULL THEN
    RETURN QUERY 
      WITH pass_1 AS (
        SELECT
          ( (center(excluded_area))[0]/z
             + v_lc_config.charge_carrier_mass)::real  AS mz,
          (5e5 * width(excluded_area)/(center(excluded_area))[0])::real
              AS delta_mz_ppm,
          z AS charge,
          (center(excluded_area))[1]::real AS rt,
          (0.5 * height(excluded_area))::real  AS delta_rt
        FROM ms_exclusion CROSS JOIN generate_series(1, 9) AS z
        WHERE consensus_id = p_consensus_id  AND iteration = v_iteration  AND
              (center(excluded_area))[0]::real < v_cons_parms.max_mass  AND
              z BETWEEN min_z and max_z
        )
      SELECT p.mz, p.delta_mz_ppm, p.charge, p.rt, p.delta_rt
      FROM pass_1 p
      WHERE p.mz < v_lc_config.theoretical_max_mz;
  END IF;
END;
$$;


--------------------------------------------------------------------------
--  Get the MS exclusions for the latest run.
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION ms_exclusions_report(
    p_consensus_id           varchar)
    RETURNS TABLE(
      mz                     real,
      delta_mz_ppm           real,
      charge                 integer,
      rt                     real,
      delta_rt               real
      )
    LANGUAGE SQL AS $$

  SELECT * FROM ms_exclusions_report($1, NULL);

$$;


--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION ms_inclusions_report(
    p_consensus_id           varchar)
    RETURNS TABLE(
      mz                     real,
      charge                 integer,
      rt                     real,
      delta_rt               real
      )
    LANGUAGE SQL AS $$

  WITH t AS (
    SELECT charge_carrier_mass
    FROM sample_consensus
    JOIN consensus_parameters USING(consensus_parameters_key)
    JOIN lc_configuration USING(lc_config_id)
    WHERE consensus_id = $1
    )
  SELECT 
      (mass/z + t.charge_carrier_mass)::real AS mz,
      z AS charge,
      rt,
      (rt_end - rt_start)/2::real  AS delta_rt
  FROM consensus_compound CROSS JOIN generate_series(1, 9) AS z,  t
  WHERE consensus_id = $1  AND (NOT exclude_compound) AND
        z BETWEEN min_z and max_z;
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

