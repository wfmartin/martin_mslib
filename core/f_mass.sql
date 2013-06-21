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
------------------------------------------------------------------------
--  Given the table tmp_cluster, containing
--    cluster_id  integer,
--    id          integer,
--
--  If there's one input dataset, then:
--    1) Singleton peaks retain their peak_id.
--    2) Quantities are summed.
--  Otherwise (multiple datasets combined):
--    1) New peak_ids are generated for all rows.
--    2) Quantities are averaged.
------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION generate_merged_peaks(
    p_old_dataset            varchar)
    RETURNS TABLE (
      cluster_id             integer,
      peak_id                integer,
      mass                   double precision,
      rt                     real,
      rt_width_at_half_ht    real,
      rt_start               real,
      rt_end                 real,
      min_z                  integer,
      max_z                  integer,
      dominant_z             integer,
      quantity               real,
      rel_quantity           real,
      merged_ids             integer[]
      )
    LANGUAGE plpgsql AS $$
DECLARE
  v_max_peak_id              integer;
BEGIN
  IF p_old_dataset IS NOT NULL THEN
    SELECT max(o.peak_id) INTO v_max_peak_id
    FROM observed_mass o
    WHERE dataset = p_old_dataset;
  ELSE
    v_max_peak_id = 0;
  END IF;

  EXECUTE 'CREATE TEMP SEQUENCE peak_id_seq START WITH ' || (v_max_peak_id+1);

  ------------------------------------------------------------------------
  --  New observed_mass rows:
  --   1) For peaks that don't cluster, preserve the peak_id,
  --      otherwise generate a new one.
  --   2) Use the quantity-weighted averages for mass and retention time.
  --   3) Take the start/end range that encompasses all start/end ranges.
  --   4) Use the maximum width of all peaks.
  ------------------------------------------------------------------------
  RETURN QUERY 
  SELECT
    c.cluster_id,
    CASE WHEN count(*) = 1 THEN min(o.peak_id)  -- min() does nothing
         ELSE nextval('peak_id_seq')::integer
    END  AS peak_id,
    avg(o.mass) AS mass,
    avg(o.rt + COALESCE(rt_adjustment_to_consensus,0::real))::real   AS rt,
    max(o.rt_width_at_half_ht) AS rt_width_at_half_ht,
    min(o.rt_start + COALESCE(rt_adjustment_to_consensus,0::real))
        AS rt_start,
    max(o.rt_end + COALESCE(rt_adjustment_to_consensus,0::real))
        AS rt_end,
    min(o.min_z) AS min_z,
    max(o.max_z) AS max_z,

    (array_agg(o.dominant_z ORDER BY o.quantity DESC))[1] AS dominant_z,

    CASE WHEN p_old_dataset IS NULL THEN avg(o.quantity)::real
                                    ELSE sum(o.quantity)::real
    END   AS quantity,

    CASE WHEN p_old_dataset IS NULL THEN avg(o.rel_quantity)::real
                                    ELSE sum(o.rel_quantity)::real
    END   AS rel_quantity,

    array_agg(o.id) AS merged_ids
  FROM tmp_cluster c
  JOIN observed_mass o ON (o.id = c.id)
  GROUP BY c.cluster_id;
  
  DROP SEQUENCE peak_id_seq;

END;
$$;


------------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION obs_mass_times_fragmented(
    p_dataset                varchar,
    p_peak_id                integer)
    RETURNS integer
    LANGUAGE SQL AS $$

  SELECT count(*)::integer
  FROM component_spectrum
  WHERE dataset = $1 AND parent_peak_id = $2

$$;


------------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION obs_mass_times_fragmented(
    p_obs_mass_id            integer)
    RETURNS integer
    LANGUAGE SQL AS $$

  WITH pass_1 AS (SELECT dataset, peak_id FROM observed_mass WHERE id = $1)
  SELECT obs_mass_times_fragmented(dataset, peak_id)
  FROM pass_1;

$$;


------------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION get_dataset_list(
    p_obs_mass_list          integer[])
    RETURNS varchar[]
    LANGUAGE SQL AS $$

  WITH ids(id) AS (SELECT unnest($1))
  SELECT array_unique_vals (array_agg(dataset))
  FROM ids JOIN observed_mass USING(id);
$$;


------------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION determine_merged_clusters(
    p_old_dataset            varchar,
    p_mass_ppm_error_limit   real,
    p_max_rt_gap             real)
    RETURNS TABLE (
      cluster_id             integer,
      id                     integer
    )
    LANGUAGE plpgsql AS $$
BEGIN
  CREATE TEMP TABLE tmp_pair AS
    SELECT o1.id AS id_1, o2.id AS id_2
    FROM observed_mass o1, observed_mass o2
    WHERE o1.dataset = p_old_dataset AND 
          o2.dataset = p_old_dataset AND
          o1.id <> o2.id   AND
          masses_within_error_limit(o1.mass, o2.mass, p_mass_ppm_error_limit)
                AND
          ranges_overlap(o1.rt_start-p_max_rt_gap, o1.rt_end+p_max_rt_gap,
                         o2.rt_start, o2.rt_end)
  ;
  
  CREATE INDEX tmp_pair_index_1 ON tmp_pair(id_1);
  CREATE INDEX tmp_pair_index_2 ON tmp_pair(id_2);

  CREATE TEMP TABLE tmp_all_ids AS 
    SELECT o.id FROM observed_mass o  WHERE dataset = p_old_dataset;

  RETURN QUERY  SELECT * FROM cluster_ids();

  DROP TABLE tmp_pair, tmp_all_ids;
END;
$$;


------------------------------------------------------------------------------
--  Create a new dataset from an existing one, merging together similar peaks
--  (kind of overlapping).
--  Note this doesn't handle component spectra (for Tandem MS).  It's probably
--  only useful for MS-only datasets.
------------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION create_dataset_from_merged_peaks(
    p_old_dataset            varchar,
    p_new_dataset            varchar,
    p_mass_ppm_error_limit   real,
    p_max_rt_gap             real)
    RETURNS void
    LANGUAGE plpgsql AS $$
DECLARE
  v_num_rows                 integer;
BEGIN
  CREATE TEMP TABLE tmp_pair AS
    WITH pass_1 AS (
      SELECT o1.id AS id_1, o2.id AS id_2
      FROM observed_mass o1, observed_mass o2
      WHERE o1.dataset = p_old_dataset AND 
            o2.dataset = p_old_dataset AND
            o1.id < o2.id   AND
            masses_within_error_limit(o1.mass, o2.mass, p_mass_ppm_error_limit)
                  AND
            ranges_overlap(o1.rt_start-p_max_rt_gap, o1.rt_end+p_max_rt_gap,
                           o2.rt_start, o2.rt_end)
      )
  SELECT id_1, id_2 FROM pass_1
  UNION ALL
  SELECT id_2, id_1 FROM pass_1;
  
  CREATE INDEX tmp_pair_index_1 ON tmp_pair(id_1);
  CREATE INDEX tmp_pair_index_2 ON tmp_pair(id_2);

  CREATE TEMP TABLE tmp_all_ids AS 
    SELECT o.id FROM observed_mass o  WHERE dataset = p_old_dataset;

  CREATE TEMP TABLE tmp_cluster AS SELECT * FROM cluster_ids();

  DROP TABLE tmp_pair, tmp_all_ids;

  ------------------------------------------------------------------------
  --  New observed_mass rows:
  --   1) For peaks that don't cluster, preserve the peak_id,
  --      otherwise generate a new one.
  --   2) Use the quantity-weighted averages for mass and retention time.
  --   3) Take the start/end range that encompasses all start/end ranges.
  --   4) Use the maximum width of all peaks.
  ------------------------------------------------------------------------
  INSERT INTO observed_mass(dataset,
      peak_id, mass, rt, rt_width_at_half_ht, rt_start, rt_end,
      quantity, rel_quantity, min_z, max_z, dominant_z, merged_ids)
    SELECT
      p_new_dataset  AS dataset,
      peak_id, mass, rt, rt_width_at_half_ht, rt_start, rt_end,
      quantity, rel_quantity, min_z, max_z, dominant_z, merged_ids
    FROM generate_merged_peaks(p_old_dataset);

  GET DIAGNOSTICS v_num_rows := ROW_COUNT;

  INSERT INTO mass_dataset(
      dataset, 
      merge_source_dataset,
      num_peaks,
      bio_context_id,
      lc_config_id,
      data_source_type,
      frag_data_is_raw,
      rt_adjustment_ts,
      consensus_id)
    SELECT
      p_new_dataset,
      p_old_dataset,
      v_num_rows,
      bio_context_id,
      lc_config_id,
      data_source_type,
      frag_data_is_raw,
      rt_adjustment_ts,
      consensus_id
    FROM mass_dataset
    WHERE dataset = p_old_dataset;

  DROP TABLE tmp_cluster;

END;
$$;
