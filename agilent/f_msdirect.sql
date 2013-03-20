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

---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION agilent_ms_directions(
    p_targeted_mode          boolean,
    p_mz_width_string        varchar,
    p_cmg_id                 varchar,
    p_num_preferred_masses   integer,
    p_preferred_mass_offset  integer,
    p_dataset                varchar, -- optional, must match later loaded data
    p_max_num_neutrons       integer,  -- used only for exclusion lists
    p_mz_ppm_margin          real,
    p_rt_margin              real,
    p_min_mz                 real,
    p_max_mz                 real,
    p_pref_min_rt            real,  -- only effects preferred list
    p_pref_max_rt            real,
    p_exclude_all_fragmented boolean,
    p_max_mz_filter          real)  -- remove all lines with mz above this value
    RETURNS SETOF varchar
    LANGUAGE plpgsql AS $$
DECLARE
  v_mz_width_string          varchar;
  v_theoretical_max_rt       real;
BEGIN
  CREATE TEMP TABLE tmp_msd_items AS 
    SELECT *
    FROM get_ms_direct_items($3, NOT p_targeted_mode,
                             $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14);

  -----------------------------------------------------------------------
  --  Create one wide exclusion entry (times number of charges) for each of
  --  minimum and/or maximum mz parameters.
  --
  --  This is instead of all the individual mz values for masses that 
  --  would otherwise be in the list.
  -----------------------------------------------------------------------
  IF p_min_mz IS NOT NULL  OR  p_max_mz IS NOT NULL  THEN
    SELECT theoretical_max_rt INTO v_theoretical_max_rt
    FROM cum_matches_group JOIN lc_configuration USING(lc_config_id)
    WHERE cmg_id = p_cmg_id;

    CREATE TEMP TABLE tmp_pass_1(
       mz                  real,
       delta_mz_ppm        real
    );

    IF p_min_mz IS NOT NULL THEN
      INSERT INTO tmp_pass_1(mz, delta_mz_ppm)
      VALUES(p_min_mz/2.0::real, 500000::real);
    END IF;

    -- Use 6000 mz as the top of mz block
    IF p_max_mz IS NOT NULL THEN
      INSERT INTO tmp_pass_1(mz, delta_mz_ppm)
      VALUES((6000::real + p_min_mz)/2.0::real, 500000::real);
    END IF;

    INSERT INTO tmp_msd_items(is_preferred, charge, mz, delta_mz_ppm,
                              rt, delta_rt)
      WITH tmp_charge(charge) AS (SELECT generate_series(1,8))
      SELECT
        False,
        charge,
        mz,
        delta_mz_ppm,
        v_theoretical_max_rt/(2::real),
        v_theoretical_max_rt/(2::real)
      FROM tmp_pass_1 CROSS JOIN tmp_charge;

    DROP TABLE tmp_pass_1;

  END IF;

  ---------------------------------------------------------------------
  --  With targeted mode, only preferred peaks are used.
  --  With Auto MS/MS mode, excluded peaks are always used and 
  --  preferred peaks are optionally used.
  ---------------------------------------------------------------------

  IF p_targeted_mode THEN
    RETURN NEXT 'TargetedMSMSTable,,,,,,,';
    RETURN NEXT 'On,Prec. m/z,Z,Ret. Time (min),Delta Ret. Time (min),' ||
                'Iso. Width,Collision Energy,Acquisition Time (ms/spec)';
  ELSE
    RETURN NEXT 'AutoPreferredExcludeMSMSTable,,,,,,,,';
    RETURN NEXT 'On,Prec. m/z,Delta m/z (ppm),Z,Prec. Type,Ret. Time (min),' ||
                'Delta Ret. Time (min),Iso. Width,Collision Energy';
  END IF;

  v_mz_width_string = COALESCE(v_mz_width_string, 'Narrow (~1.3 m/z)');

  IF p_targeted_mode THEN
    RETURN QUERY
    SELECT concat_ws(',', 'True',
      mz::numeric(9,4),
      charge,
      rt,
      delta_rt,
      v_mz_width_string,
      '',
      ''
      )::varchar
    FROM tmp_msd_items
    WHERE p_max_mz_filter IS NULL OR mz < p_max_mz_filter;
  ELSE
    RETURN QUERY
    SELECT concat_ws(',', 'True',
      mz::numeric(9,4),
      delta_mz_ppm::numeric(9,4),
      charge,
      CASE WHEN is_preferred THEN 'Preferred' ELSE 'Exclude' END,
      rt,
      delta_rt,
      v_mz_width_string,
      ''
      )::varchar
    FROM tmp_msd_items
    WHERE p_max_mz_filter IS NULL OR mz < p_max_mz_filter;

  END IF;

  DROP TABLE tmp_msd_items;
END;
$$;
