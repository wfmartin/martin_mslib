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
--  Return the database ID (compound_w_origin_id) for a particular
--  peptide compound (with known protein of origin) from table
--  'compound_w_origin'.
--
--  If no such entity exists in the database, create one.  Also create other
--  dependent database entities as required.
---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gpm_get_compound(
    p_protein_name           varchar,
    p_aa_seq                 varchar,
    p_start_pos              integer,
    p_mods_list_str          varchar)
    RETURNS integer
    LANGUAGE plpgsql AS $$
DECLARE
  v_peptide_id               integer;
  v_compound_id              integer;
  v_compound_w_origin_id     integer;
BEGIN
  ----------------------------------------------------------------------
  --  Read the peptide from the database.  If it's not there, insert it.
  --
  --  Then read the corresponding compound (with no glycan) without
  --  information about the protein of origin from the database.
  --  If it's not there, insert it.
  ----------------------------------------------------------------------
  SELECT peptide_id INTO v_peptide_id
  FROM peptide
  WHERE aa_seq = p_aa_seq AND ng_mods_str = p_mods_list_str;
  
  IF NOT FOUND THEN
    INSERT INTO peptide(aa_seq, ng_mods_str)
      VALUES(p_aa_seq, p_mods_list_str)
    RETURNING peptide_id INTO v_peptide_id;
  END IF;
  
  SELECT compound_id INTO v_compound_id
  FROM compound
  WHERE peptide_id = v_peptide_id AND glycan_id = 0;
  
  IF NOT FOUND THEN
    INSERT INTO compound(peptide_id)
    VALUES(v_peptide_id)
    RETURNING compound_id INTO v_compound_id;
  END IF;
  
  SELECT compound_w_origin_id INTO v_compound_w_origin_id
  FROM compound_w_origin
  WHERE protein_name = p_protein_name      AND
        start_pos_on_protein = p_start_pos AND
        saps = array[]::integer[]              AND
        compound_id = v_compound_id;
  
  IF NOT FOUND THEN
    INSERT INTO compound_w_origin(compound_id, protein_name,
                                  start_pos_on_protein)
    VALUES(v_compound_id, p_protein_name, p_start_pos)
    RETURNING compound_w_origin_id INTO v_compound_w_origin_id;
  END IF;

  RETURN v_compound_w_origin_id;
END;
$$;


--------------------------------------------------------------------------
--  Parse relationship between the internal GPM id and the molf spectrum_id.
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gpm_parse_support_groups(
    p_xml                    xml)
    RETURNS TABLE(
      gpm_id                 integer,
      mgf_title              varchar)
    LANGUAGE SQL AS $$
     
  WITH sup_grps(grp_xml) AS (
    SELECT unnest( XPATH(
      '//group[@type="support" and @label="fragment ion mass spectrum"]', $1))
    )
  SELECT
    (XPATH('//gaml_trace/@id', grp_xml))[1]::varchar::integer  AS gpm_id,
    trim(both ' ' from translate(
        (XPATH('//note/text()', grp_xml))[1]::varchar, chr(13) || chr(10), ''))
      AS spectrum_id
  FROM sup_grps;
$$;


--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gpm_parse_protein_name(
    p_group_label            varchar)
    RETURNS varchar
    LANGUAGE SQL AS $$

  SELECT regexp_replace(regexp_replace($1, '^[^ ]+\|', ''), ' .*$', '');
$$;


--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gpm_get_mods_list(
    p_domain_xml             xml)
    RETURNS varchar 
    LANGUAGE SQL AS $$

  WITH type_table(mod_mass, mod_label) AS (VALUES
    (79.96633::numeric(8,5), 'phosphate'),
    (15.99492::numeric(8,5), 'oxidation'),
    (0.98402::numeric(8,5),  'deamidation'),
    (-17.02655::numeric(8,5),'ammonia_loss'),   -- Glu => pyroglutamate
    (42.01057::numeric(8,5), 'acetate'),
    (57.02146::numeric(8,5), 'carbamidomethyl')
    ),
  aa_nodes(aa_xml) AS (SELECT unnest(XPATH('//aa', $1))),
  aa_data(mod_mass, position) AS (
    SELECT
      (XPATH('./@modified', aa_xml))[1]::varchar::numeric(8,5),
      (XPATH('./@at',       aa_xml))[1]::varchar::integer -
         (XPATH('./@start', $1))[1]::varchar::integer + 1
    FROM aa_nodes
  ),
  mods AS (
    SELECT mod_label, position
    FROM aa_data JOIN type_table USING(mod_mass)
  )
  SELECT string_agg(mod_label || '@' || position, ',')
  FROM mods;
$$;


--------------------------------------------------------------------------
--  Parse info from one peptide group a "group" element.
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gpm_parse_peptide_group(
    p_pep_grp_xml            xml,
    OUT gpm_id               integer,
    OUT protein_name         varchar,
    OUT start_pos            integer,
    OUT end_pos              integer,
    OUT aa_seq               varchar,
    OUT gpm_hits             gpm_hit[])
    LANGUAGE plpgsql AS $$
DECLARE
  v_domains                  xml[];
BEGIN
  gpm_id := (XPATH('/group/@id', p_pep_grp_xml))[1]::varchar::integer;
  protein_name := gpm_parse_protein_name(
      (XPATH('/group/@label', p_pep_grp_xml))[1]::varchar);

  v_domains := XPATH('//domain', p_pep_grp_xml);
  start_pos := (XPATH('./@start', v_domains[1]))[1]::varchar::integer;
  end_pos   := (XPATH('./@end', v_domains[1]))[1]::varchar::integer;
  aa_seq    := (XPATH('./@seq', v_domains[1]))[1]::varchar;

  WITH domain_nodes(d_xml) AS (SELECT unnest(v_domains)),
  domain_data AS (
    SELECT 
      (XPATH('./@expect', d_xml))[1]::varchar::real  AS expect,
      COALESCE(gpm_get_mods_list(d_xml), '')         AS mods_list_str
    FROM domain_nodes)
  SELECT array_agg( (mods_list_str, expect)::gpm_hit )
  INTO gpm_hits
  FROM domain_data;
END;
$$;


--------------------------------------------------------------------------
--  Iterate through all the "group" elements, parsing the peptide info.
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gpm_parse_peptide_groups(
    p_xml                    xml)
    RETURNS TABLE(
      gpm_id                 integer,
      protein_name           varchar,
      start_pos              integer,
      end_pos                integer,
      aa_seq                 varchar,
      gpm_hits               gpm_hit[]
    )
    LANGUAGE SQL AS $$

  WITH pep_grp_nodes(pep_xml) AS (SELECT unnest(XPATH('//group[@id]', $1)))
  SELECT (gpm_parse_peptide_group(pep_xml)).*
  FROM pep_grp_nodes
$$;
  

---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gpm_handle_xml(
    p_gpm_xml                xml,
    p_dataset                varchar)
    RETURNS TABLE(
      obs_mass_id            integer,
      compound_w_origin_id   integer,
      score                  real,
      spectrum_id            integer
    )
    LANGUAGE SQL AS $$
    
  WITH 
    pep  AS (SELECT (gpm_parse_peptide_groups($1)).*),
    spec AS (SELECT (gpm_parse_support_groups($1)).*),
    combined AS (SELECT * FROM pep JOIN spec USING(gpm_id)),
    expanded AS (
      SELECT
        o.id AS obs_mass_id,
        protein_name,
        aa_seq,
        start_pos,
        unnest(gpm_hits) AS gpm_hit,
        cs.spectrum_id
      FROM combined c
      JOIN molf_spectrum s USING(mgf_title)
      JOIN component_spectrum cs  ON 
          (cs.dataset = $2 AND cs.spectrum_id = s.spectrum_id)
      JOIN observed_mass o ON
          (o.dataset = $2 AND o.peak_id = parent_peak_id)
      ),
    by_spectrum AS (
      SELECT
        e.obs_mass_id,
        gpm_get_compound(protein_name, aa_seq, start_pos,
                         (gpm_hit).mods_list_str)
            AS compound_w_origin_id,
        (-log((gpm_hit).expect))::real AS score,
        e.spectrum_id
      FROM expanded e
      )
    --  Choose the best scoring spectrum:
    SELECT
      s.obs_mass_id,
      s.compound_w_origin_id,
      max(s.score) AS score,
      (array_agg(s.spectrum_id ORDER BY s.score DESC))[1] AS spectrum_id
    FROM by_spectrum s
    GROUP BY s.obs_mass_id, s.compound_w_origin_id;
$$;


---------------------------------------------------------------------------
--  Find peptide compound matches using the GPM (Global Proteome Machine).
---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gpm_parse_results(
    p_gpm_output_fname       varchar,
    p_dataset                varchar)
    RETURNS TABLE(
      obs_mass_id            integer,
      compound_w_origin_id   integer,
      score                  real,
      spectrum_id            integer
    )
    LANGUAGE plpgsql AS $$
DECLARE
  v_gpm_output               xml;
BEGIN
  ----------------------------------------------------------------------
  --  Parse the output from the GPM process and store the peptide matches
  --  into a standard representation (in 'match' table).
  ----------------------------------------------------------------------
  v_gpm_output := XMLPARSE(DOCUMENT
      replace( read_text(p_gpm_output_fname), 'GAML:', 'gaml_') );

  RETURN QUERY SELECT * FROM gpm_handle_xml(v_gpm_output, p_dataset);

END;
$$;


---------------------------------------------------------------------------
--  This function takes files with the raw output from the GPM program,
--  parses them, and stores the results as 'match' records in the database.
---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION process_gpm_output(
    p_gpm_parameters_key     varchar,
    p_dataset                varchar,
    p_gpm_outfile            varchar,
    p_gpm_decoys_outfile     varchar,
    p_min_score              real)
    RETURNS integer
    LANGUAGE plpgsql AS $$
DECLARE
  v_match_run_id             integer;
  v_num_matches              integer;
  v_num_decoy_matches        integer;
BEGIN
  INSERT INTO match_run(dataset, match_method, match_params, score_method)
  VALUES(p_dataset, 'GPM', p_gpm_parameters_key, '-log(expect)')
  RETURNING match_run_id INTO v_match_run_id;

  ----------------------------------------------------------------------
  --  Parse the output from the GPM process and store the peptide matches
  --  into a standard representation (in 'match' table).
  ----------------------------------------------------------------------
  INSERT INTO match(match_run_id, obs_mass_id, compound_w_origin_id, score,
                    spectrum_id, is_decoy, misc)
    WITH pass_1 AS (
        SELECT *
        FROM gpm_parse_results(p_gpm_outfile, p_dataset)
        WHERE p_min_score IS NULL OR score >= p_min_score
      )
      SELECT
        v_match_run_id,
        obs_mass_id,
        compound_w_origin_id,
        score,
        spectrum_id,
        False,
        'GPM'
      FROM pass_1;

  GET DIAGNOSTICS v_num_matches := ROW_COUNT;

  v_num_decoy_matches := 0;

  IF p_gpm_decoys_outfile  IS NOT NULL  THEN
    INSERT INTO match(match_run_id, obs_mass_id, compound_w_origin_id, score,
                      spectrum_id, is_decoy, misc)
      WITH pass_1 AS (
          SELECT *
          FROM gpm_parse_results(p_gpm_decoys_outfile, p_dataset)
          WHERE p_min_score IS NULL OR score >= p_min_score
          )
        SELECT
          v_match_run_id,
          obs_mass_id,
          compound_w_origin_id,
          score,
          spectrum_id,
          True,
          'GPM'
        FROM pass_1;

    GET DIAGNOSTICS v_num_decoy_matches := ROW_COUNT;
  END IF;

  UPDATE match_run
  SET num_matches       = v_num_matches,
      num_decoy_matches = v_num_decoy_matches
  WHERE match_run_id = v_match_run_id;

  RETURN v_match_run_id;
END;
$$;

