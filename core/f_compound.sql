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
--
--**  Public COMPOUNDS (functions)
--**************************************************************************

-----------------------------------------------------------------------------
-- Compute the mass of a list of modifications.
-----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION compute_ng_mods_mass(
    p_ng_mods                ng_mod[])
    RETURNS double precision
    LANGUAGE SQL AS $$

  WITH pass_1 AS (SELECT * FROM unnest($1)) 
  SELECT COALESCE(sum(mass), 0)::double precision
  FROM pass_1 JOIN ng_mod_type USING(label);

$$;


---------------------------------------------------------------------------
-- Convert from a list of structs of modifications to string representation.
---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION mods_list_array_to_str(
    p_mods_list              ng_mod[])
    RETURNS varchar
    LANGUAGE SQL AS $$

  SELECT string_agg(label || '@' || position::varchar, ',')
  FROM unnest($1);

$$;


---------------------------------------------------------------------------
-- Convert from string representation to a list of structs of modifications.
---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION mods_list_to_array(
    p_mods_list_str          varchar)
    RETURNS ng_mod[]
    LANGUAGE SQL AS $$

  WITH matched(tokens) AS (SELECT regexp_matches($1, '([a-z]\w+)@(\d+)', 'g'))
  SELECT array_agg((tokens[2]::integer, tokens[1])::ng_mod
                    ORDER BY tokens[2]::integer)
  FROM matched;

$$;


---------------------------------------------------------------------------
--  Shift the positions of the modifications, for translating between
--  positions on the origin protein and positions on the peptides.
---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION mods_list_shift_positions(
    p_mods_list              ng_mod[],
    p_pos_shift              integer)
    RETURNS ng_mod[] 
    IMMUTABLE
    LANGUAGE SQL AS $$
  WITH recs(mods_list) AS (SELECT unnest($1))
  SELECT array_agg( ((mods_list).position+$2, (mods_list).label)::ng_mod )
  FROM recs;
$$;


---------------------------------------------------------------------------
--  Return the database ID (compound_w_origin_id) for a particular
--  peptide compound (with known protein of origin) from table
--  'compound_w_origin'.
--
--  If no such entity exists in the database, create one.  Also create other
--  dependent database entities as required.
---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION get_compound_for_match(
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

