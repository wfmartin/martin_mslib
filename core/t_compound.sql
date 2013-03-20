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

--************************************************************************
--**  Public COMPOUNDS  (tables)
--************************************************************************

-----------------------------------------------------------------------------
CREATE TABLE ng_mod_type(
  label                      varchar PRIMARY KEY,
  mass                       double precision
);

COPY ng_mod_type FROM STDIN WITH delimiter ',';
phosphate,79.966331
oxidation,15.994915
deamidation,0.984016
ammonia_loss,-17.02655
acetate,42.01057
carbamidomethyl,57.021464
\.


-----------------------------------------------------------------------------
-- Non-glycan modification of a peptide.
-----------------------------------------------------------------------------
CREATE TYPE ng_mod AS (
  position                   integer,
  label                      varchar
);


-----------------------------------------------------------------
CREATE TABLE sap(
  sap_id                     SERIAL PRIMARY KEY,

  protein_name               varchar,
  position                   integer,
  ref_aa                     char,
  alt_aa                     char,
  dbsnp_id                   varchar,
  mass_adjustment            double precision,

  UNIQUE(protein_name, position, ref_aa)
);


-----------------------------------------------------------------------------
--  Global library of peptides.
--  Peptide residues (without context information,
--                    what larger protein they come from).
-----------------------------------------------------------------------------
CREATE TABLE peptide(
  peptide_id                 SERIAL PRIMARY KEY,
  aa_seq                     varchar NOT NULL,   -- always uppercase
  ng_mods_str                varchar   DEFAULT '',
  ng_mods                    ng_mod[]  DEFAULT array[]::ng_mod[],

  seq_len                    integer,

  mass                       double precision,

  UNIQUE (aa_seq, ng_mods_str)
);



-----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION peptide_insert()
    RETURNS TRIGGER
    LANGUAGE plpgsql AS $$
BEGIN
  IF NEW.aa_seq IS NULL THEN
    RAISE EXCEPTION 'peptide.aa_seq cannot be null.';
  END IF;

  IF NEW.ng_mods IS NULL THEN
    NEW.ng_mods := ARRAY[]::ng_mod[];
  END IF;

  IF NEW.ng_mods_str IS NULL THEN
    NEW.ng_mods_str := '';
  END IF;

  -------------------------------------------------------------------------
  -- If modification information is provided either way:
  -------------------------------------------------------------------------
  IF char_length(NEW.ng_mods_str) > 0  OR 
     array_length(NEW.ng_mods,1) IS NOT NULL  THEN
  
    IF char_length(NEW.ng_mods_str) = 0 THEN
      NEW.ng_mods_str := mods_list_array_to_str(NEW.ng_mods);

    ELSIF array_length(NEW.ng_mods,1) IS NULL THEN
      NEW.ng_mods := mods_list_to_array(NEW.ng_mods_str);
    END IF;
  END IF;

  IF NEW.seq_len IS NULL THEN
    NEW.seq_len := char_length(NEW.aa_seq);
  END IF;

  IF NEW.mass IS NULL THEN
    NEW.mass := compute_peptide_mass(NEW.aa_seq) +
                compute_ng_mods_mass(NEW.ng_mods);
  END IF;

  RETURN NEW;
END;
$$;

CREATE TRIGGER peptide_insert_trigger
  BEFORE INSERT ON peptide
  FOR EACH ROW EXECUTE PROCEDURE peptide_insert();


-------------------------------------------------------------------------------
----
-------------------------------------------------------------------------------
--CREATE TYPE glycan_mod AS (
--  placeholder integer;
--);
--
--
-------------------------------------------------------------------------------
----
-------------------------------------------------------------------------------
--CREATE TABLE glycan(
--  glycan_id                  SERIAL PRIMARY KEY,
--
--  composition                integer,
--  g_mods                     glycan_mod[],
--
--  taxonomy_id                integer,
--
--  mass                       double precision
--);


-----------------------------------------------------------------------------
--  Context-free compound; that is to say that the information about which
--  protein it's part of is not included here.
-----------------------------------------------------------------------------
CREATE TABLE compound(
  compound_id                SERIAL PRIMARY KEY,
  peptide_id                 integer DEFAULT 0, 
  glycan_id                  integer DEFAULT 0,
  glycan_position            integer,

  mass                       numeric(10,5)
);



-----------------------------------------------------------------
CREATE OR REPLACE FUNCTION compound_insert()
    RETURNS TRIGGER
    LANGUAGE plpgsql AS $$
DECLARE
  v_peptide_mass             double precision;
  v_glycan_mass              double precision;
BEGIN
  IF NEW.glycan_id IS NULL THEN
    NEW.glycan_id := 0;
  END IF;

  IF NEW.peptide_id IS NULL THEN
    NEW.peptide_id := 0;
  END IF;

  IF NEW.mass IS NULL THEN
    IF NEW.peptide_id IS NOT NULL THEN
      SELECT mass INTO v_peptide_mass
      FROM peptide
      WHERE peptide_id = NEW.peptide_id;
    ELSE
      v_peptide_mass := 0.0;
    END IF;

    -- ***************  NEED CODE FOR GLYCANS  *******************
    v_glycan_mass := 0.0;

    NEW.mass = v_peptide_mass + v_glycan_mass + water_mass();
  END IF;

  RETURN NEW;
END;
$$;

CREATE TRIGGER compound_insert_trigger
  BEFORE INSERT ON compound
  FOR EACH ROW EXECUTE PROCEDURE compound_insert();


-----------------------------------------------------------------------------
--  One compound from a particular protein origin.
-----------------------------------------------------------------------------
CREATE TABLE compound_w_origin(
  compound_w_origin_id       SERIAL PRIMARY KEY,
  
  compound_id                integer  NOT NULL,
  protein_name               varchar  NOT NULL,
  start_pos_on_protein       integer  NOT NULL,

  saps                       integer[]  DEFAULT array[]::integer[],

  mass                       numeric(10,5)
);

CREATE OR REPLACE FUNCTION compound_w_origin_insert()
    RETURNS TRIGGER
    LANGUAGE plpgsql AS $$
DECLARE
  v_mass                     double precision;
BEGIN
  IF NEW.mass IS NULL THEN
    SELECT mass INTO v_mass
    FROM compound
    WHERE compound_id = NEW.compound_id;

    NEW.mass = v_mass;
  END IF;

  RETURN NEW;
END;
$$;

CREATE TRIGGER compound_w_origin_insert_trigger
  BEFORE INSERT ON compound_w_origin
  FOR EACH ROW EXECUTE PROCEDURE compound_w_origin_insert();

