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

--**************************************************************************
--**  MATCHES of observed (MS) masses to those of theoretical compounds (tables)
--
--**************************************************************************

-------------------------------------------------------------------------
--  Summary of one analysis run matching observed masses to the masses of
--  theoretical compounds (regardless of how the matches were found).
-------------------------------------------------------------------------
CREATE TABLE match_run(
  match_run_id               SERIAL PRIMARY KEY,
  dataset                    varchar,

  match_method               varchar,
  match_source               varchar,
  match_params               varchar,

  num_matches                integer,
  num_decoy_matches          integer  DEFAULT 0,

  score_method               varchar,

  error_ppm_threshold        real,

  date_ts                    timestamp default now()
);


-------------------------------------------------------------------------
--  One match of an observed mass to the mass of a theoretical compound.
-------------------------------------------------------------------------
CREATE SEQUENCE match_id_seq;
CREATE TABLE match(
  id                         integer
                               PRIMARY KEY default nextval('match_id_seq'),
  match_run_id               integer  NOT NULL,
  obs_mass_id                integer  NOT NULL,
  compound_w_origin_id       integer  NOT NULL,
  error_ppm                  real,
  score                      real,
  spectrum_id                integer,
  is_decoy                   boolean DEFAULT False,

  misc                       varchar
);


-------------------------------------------------------------------------
--  Trigger to fill in the error_ppm if not already assigned.
-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION match_insert()
    RETURNS TRIGGER AS $$
DECLARE
  v_observed_mass            double precision;
  v_compound_mass            double precision;
BEGIN
  IF NEW.error_ppm IS NULL THEN
    SELECT mass INTO v_observed_mass 
    FROM observed_mass
    WHERE id = NEW.obs_mass_id;

    SELECT mass INTO v_compound_mass 
    FROM compound_w_origin
    WHERE compound_w_origin_id = NEW.compound_w_origin_id;

    NEW.error_ppm := 
        (1e6 * (v_observed_mass - v_compound_mass)/v_compound_mass)::real;
  END IF;

  RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER match_insert_trigger
  BEFORE INSERT ON match
  FOR EACH ROW EXECUTE PROCEDURE match_insert();


CREATE TABLE deleted_bad_match (LIKE match EXCLUDING ALL);

-------------------------------------------------------------------------
--  One match of an observed fragment mass to the mass of a theoretical
--  fragment compound.
-------------------------------------------------------------------------
CREATE TABLE fragment_match(
  id                         SERIAL PRIMARY KEY,
  match_run_id               integer,
  fragment_mass_id           integer,
  compound_w_origin_id       integer,
  error_ppm                  real,

  UNIQUE(match_run_id, fragment_mass_id, compound_w_origin_id)
);

-------------------------------------------------------------------------
--  Trigger to fill in the error_ppm if not already assigned.
-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION fragment_match_insert()
    RETURNS TRIGGER AS $$
DECLARE
  v_fragment_mass            double precision;
  v_compound_mass            double precision;
BEGIN
  IF NEW.error_ppm IS NULL THEN
    SELECT mass INTO v_fragment_mass 
    FROM fragment_mass
    WHERE fragment_mass_id = NEW.fragment_mass_id;

    SELECT mass INTO v_compound_mass 
    FROM compound_w_origin
    WHERE compound_w_origin_id = NEW.compound_w_origin_id;

    NEW.error_ppm := 
        (1e6 * (v_fragment_mass - v_compound_mass)/v_compound_mass)::real;
  END IF;

  RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER fragment_match_insert_trigger
  BEFORE INSERT ON fragment_match
  FOR EACH ROW EXECUTE PROCEDURE fragment_match_insert();
