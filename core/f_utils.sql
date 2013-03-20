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

--***********************************************************************
--**  General Utility Functions.
--***********************************************************************

--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION proton_mass()
    RETURNS double precision
    LANGUAGE SQL AS $$
  SELECT 1.00727::double precision;
$$;


--------------------------------------------------------------------------
--  Mass of carbon-13 minus carbon-12, rounded up.
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION effective_neutron_mass()
    RETURNS real
    LANGUAGE SQL AS $$
  SELECT 1.0034::real
$$;


--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION water_mass()
    RETURNS double precision
    LANGUAGE SQL AS $$
  SELECT 18.0105::double precision;
$$;


-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION read_text(p_fname varchar)
    RETURNS varchar
    LANGUAGE plperlu AS $$
  my $fname = shift;

  -f $fname or elog(ERROR, "File $fname doesn't exist.");
  -r $fname or elog(ERROR, "File $fname doesn't have read permission.");

  open(INF, $fname)  or  elog(ERROR, "Can't open $fname");
  local($/) = undef;
  my $text = <INF>;
  close(INF);
  
  $text =~ s/\x00//gs; # remove nulls since they would cause function to crash

  $text;
$$;


-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION write_text(
    p_fname varchar,
    p_text  varchar)
    RETURNS void
    LANGUAGE plperlu AS $$
  my ($fname, $text) = @_;
  open(OUTF, ">$fname")  or  elog(ERROR, "Can't create $fname");
  print OUTF $text;
  close(OUTF);
$$;


-----------------------------------------------------------------------
--  Returns lines of a file with CR/LF removed.
-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION read_lines(p_fname varchar)
    RETURNS SETOF varchar
    LANGUAGE plperlu AS $$

  my $fname = shift;
  open(INF, $fname)  or  elog(ERROR, "Can't open $fname");
  while (<INF>) {
    s/[\r\n]//g;
    return_next $_;
  }
  close(INF);

  [];
$$;


-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION read_lines_to_array(p_fname varchar)
    RETURNS varchar[]
    LANGUAGE SQL AS $$

  WITH x(line) AS (SELECT read_lines($1))
  SELECT array_agg(line) FROM x;
$$;


-----------------------------------------------------------------------
--  Return the mass for an amino acid sequence.
-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION compute_peptide_mass(p_seq varchar)
    RETURNS double precision
    LANGUAGE plperl AS $$

  my ($seq) = @_;
  my $mass = 0.0;

  my %aa_masses = (
   A => 71.037114,
   R => 156.101111,
   N => 114.042927,
   D => 115.026943,
   C => 103.009185,
   E => 129.042593,
   Q => 128.058578,
   G => 57.021464,
   H => 137.058912,
   I => 113.084064,
   L => 113.084064,
   K => 128.094963,
   M => 131.040485,
   F => 147.068414,
   P => 97.052764,
   S => 87.032028,
   T => 101.047679,
   W => 186.079313,
   Y => 163.06332,
   V => 99.068414,
  );

  foreach (0..length($seq)-1) {
    $mass += $aa_masses{substr($seq, $_,1)};
  }

  $mass;
$$;


-----------------------------------------------------------------------
--  Given two pairs of coordinates, return whether the intervals specified
--  by those two pairs of coordinates overlap.
-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION ranges_overlap(
    w1_start                 anyelement,
    w1_end                   anyelement,
    w2_start                 anyelement,
    w2_end                   anyelement)
    RETURNS boolean
    LANGUAGE SQL AS $$

  SELECT $1 BETWEEN $3 AND $4  OR
         $2 BETWEEN $3 AND $4  OR
         $3 BETWEEN $1 AND $2  OR
         $4 BETWEEN $1 AND $2;

$$;


-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gap_between_ranges(
    w1_start                 anyelement,
    w1_end                   anyelement,
    w2_start                 anyelement,
    w2_end                   anyelement)
    RETURNS anyelement
    LANGUAGE SQL AS $$

  SELECT CASE
    WHEN $1 > $4 THEN $1-$4
    WHEN $3 > $2 THEN $3-$2
    ELSE 0
  END;
$$;


-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION ppm_error(
    p_theoretical_mass       double precision,
    p_observed_mass          double precision)
    RETURNS real
    LANGUAGE SQL AS $$

  SELECT (1e6 * ($2-$1)/$1)::real;
$$;


-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION masses_within_error_limit(
    p_theoretical_mass       double precision,
    p_observed_mass          double precision,
    p_error_ppm_limit        real)
    RETURNS boolean
    LANGUAGE SQL AS $$

  SELECT ABS(ppm_error($1, $2)) <= $3;
$$;


-----------------------------------------------------------------------
--  Aggregate function to combine (concatenate) arrays.
-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION combine_arrays_func(
    p_state                  anyarray,
    p_parm                   anyarray)
    RETURNS anyarray
    LANGUAGE SQL AS $$

  SELECT CASE WHEN COALESCE(array_length($1,1),0) = 0  THEN $2
              ELSE $1 || $2
         END;

$$;

CREATE AGGREGATE combine_arrays(anyarray) (
  SFUNC=combine_arrays_func,
  STYPE=anyarray
);


----------------------------------------------------------------------
--  Return a list with redundant values removed.
--------------------------------------------------------------------
CREATE OR REPLACE FUNCTION array_unique_vals(
    p_array                    anyarray)
    RETURNS anyarray
    LANGUAGE SQL AS $$

   WITH pass_1 AS (SELECT DISTINCT z FROM unnest($1) AS z)
   SELECT array_agg(z ORDER BY z) FROM pass_1

$$;


----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION change_session_variable(
    p_variable_name          varchar,
    p_new_value              varchar)
    RETURNS varchar 
    LANGUAGE plpgsql AS $$
DECLARE
  v_old_value                varchar;
BEGIN
  EXECUTE format('SHOW %s', p_variable_name) INTO v_old_value;
  
  EXECUTE format('SET SESSION %s = %s', p_variable_name, p_new_value);

  RETURN v_old_value;
END;
$$;


----------------------------------------------------------------------
-- Reads from 2 tables (or views):
--  tmp_pair(id_1, id_2)   AND tmp_all_ids(id)
--  (all columns are integers)
--
-- Recursively links together paired relationships to generate clusters.
----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION cluster_ids()
    RETURNS TABLE (
      cluster_id             integer,
      id                     integer
    )
    LANGUAGE plpgsql AS $$
BEGIN

  CREATE TEMP TABLE tmp_seed_id AS 
    SELECT DISTINCT id_1 AS id
    FROM tmp_pair t1
    WHERE NOT EXISTS
      (SELECT * FROM tmp_pair t2
       WHERE t1.id_1 = t2.id_2  AND t2.id_1 > t2.id_2);

--CREATE TEMP TABLE tmp_pair2 AS
--  SELECT id_1, id_2
--  FROM tmp_pair
--  WHERE NOT EXISTS (SELECT s.id FROM tmp_seed_id s WHERE s.id = id_2);
    

  RETURN QUERY
  WITH RECURSIVE recurse(cluster_id, id) AS (
    -------------------------------------------------------------------
    --  Seed clusters with peaks that don't overlap with any peak in
    --  a lower sorting order and assign new cluster ID to each seed.
    -------------------------------------------------------------------
    SELECT s.id AS cluster_id, s.id  FROM tmp_seed_id s
      UNION
    -------------------------------------------------------------------
    --  Pull overlapping peaks into the clusters.
    -------------------------------------------------------------------
    SELECT r.cluster_id, id_2 AS id
    FROM recurse r
    JOIN tmp_pair t ON (r.id = t.id_1)
  ),
  tmp_singleton AS (
    SELECT a.id
    FROM tmp_all_ids a
      EXCEPT
    SELECT r.id
    FROM recurse r
  )
  SELECT min(r.cluster_id) AS cluster_id, r.id
  FROM recurse r GROUP BY r.id
  UNION ALL
  SELECT t.id AS cluster_id, t.id AS id
  FROM tmp_singleton t;

  DROP TABLE tmp_seed_id;

END;
$$;
