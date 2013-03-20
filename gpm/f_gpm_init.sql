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

-----------------------------------------------------------------------------
--  For each validation, there's a hash keyed by label
--  with any of the following:
--    1) not_empty
--    2) function 
--    3) min_value and/or max_value
--    4) valid_values a hash with values as key
--    5) pattern  (regular expression)
-----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION init_gpm_validations()
    RETURNS void
    LANGUAGE plperl AS $$

  $_SHARED{gpm_parm_validations} = {

    'spectrum, fragment monoisotopic mass error' => {
      not_empty => 1,
    },

    'spectrum, parent monoisotopic mass error plus' => {
      not_empty => 1,
    },

    'spectrum, parent monoisotopic mass error minus' => {
      not_empty => 1,
    },

    'spectrum, parent monoisotopic mass isotope error' => {
      not_empty => 1,
    },

    'spectrum, fragment monoisotopic mass error units' => {
      not_empty => 1,
    },

    'spectrum, maximum parent charge' => {
      not_empty => 1,
    },

    'spectrum, use noise suppression' => {
      not_empty => 1,
    },

    'spectrum, minimum parent m+h' => {
      not_empty => 1,
    },

    'spectrum, minimum fragment mz' => {
      not_empty => 1,
    },

    'spectrum, minimum peaks' => {
      not_empty => 1,
    },

    'residue, modification mass' => {
      #not_empty => 1,
    },

    'residue, potential modification mass' => {
      not_empty => 1,
    },

    'protein, taxon' => {
      not_empty => 1,
    },

    'protein, cleavage site' => {
      not_empty => 1,
    },

    'refine, spectrum synthesis' => {
      not_empty => 1,
    },

    'refine, maximum valid expectation value' => {
      not_empty => 1,
    },
  };
$$;


-----------------------------------------------------------------------------
--  The function init_gpm_validations must be called (within the session)
--  before this function is called.
-----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION validate_gpm_parameter(
    p_label                  varchar,
    p_value                  varchar)
    RETURNS varchar            -- error message (null for OK)
    LANGUAGE plperl AS $$

  my $label = shift;
  local($_) = shift;

  my $validation = $_SHARED{gpm_parm_validations}{$label};

  return "Parameter '$label' is either non-existent or can't be changed."
    unless ($validation);

  if (exists $validation->{not_empty}) {
    return "Parameter '$label' cannot be empty."
      unless (defined and $_ ne '');
  }

  if (exists $validation->{function}) {
    return "Parameter '$label' invalid"
      unless &{$validation->{function}}($_);
  }

  if (exists $validation->{min_value}) {
    return "Parameter '$label' must be greater than or equal to " .
           $validation->{min_value}
      unless ($_ >= $validation->{min_value});
  }

  if (exists $validation->{max_value}) {
    return "Parameter '$label' must be less than or equal to " .
           $validation->{max_value}
      unless ($_ <= $validation->{max_value});
  }

  if (exists $validation->{pattern}) {
    return "Parameter '$label' must match pattern '$validation->{pattern}'"
      unless (/$validation->{pattern}/);
  }

  if (exists $validation->{valid_values}) {
    return "Parameter '$label' must be one of the following values: " .
           join(', ', map("'$_'", keys %{$validation->{valid_values}}))
      unless (exists $validation->{valid_values}{$_});
  }

  return undef;
$$;


-----------------------------------------------------------------------------
--  Validate each of the GPM parameters (in gpm_parameter_mods) table,
--  returning a message for each invalid parameter.
-----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gpm_validate(
    p_gpm_parameters_key     varchar)
    RETURNS varchar          -- error message (null for OK)
    LANGUAGE SQL AS $$

  WITH pass_1 AS (
    SELECT
      validate_gpm_parameter(label, value) AS message
    FROM gpm_parameter_mods
    WHERE gpm_parameters_key = $1
    )
  SELECT message
  FROM pass_1
  WHERE message IS NOT NULL
  LIMIT 1;
$$;


-----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gpm_gen_input_xml(
    p_gpm_parameters_key     varchar,
    p_dataset                varchar,
    p_mgf_input_file         varchar)
    RETURNS varchar
    LANGUAGE plpgsql AS $$
DECLARE
  v_gpm_default_input        varchar;
  v_parms_str                varchar;
BEGIN
  SELECT value INTO v_gpm_default_input
  FROM configuration
  WHERE tag = 'gpm_default_input';

  WITH parms(label,value) AS ( VALUES
      ('list path, default parameters', v_gpm_default_input),
      ('protein, taxon', 'protein_lib'),
      ('output, path', 'output.xml'),
      ('spectrum, path', p_mgf_input_file ),
      ('spectrum, path type', 'mgf')
    ),
  all_parms AS (
    SELECT label, value
    FROM gpm_parameter_mods
    WHERE gpm_parameters_key = p_gpm_parameters_key
    UNION
    SELECT label, value
    FROM parms
    )
  SELECT string_agg(format('  <note type="input" label="%s">%s</note>',
                           label, value), chr(10))
  INTO v_parms_str
  FROM all_parms;

  RETURN '<?xml version="1.0"?>' || chr(10) ||
         '<bioml>' || chr(10) ||
         v_parms_str || chr(10) ||
         '</bioml>' || chr(10);
END;
$$;

