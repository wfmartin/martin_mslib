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

--*************************************************************************
--**  MODULE (tables) related to the Global Proteome Machine.
--*************************************************************************
--------------------------------------------------------------------------
CREATE TYPE gpm_hit AS (
  mods_list_str              varchar,
  expect                     real
);


---------------------------------------------------------------------------
CREATE OR REPLACE VIEW gpm_data_info AS 
  SELECT
    dataset,
    csv_filename,
    cef_filename,
    mgf_filename,
    match_source AS gpm_database
  FROM molf_dataset 
  JOIN match_run  USING(dataset);


-----------------------------------------------------------------------------
--  Modifications from the default parameters.
-----------------------------------------------------------------------------
CREATE TABLE gpm_parameter_mods(
  gpm_parameters_key         varchar,
  label                      varchar,
  value                      varchar,

  PRIMARY KEY(gpm_parameters_key, label)
);

COPY gpm_parameter_mods(gpm_parameters_key, label, value) FROM STDIN;
pep_defaults	spectrum, fragment monoisotopic mass error	40
pep_defaults	spectrum, parent monoisotopic mass error plus	20
pep_defaults	spectrum, parent monoisotopic mass error minus	20
pep_defaults	spectrum, parent monoisotopic mass isotope error	yes
pep_defaults	spectrum, fragment monoisotopic mass error units	ppm
pep_defaults	spectrum, maximum parent charge	9
pep_defaults	spectrum, use noise suppression	no
pep_defaults	spectrum, minimum parent m+h	300.0
pep_defaults	spectrum, minimum fragment mz	50.0
pep_defaults	spectrum, minimum peaks	2
pep_defaults	residue, modification mass	
pep_defaults	residue, potential modification mass	+15.994915@M
pep_defaults	protein, cleavage site	[X]|[X]
pep_defaults	refine, spectrum synthesis	yes
pep_defaults	refine, maximum valid expectation value	0.1
pep_defaults	output, maximum valid expectation value	1
pep_defaults	output, maximum valid protein expectation value	1
pep_defaults	protein, quick acetyl	yes
pep_defaults	refine	yes
pep_defaults	refine, potential modification mass	+15.994915@M,+79.966331@S,+79.966331@T,+79.966331@Y
\.


---------------------------------------------------------------------------
CREATE TABLE gpm_run(
  id                         SERIAL PRIMARY KEY,
  dataset                    varchar,
  proteins_file              varchar,
  gpm_parameters_key         varchar,
  decoy_is_random            boolean,
  output_text                varchar,
  decoy_output_text          varchar,
  ts                         timestamp default now()
);

