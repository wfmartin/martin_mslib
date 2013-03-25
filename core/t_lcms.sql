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
--  This is a list of compounds (which could potentially be groups of 
--  isomeric compounds that elute similarly) as exported from cumulative 
--  matches generated from technical replicates of a given sample.
--
--  These libraries can also be used to initiate a sample_consensus so that
--  known compounds can be immediately put on the MS/MS exclusion list.
-----------------------------------------------------------------------------
CREATE TABLE lcms_library(
  lcms_library_id            varchar PRIMARY KEY,

  bio_context_id             varchar,
  lc_config_id               varchar,
  
  consensus_ids              varchar[],

  UNIQUE(bio_context_id, lc_config_id) 
  -- It would be too confusing to allow duplicates of this combination.
);


-----------------------------------------------------------------------------
--  This is a library of compounds for a particular LC solvent gradient.
-----------------------------------------------------------------------------
CREATE TABLE lcms_library_compound(
  lcms_library_compound_id   SERIAL PRIMARY KEY,
  lcms_library_id            varchar,
  compound_w_origin_ids      integer[],  -- isomeric compounds
                                         -- NULL for excluded compounds
  source_consensus_cpd_id    integer,

  rt_start                   real,
  rt_end                     real,

  min_z                      integer,
  max_z                      integer,

  mass                       double precision,

  ---------------------------------------------------------------------------
  -- The box type encapsulates a range of both mass (x-axis) and
  -- retention time (y-axis).  There are convenient and efficient operators
  -- for checking overlaps.
  ---------------------------------------------------------------------------
  mass_rt_rectangle          box
);
