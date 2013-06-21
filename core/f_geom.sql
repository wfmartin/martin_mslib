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
-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION encompassing_box(
    p_boxes                  box[])
    RETURNS box
    LANGUAGE SQL AS $$

  SELECT box(
    point( min( least(    (b[0])[0], (b[1])[0] ) ),
           min( least(    (b[0])[1], (b[1])[1] ) ) ),
    point( max( greatest( (b[0])[0], (b[1])[0] ) ),
           max( greatest( (b[0])[1], (b[1])[1] ) ) ) )
  FROM unnest($1) AS b;

$$;


-------------------------------------------------------------------------
--  The proportion of box #1 is overlapped by box #2.
-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION box_overlap_proportion(
    p_box_1                  box,
    p_box_2                  box)
    RETURNS double precision
    LANGUAGE SQL AS $$
  SELECT CASE
    WHEN area($1)>0 AND $1 && $2 THEN  area($1 # $2)/area($1)
    ELSE 0::real
  END;
$$;

