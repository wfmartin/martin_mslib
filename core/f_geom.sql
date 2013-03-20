-- -------------------------------------------------------------------------
-- CREATE OR REPLACE FUNCTION extend_box(
--     p_box                    box,
--     p_point                  point)
--     RETURNS box
--     LANGUAGE plpgsql AS $$
-- DECLARE
--   v_box                      box;
-- BEGIN
--   v_box := p_box;
-- 
--   --  Shift to X coordinate:
-- 
--   --  First point of box is (relatively) pointing in the direction of shift:
--   IF sign(p_point[0]) = sign((v_box[0])[0] - (v_box[1])[0])  THEN
--     v_box[0] := point( (v_box[0])[0] + p_point[0], (v_box[0])[1] );
--     v_box := box( point((v_box[0])[0], (v_box[0])[0] + p_point[0]), v_box[1]);
-- 
--   ELSIF sign(p_point[0]) <> 0 THEN
--     v_box[1] := point( (v_box[1])[0] + p_point[0], (v_box[1])[1] );
--     v_box := box( v_box[0], point(((v_box[1])[1] + p_point[1], v_box[1])[0]) );
--   END IF;
-- 
-- 
--   --  Shift to Y coordinate:
-- 
--   --  First point of box is (relatively) pointing in the direction of shift:
--   IF sign(p_point[1]) = sign((v_box[0])[1] - (v_box[1])[1])  THEN
--     v_box := box( point((v_box[0])[0], (v_box[0])[1] + p_point[1]), v_box[1] );
-- 
--   ELSIF sign(p_point[0]) <> 0 THEN
--     v_box := box( v_box[0], point((v_box[1])[0], (v_box[1])[1] + p_point[1]) );
--   END IF;
-- 
--   RETURN v_box;
-- 
-- END;
-- $$;


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

