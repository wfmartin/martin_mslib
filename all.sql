SET SESSION client_min_messages = WARNING;
CREATE LANGUAGE plperl;
CREATE LANGUAGE plperlu;

\i core/t_utils.sql
\i core/f_utils.sql
\i core/f_geom.sql
\i core/f_aa.sql
\i core/t_compound.sql
\i core/f_compound.sql
\i core/t_mass_data.sql
\i core/t_match.sql
\i core/f_match.sql
\i core/t_lcms.sql
\i core/f_lcms.sql
\i core/t_consensus.sql
\i core/f_consensus.sql

\i core/f_mass.sql

\i agilent/t_molf.sql
\i agilent/f_molf.sql

\i gpm/t_gpm.sql
\i gpm/f_gpm_init.sql
\i gpm/f_gpm_parse.sql


SET SESSION client_min_messages = INFO;
