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

--------------------------------------------------------------------------
--  Parse the ".csv" file produced by "Find compounds by Molecular Feature".
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION parse_molf_csv(
    p_molf_csv_file          varchar)
    RETURNS TABLE(
      cpd                        integer,
      base_peak                  numeric(10,5),
      rt                         real,
      rt_start                   real,
      rt_end                     real,
      rt_width_at_half_ht        real,
    
      min_z                      integer,
      max_z                      integer,
      ions                       integer,
      saturated                  varchar,
      ms_ms_count                integer,
    
      mz                         numeric(10,4),
      mass                       numeric(10,4),
      volume                     real,
      volume_pct                 real,
      height                     real
    )
    STABLE
    LANGUAGE plperlu AS $$

  my ($fname) = @_;

  local $/ = undef;
  open(INF, $fname)  or  elog(ERROR, "Can't open file $fname");
  my $text = <INF>;
  close(INF);

  my (@lines) =  split /[\r\n]+/, $text;
  shift @lines   if ($lines[0] =~ /^Version/);
  shift @lines   if (length($lines[0]) == 0);

  #s/,"[^"]+,[^"]+",/,unused,/g  foreach @lines;
  #s/"//g   foreach @lines;    # " (match it)

  my $header = shift @lines;
  $header =~ s/"(Diff[_ ]\(\w+),\s*(\w+\))"/\1;\2/g;
                              # remove commas in hdr names
  my @in_cols = split /,/, $header;
  my $space_or_underscore = ($header =~ /Min Z/) ? ' ' : '_';

  my %map_cols = (
    'cpd'            => "Cpd",
    'base_peak'      => "Base${space_or_underscore}Peak",
    'ms_ms_count'    => "MS/MS${space_or_underscore}Count",
    'mz'             => "m/z",
    'rt_end'         => "End",
    'ions'           => "Ions",
    'height'         => "Height",
    'max_z'          => "Max${space_or_underscore}Z",
    'min_z'          => "Min${space_or_underscore}Z",
    'mass'           => "Mass",
    'rt'             => "RT",
    'rt_width_at_half_ht' => "Width",
    'rt_start'       => "Start",
    'saturated'      => "Saturated",
    'volume'         => "Vol",
    'volume_pct'     => "Vol \%",
  );

  my ($key, $val);
  foreach (@lines) {
    my $rec = {};
    @{$rec}{@in_cols} = split /,/;

    my $out_rec = { };
    while ( ($key, $val) = each %map_cols) {
      $out_rec->{$key} = $rec->{$val};
    }
    $out_rec->{ms_ms_count} = 0  if ($out_rec->{ms_ms_count} == '');
    return_next($out_rec);
  }

  return [];
$$;


--------------------------------------------------------------------------
--
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION parse_molf_cef(
    p_molf_cef_fname         varchar)
    RETURNS TABLE(
      csv_mass               numeric(10,4),
      csv_rt                 real,
      ce                     real,
      component_mz           numeric(10,4),
      component_rt           real,
      component_charge       integer,
      component_num_neutrons integer,
      component_recipe       varchar,
      component_volume       real
    )
    LANGUAGE plpgsql AS $$
DECLARE
  v_text                     varchar;
  v_root_xml                 xml;
  v_compound                 xml;
  v_location                 xml;
  v_peak                     xml;
  v_spec                     xml;
  v_ce                       varchar;
BEGIN
  ------------------------------------------------------------------------
  --  For some unknown reason the file has 3 unprintable characters at the
  --  start of it which would upset the parser if not removed.
  --  Also, carriage returns (legacy of MS-DOS) have to be removed.
  ------------------------------------------------------------------------
  v_text :=
    regexp_replace(
      regexp_replace(
         replace(read_text(p_molf_cef_fname), chr(13), ''),
      '^[\177-\377]+', ''),
    '<XYData>.*?</XYData>', '', 'g');

  IF v_text NOT LIKE '<CEF%' THEN
    RAISE EXCEPTION 'CEF file (%) not of expected format.', p_molf_cef_fname;
  END IF;

  v_root_xml := XMLPARSE( DOCUMENT v_text );

  FOR v_compound IN (SELECT unnest( XPATH('//Compound', v_root_xml)) ) LOOP
    v_spec := ( XPATH('//Spectrum[@type="TOF-MS2"]', v_compound) )[1];
    v_ce := ( XPATH('//MSDetails/@ce', v_spec) )[1]::varchar;
    ce := substring(v_ce FROM 1 FOR char_length(v_ce)-1);

    v_location := ( XPATH('//Location', v_compound ) )[1];
    csv_mass := ( XPATH('@m', v_location ) )[1];
    csv_rt := ( XPATH('@rt', v_location ) )[1];

    FOR v_peak IN (SELECT unnest( XPATH('//p', v_compound) ) )  LOOP
      component_mz     := ( XPATH('@x',  v_peak ) )[1];
      component_rt     := ( XPATH('@rt', v_peak ) )[1];
      component_recipe := ( XPATH('@s',  v_peak ) )[1];
      component_charge := ( XPATH('@z',  v_peak ) )[1];
      component_volume := ( XPATH('@v',  v_peak ) )[1];

      component_num_neutrons := COALESCE(
        (WITH pass_1(match) AS
           (SELECT regexp_matches(component_recipe, '(\d)$'))
         SELECT match[1]::integer FROM pass_1),
          0);

      RETURN NEXT;
    END LOOP;
  END LOOP;
END;
$$;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION parse_molf_mgf(
    p_molf_mgf_fname         varchar)
    RETURNS TABLE(
      pepmass                double precision,
      rtinseconds            real,
      charge                 integer,
      mgf_index              integer,
      title                  varchar,
      peaks                  molf_peak[]
    )
    LANGUAGE plperlu AS $$

  my ($fname) = @_;

  $/ = undef;
  open(INF, $fname)  or  elog(ERROR, "Can't open file $fname");
  my $text = <INF>;
  close(INF);

  elog(ERROR, "MGF file ($fname) not of expected format.")
      unless ($text =~ /^BEGIN IONS/);

  my ($ignore, @sections) = split /BEGIN IONS/s, $text;

  my $peaks_str;

  my @pattern_lines = (
      'PEPMASS=(\d+\.\d+)',
      'CHARGE=(\d{1,2}[-\+])',
      'TITLE=([^\n\r]+)',
      'RTINSECONDS=(\d+(?:\.\d+)?)',
      '(.*)\z'
      );

  my $pattern_w_charge = join('\s+', @pattern_lines);
  my $pattern_wo_charge = join('\s+', @pattern_lines[0,2,3,4]);

  my $mgf_index = 0;
  my @rows = ();
  foreach my $section (@sections) {
    $section =~ s/END IONS\s*$//s;
    my $rec = { mgf_index => ++$mgf_index };

    if ($section =~ /CHARGE=/) {
      (@{$rec}{qw(pepmass charge title rtinseconds)}, $peaks_str) =
          ($section =~ /$pattern_w_charge/s);
      $rec->{charge} = substr($rec->{charge}, 1,1) . substr($rec->{charge},0,1);
    }
    else {
      (@{$rec}{qw(pepmass title rtinseconds)}, $peaks_str) =
          ($section =~ /$pattern_wo_charge/s);
      $rec->{charge} = undef;
    }

    my @peaks = map {
      s/^\s+//;
      s/\s+\z//;
      my $peak = {};
      @{$peak}{qw(mass ion_count)} = split /\t/;
      $peak;
    }  split /[\r\n]+/, $peaks_str;

    $rec->{peaks} = [ @peaks ];
    return_next $rec;
  }

  [];
$$;


-----------------------------------------------------------------------
--  Input temp tables:  tmp_cef, tmp_molf_rec, tmp_mgf
-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION combine_molf_data(
    p_charge_carrier_mass    double precision)
    RETURNS TABLE(
      cpd                    integer,
      spectrum_id            integer,
      csv_mass               double precision,
      component_charge       integer,
      component_recipe       varchar,
      component_volume       real,
      min_z                  integer,
      max_z                  integer,
      mz                     double precision,
      mgf_rt                 real,
      csv_rt                 real,
      rt_width_at_half_ht    real,
      rt_start               real,
      rt_end                 real,
      num_neutrons           integer,
      volume                 real,
      volume_pct             real,
      height                 real,
      peaks                  molf_peak[],
      collision_energy       real,
      mgf_index              integer,
      mgf_title              varchar
    )
    LANGUAGE plpgsql AS $$
BEGIN
  RETURN QUERY 
  WITH pass_1 AS (
    SELECT
      m.cpd,
      COALESCE(m.csv_mass,
              (tmp_mgf.pepmass - p_charge_carrier_mass) * tmp_mgf.charge)
          AS csv_mass,
      COALESCE(m.component_charge, tmp_mgf.charge) 
          AS component_charge,
      m.component_recipe,
      m.component_volume,
      m.min_z,
      m.max_z,
      COALESCE(m.component_mz, tmp_mgf.pepmass) AS mz,
      (tmp_mgf.rtinseconds/60.0)::real AS mgf_rt,
      m.csv_rt,
      m.rt_width_at_half_ht,
      m.rt_start,
      m.rt_end,
      COALESCE(m.component_num_neutrons,0) AS num_neutrons,
      m.volume,
      m.volume_pct,
      m.height,
      tmp_mgf.peaks,
      m.collision_energy,
      tmp_mgf.mgf_index,
      TRIM(both ' ' from tmp_mgf.title)::varchar AS mgf_title
    FROM tmp_molf_rec AS m
    FULL OUTER JOIN tmp_mgf ON (
       ABS(tmp_mgf.pepmass - m.component_mz)/tmp_mgf.pepmass < 3e-5
         AND
       tmp_mgf.charge = m.component_charge  AND
       tmp_mgf.rtinseconds/60.0 BETWEEN m.rt_start AND m.rt_end+0.2)
    )
  SELECT
    p.cpd,
    nextval('molf_spectrum_seq')::integer AS spectrum_id,
    p.csv_mass,
    p.component_charge,
    p.component_recipe,
    p.component_volume,
    p.min_z,
    p.max_z,
    p.mz,
    p.mgf_rt,
    COALESCE(p.csv_rt, p.mgf_rt)    AS csv_rt,
    p.rt_width_at_half_ht,
    COALESCE(p.rt_start, p.mgf_rt)  AS rt_start,
    COALESCE(p.rt_end, p.mgf_rt)    AS rt_end,
    p.num_neutrons,
    p.volume,
    p.volume_pct,
    p.height,
    p.peaks,
    p.collision_energy,
    p.mgf_index,
    p.mgf_title
  FROM pass_1 p;

END;
$$;


-----------------------------------------------------------------------
CREATE OR REPLACE FUNCTION load_molf_data(
    p_dataset                varchar,
    p_lc_config_id           varchar,
    p_molf_csv_file          varchar,
    p_molf_cef_file          varchar,
    p_molf_mgf_file          varchar)
    RETURNS void AS $$
DECLARE
  v_cef_count                integer;
  v_molf_count               integer;
  v_mgf_count                integer;
  v_net_mgf_count            integer;
  v_charge_carrier_mass      double precision;
BEGIN

  IF p_molf_csv_file IS NULL  OR p_molf_cef_file IS NULL  THEN
    RAISE EXCEPTION 'load_molf_data: You must provide both CSV and CEF files.';
  END IF;

  SELECT charge_carrier_mass INTO v_charge_carrier_mass
  FROM lc_configuration
  WHERE lc_config_id = p_lc_config_id;

  IF NOT FOUND THEN
    RAISE EXCEPTION 'Unknown lc_config_id (%)', lc_config_id;
  END IF;

  INSERT INTO molf_dataset(dataset, lc_config_id,
      csv_filename, cef_filename, mgf_filename)
  VALUES(p_dataset, p_lc_config_id, 
         p_molf_csv_file, p_molf_cef_file, p_molf_mgf_file);

  CREATE TEMP TABLE tmp_cef AS 
    SELECT * FROM parse_molf_cef(p_molf_cef_file);

  GET DIAGNOSTICS v_cef_count := ROW_COUNT;
  
  -----------------------------------------------------------------------
  -- Read molf CSV and CEF files.
  -----------------------------------------------------------------------
  CREATE TEMP TABLE tmp_molf_rec AS
    WITH
      csv AS (SELECT * FROM parse_molf_csv(p_molf_csv_file))
    SELECT
      p_dataset AS dataset,
      cpd,
      csv_mass,
      component_charge,
      component_recipe,
      component_volume,
      component_mz,
      min_z,
      max_z,
      csv_rt,
      rt_start,
      rt_end,
      rt_width_at_half_ht,
      component_num_neutrons,
      volume,
      volume_pct,
      height,
      ce AS collision_energy
    FROM tmp_cef
    JOIN csv ON (csv_mass = csv.mass AND csv_rt = csv.rt);

  GET DIAGNOSTICS v_molf_count := ROW_COUNT;

  IF v_cef_count <> v_molf_count  THEN
    RAISE EXCEPTION 'The CSV and CEF files do not synch up (%,%).',
        p_molf_csv_file, p_molf_cef_file;
  END IF;
      
  CREATE TEMP SEQUENCE molf_spectrum_seq;

  -----------------------------------------------------------------------
  --  When a tandem data file is provided:
  -----------------------------------------------------------------------
  IF p_molf_mgf_file IS NOT NULL THEN

    CREATE TEMP TABLE tmp_mgf AS
      SELECT * FROM parse_molf_mgf(p_molf_mgf_file);

    GET DIAGNOSTICS v_mgf_count := ROW_COUNT;

    ---------------------------------------------------------------------
    CREATE TEMP TABLE tmp_combined AS 
      SELECT * FROM combine_molf_data(v_charge_carrier_mass);

    ---------------------------------------------------------------------
    -- Cluster the MGF recs that don't have a corresponding CSV rec.
    ---------------------------------------------------------------------
    CREATE TEMP TABLE tmp_mgf_orphan AS
      SELECT spectrum_id, csv_mass, mgf_rt
      FROM tmp_combined t
      WHERE cpd IS NULL;

    CREATE TEMP TABLE tmp_pair AS
      WITH pass_1 AS (
        SELECT t1.spectrum_id AS id_1, t2.spectrum_id AS id_2
        FROM tmp_mgf_orphan t1, tmp_mgf_orphan t2
        WHERE t1.spectrum_id < t2.spectrum_id  AND
              masses_within_error_limit(t1.csv_mass, t2.csv_mass, 5::real) AND
              ABS(t1.mgf_rt - t2.mgf_rt) < 2::real
        )
      SELECT id_1, id_2
      FROM pass_1
      UNION
      SELECT id_2, id_1
      FROM pass_1;

    CREATE TEMP TABLE tmp_all_ids AS
      SELECT spectrum_id AS id FROM tmp_mgf_orphan;

    CREATE TEMP TABLE tmp_cluster AS SELECT * FROM cluster_ids();
    CREATE INDEX tmp_clust_ind ON tmp_cluster(id);
    DROP TABLE tmp_mgf, tmp_mgf_orphan, tmp_pair, tmp_all_ids;

    ---------------------------------------------------------------------
    INSERT INTO molf_spectrum(
        dataset,
        spectrum_id,
        cpd,
        csv_mass,
        component_charge,
        component_recipe,
        component_volume,
        min_z,
        max_z,
        mz,
        mgf_rt,
        csv_rt,
        rt_width_at_half_ht,
        rt_start,
        rt_end,
        num_neutrons,
        volume,
        volume_pct,
        height,
        peaks,
        collision_energy,
        mgf_index,
        mgf_title
      )
      SELECT
        p_dataset,
        spectrum_id,
        COALESCE(cpd,
               (SELECT -cluster_id FROM tmp_cluster WHERE id = spectrum_id))
            AS cpd,
        csv_mass,
        component_charge,
        component_recipe,
        component_volume,
        min_z,
        max_z,
        mz,
        mgf_rt,
        csv_rt,
        rt_width_at_half_ht,
        rt_start,
        rt_end,
        num_neutrons,
        volume,
        volume_pct,
        height,
        peaks,
        collision_energy,
        mgf_index,
        mgf_title
      FROM tmp_combined;

    DROP TABLE tmp_cluster;

    CREATE TEMP SEQUENCE molf_fragment_seq;
    INSERT INTO molf_fragment_peak(
        dataset,
        spectrum_id,
        peak_index,
        mass,
        ion_count
        )
      WITH pass_1 AS (
        SELECT
          dataset,
          spectrum_id,
          (unnest(peaks)).*
        FROM molf_spectrum
        WHERE dataset = p_dataset AND peaks IS NOT NULL
        )
      SELECT 
        dataset, 
        spectrum_id,
        nextval('molf_fragment_seq') AS peak_index,
        mass,
        ion_count
      FROM pass_1
      ORDER BY spectrum_id, mass;

    DROP SEQUENCE molf_fragment_seq;

  -------------------------------------------------------------------------
  --  No tandem data.
  -------------------------------------------------------------------------
  ELSE
    INSERT INTO molf_spectrum(
          dataset,
          spectrum_id,
          cpd,
          csv_mass,
          component_charge,
          component_recipe,
          min_z,
          max_z,
          mz,
          csv_rt,
          rt_width_at_half_ht,
          rt_start,
          rt_end,
          num_neutrons,
          volume,
          volume_pct,
          height
      )
      SELECT
        p_dataset AS dataset,
        nextval('molf_spectrum_seq') AS spectrum_id,
        cpd,
        csv_mass,
        component_charge,
        component_recipe,
        min_z,
        max_z,
        component_mz AS mz,
        csv_rt,
        rt_width_at_half_ht,
        rt_start,
        rt_end,
        component_num_neutrons AS num_neutrons,
        volume,
        volume_pct,
        height
      FROM tmp_molf_rec;
  END IF;

  DROP SEQUENCE molf_spectrum_seq;
  DROP TABLE tmp_molf_rec, tmp_cef;

END;
$$ LANGUAGE plpgsql;


----------------------------------------------------------------------------
--  Convert the data parsed from the "Find Compounds by Molecular Features"
--  data into the standard representation in the tables:
--    * mass_dataset         one row summarizing the whole dataset
--    * observed_mass        peaks in the CSV file
--    * component_spectrum   individual spectra combined into CSV records
--    * fragment_mass        individual fragment peaks
----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION import_molf_data(
    p_dataset                varchar,
    p_bio_context_id         varchar)
    RETURNS void AS $$
DECLARE
  v_num_peaks               integer;
  v_lc_config_id            varchar;
  v_max_volume              real;
BEGIN
  PERFORM * FROM bio_context
  WHERE bio_context_id = p_bio_context_id;
  IF NOT FOUND THEN
    RAISE EXCEPTION 'Non-existent bio_context: %', p_bio_context_id;
  END IF;

  SELECT max(volume)::real INTO v_max_volume
  FROM molf_spectrum
  WHERE dataset = p_dataset;

  --------------------------------------------------------------------------
  --  Populate a generic table (not specific to Agilent data) with only 
  --  one row per feature (CSV rec) from "Find Compounds by Mol. Feature".
  --------------------------------------------------------------------------
  INSERT INTO observed_mass(dataset, peak_id, mass, rt, rt_width_at_half_ht, 
      rt_start, rt_end, quantity, rel_quantity, min_z, max_z, dominant_z)
    SELECT
      dataset,
      cpd                              AS peak_id,
      avg(csv_mass)::double precision  AS mass,
      avg(csv_rt)::real                AS rt,
      avg(rt_width_at_half_ht)::real   AS rt_width_at_half_ht,
      min(rt_start)::real              AS rt_start,
      max(rt_end)::real                AS rt_end,
      sum(volume)::real                AS quantity,
      sum(volume)::real/v_max_volume   AS rel_quantity,
      min(min_z)                       AS min_z,
      max(max_z)                       AS max_z,
      (array_agg(component_charge ORDER BY component_volume DESC))[1]
          AS dominant_z
    FROM molf_spectrum
    WHERE dataset = p_dataset
    GROUP BY dataset, cpd;

  GET DIAGNOSTICS v_num_peaks := ROW_COUNT;

  SELECT lc_config_id INTO v_lc_config_id
  FROM molf_dataset
  WHERE dataset = p_dataset;

  IF NOT FOUND THEN
    RAISE EXCEPTION 'Non-existent dataset: %', p_dataset;
  END IF;

  INSERT INTO mass_dataset(dataset, bio_context_id, lc_config_id, 
                           data_source_type, num_peaks)
  VALUES(p_dataset, p_bio_context_id, v_lc_config_id, 'molf', v_num_peaks);

  INSERT INTO component_spectrum(dataset, spectrum_id, parent_peak_id, charge,
          num_neutrons, collision_energy, num_peaks)
    WITH tmp_pk_cnt AS (
      SELECT
        spectrum_id,
        count(*) AS num_peaks
      FROM molf_fragment_peak
      WHERE dataset = p_dataset
      GROUP BY spectrum_id
      )
    SELECT 
      p_dataset AS dataset,
      spectrum_id,
      cpd  AS parent_peak_id,
      component_charge AS charge,
      num_neutrons,
      collision_energy,
      num_peaks
    FROM molf_spectrum JOIN tmp_pk_cnt USING(spectrum_id)
    WHERE dataset = p_dataset;

  INSERT INTO fragment_mass(dataset, spectrum_id, fragment_peak_id, mass,
                            quantity)
    SELECT
      dataset,
      spectrum_id,
      peak_index AS fragment_peak_id,
      mass,
      ion_count::real AS quantity
    FROM molf_fragment_peak
    WHERE dataset = p_dataset;

END;
$$ LANGUAGE plpgsql;


----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION delete_molf_data(
    p_dataset                varchar)
    RETURNS void
    LANGUAGE SQL AS $$

  DELETE FROM molf_dataset WHERE dataset = $1;
  DELETE FROM molf_spectrum WHERE dataset = $1;
  DELETE FROM molf_fragment_peak WHERE dataset = $1;

$$;

