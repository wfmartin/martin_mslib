
perl ms_config.pl --db_name mst  --opts_file ms_config_opts

./cons.sh

./exclusions_report.pl --db_name mst  --consensus_id mso \
  > ex_init.csv

./process_msms_run.pl --opts_file msms_opts \
   --dataset cmg_00_A  --trace trc_00_A.csv
./exclusions_report.pl --db_name mst  --consensus_id mso \
  > ex_00_A.csv

./process_msms_run.pl --opts_file msms_opts \
   --dataset cmg_01_A  --trace trc_01_A.csv
./exclusions_report.pl --db_name mst  --consensus_id mso \
  > ex_01_A.csv

./process_msms_run.pl --opts_file msms_opts \
   --dataset cmg_02_A  --trace trc_02_A.csv
./exclusions_report.pl --db_name mst  --consensus_id mso \
  > ex_02_A.csv

./process_msms_run.pl --opts_file msms_opts \
   --dataset cmg_03_A  --trace trc_03_A.csv
./exclusions_report.pl --db_name mst  --consensus_id mso \
  > ex_03_A.csv
