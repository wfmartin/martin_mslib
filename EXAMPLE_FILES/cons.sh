#---------------------------------------------------------------------------
#  Combine MS-only datasets into a consensus, which is loaded into database.
#  The file $consensus_id.consensusXML is created.
#---------------------------------------------------------------------------
perl  ~/martin_mslib/RUN/gen_mso_consensus.pl  \
  --db_name mst  \
  --consensus_id mso  \
  --consensus_parameters_key def_cons_parms  \
  --cons_matching_params_key def_match_params  \
  --bio_context_id human_milk  \
  --inputs_folder /home/bmartin/martin_mslib/test_data  \
  --trace trc_cons.csv  \
  cmg_ms_only_A  cmg_ms_only_B  cmg_ms_only_C

