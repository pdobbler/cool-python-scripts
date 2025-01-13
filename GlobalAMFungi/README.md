### BINNING

/mnt/DATA/projects/avetrot/AMF_DATABASE/FINAL_DATASETS/AMF_RELEASE1_CORRECTED  

`awk -F'\t' '{ if ($4 >= 97.0 && $5 >= 98.0) print $0 }' SSU_AMF_at_least_90sim_95cov_PROCESSED_VARIANTS.txt > Nuland_tables/SSU_AMF_at_least_97sim_98cov_PROCESSED_VARIANTS_FOR_CLASSIFICATION.txt`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalAMFungi/get_VTX_table_from_processed_blast.py`

`python2.7 get_VTX_table_from_processed_blast.py SSU_qm20_renamed_PRIMARY_FORBIN_NO_DUPL_AND_AMB.fa Nuland_tables/SSU_AMF_at_least_97sim_98cov_PROCESSED_VARIANTS_FOR_CLASSIFICATION.txt Nuland_tables/SSU_VTXTABLE_CLASSIFICATION_97sim.txt`

