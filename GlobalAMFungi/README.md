### BINNING

`awk -F'\t' '{ if ($4 >= 97.0 && $5 >= 98.0) print $0 }' SSU_AMF_at_least_90sim_95cov_PROCESSED_VARIANTS.txt > Nuland_tables/SSU_AMF_at_least_97sim_98cov_PROCESSED_VARIANTS_FOR_CLASSIFICATION.txt`
