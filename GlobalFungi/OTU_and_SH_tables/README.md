
### get SH community matrix

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/OTU_and_SH_tables/TABLE_BY_SAMPLES_AND_SHS.py`

- get complete table (WARNING - could be very large):
`python2.7 TABLE_BY_SAMPLES_AND_SHS.py GF5_RAW_TABLE_TAB_SH.txt.gz GF5_RAW_TABLE_SAMPLES.txt.gz - -`

- get all SH based on samples:
`python2.7 TABLE_BY_SAMPLES_AND_SHS.py GF5_RAW_TABLE_TAB_SH.txt.gz GF5_RAW_TABLE_SAMPLES.txt.gz test_samples.txt -`

- get all samples based on SHs:
`python2.7 TABLE_BY_SAMPLES_AND_SHS.py GF5_RAW_TABLE_TAB_SH.txt.gz GF5_RAW_TABLE_SAMPLES.txt.gz - test_shs.txt`

- get selection of samples and SH:
`python2.7 TABLE_BY_SAMPLES_AND_SHS.py GF5_RAW_TABLE_TAB_SH.txt.gz GF5_RAW_TABLE_SAMPLES.txt.gz test_samples.txt test_shs.txt`

