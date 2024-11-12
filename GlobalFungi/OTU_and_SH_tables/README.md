
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


### GET OTUTABLE FROM CLUSTERED FASTA

header e.g.: >GF05023557S|Fernan_2020_JZ12|e7c017be559f77ec|OTU0027655

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/OTU_and_SH_tables/GET_OTUTAB_FROM_SEQUENCES.py`

`python2.7 GET_OTUTAB_FROM_SEQUENCES.py GF5_ALL_SAMPLES_ITS2_CLUSTERED_ECTO.fa GF5_ALL_SAMPLES_ITS2_CLUSTERED_ECTO_OTUTABLE.txt`


### GET REDUCED OTUTABLE FROM CLUSTERED FASTA

header e.g.: >GF05028373S|Chen12_2022_AKD1|c4327a8b28d158a0|OTU0096999

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/OTU_and_SH_tables/GET_OTUTAB_FROM_SEQUENCES_REDUCED_FORMATE.py`

all clusters and singletons included (no_singletons:false;-:means no selected OTUs list provided)
------------------------------------
`python2.7 GET_OTUTAB_FROM_SEQUENCES_REDUCED_FORMATE.py test_cl_fasta.fa false -`

