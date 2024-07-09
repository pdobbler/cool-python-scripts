### GET UNITE 10 All eukaryotes sh dynamic 04042024 complete taxonomy

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/unite10_complete_taxonomy.txt`


### UPDATE ORIGINAL TABLE

BLASTn results processed - GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITE10_PROCESSED.txt.gz

QUERY   HIT     SIMILARITY      COVERAGE        EVALUE  BITSCORE
3d9468560d2615a3ee0c1a00c1a143b7        SH0953900.10FU  98.718  100.0   6.35e-116       416
ba029c6222f7dafd1846ff692af82947        SH0980029.10FU  96.622  100.0   1.84e-64        244
050932ce68b26bfed574a42f9b66cd85        NO_HIT  -       -       -       -
986a02bb45ae7c4fee4d5e78797e3252        SH0982549.10FU  96.266  100.0   8.54e-110       396


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/update_table_processed.py`

`python2.7 update_table_processed.py GF5_RAW_TABLE_PROCESSED.txt.gz GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITE10_PROCESSED.txt.gz GF5_RAW_TABLE_PROCESSED_UNITE10.txt 98.5 90.0 unite10_complete_taxonomy.txt`

It generates "UPDATED_TAX_TABLE.txt" containing reduced taxonomy table...

### GET TAXONOMY TABLES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/get_taxa_table.py`

`python2.7 get_taxa_table.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt UPDATED_TAX_TABLE.txt GF5_RAW_TABLE_TAB_UNITE10_SH.txt 0`

`python2.7 get_taxa_table.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt UPDATED_TAX_TABLE.txt GF5_RAW_TABLE_TAB_UNITE10_GENUS.txt 6`

`python2.7 get_taxa_table.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt UPDATED_TAX_TABLE.txt GF5_RAW_TABLE_TAB_UNITE10_SPECIES.txt 7`
