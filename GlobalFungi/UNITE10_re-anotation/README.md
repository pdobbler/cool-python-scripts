### GET UNITE 10 All eukaryotes sh dynamic 04042024 complete taxonomy

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/unite10_complete_taxonomy.txt`


### UPDATE ORIGINAL TABLE

BLASTn results processed - GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITE10_PROCESSED.txt.gz

| QUERY                          | HIT             | SIMILARITY | COVERAGE | EVALUE     | BITSCORE |
|--------------------------------|-----------------|------------|----------|------------|----------|
| 3d9468560d2615a3ee0c1a00c1a143b7 | SH0953900.10FU  | 98.718     | 100.0    | 6.35e-116  | 416      |
| ba029c6222f7dafd1846ff692af82947 | SH0980029.10FU  | 96.622     | 100.0    | 1.84e-64   | 244      |
| 050932ce68b26bfed574a42f9b66cd85 | NO_HIT          | -          | -        | -          | -        |
| 986a02bb45ae7c4fee4d5e78797e3252 | SH0982549.10FU  | 96.266     | 100.0    | 8.54e-110  | 396      |



`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/update_table_processed.py`

`python2.7 update_table_processed.py GF5_RAW_TABLE_PROCESSED.txt.gz GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITE10_PROCESSED.txt.gz GF5_RAW_TABLE_PROCESSED_UNITE10.txt 98.5 90.0 unite10_complete_taxonomy.txt`


**Output table looks like this...(HEADER IS NOT PART OF THE TABLE)**


| QUERY                          | SAMPLES              | ABUNDANCES      | MARKER | SH       | SEQUENCE |
|--------------------------------|----------------------|-----------------|--------|----------|----------|
| 95e494b4061325c3d39e5f304528191d | 29                   | 1               | ITS1   | 0        | AAAAA... |
| 0fa49bed213860e0db1c33527cfd098f | 30;31;32;33;34;35;36;37;38;39 | 1;1;1;1;1;1;3;1;1;1 | ITS1   | 5        | CCGAG... |
| d2a2561d72610a4dea34c77bd6f3cfc2 | 40                   | 1               | ITS2   | 0        | ACACC... |
| 560dbd59fca8359abdf30da51851ae83 | 41                   | 2               | ITS1   | 0        | CCGAA... |
| 9c783f603a09389a9a2c09aacb3fe250 | 42;43                | 1;1             | ITS2   | 0        | CCACC... |
| df5271d8c7128f9a7e8bbe6e644107da | 44                   | 1               | ITS2   | 6        | AGCCT... |




It generates "UPDATED_TAX_TABLE.txt" containing reduced taxonomy table...


### GET TAXONOMY TABLES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/get_taxa_table.py`

`python2.7 get_taxa_table.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt UPDATED_TAX_TABLE.txt GF5_RAW_TABLE_TAB_UNITE10_SH.txt 0`

`python2.7 get_taxa_table.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt UPDATED_TAX_TABLE.txt GF5_RAW_TABLE_TAB_UNITE10_GENUS.txt 6`

`python2.7 get_taxa_table.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt UPDATED_TAX_TABLE.txt GF5_RAW_TABLE_TAB_UNITE10_SPECIES.txt 7`


### GET TAXONOMY TABLES

**`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/get_taxa_table_for_marker.py`**
