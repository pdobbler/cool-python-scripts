### GREENGENES TAXONOMY BREAKDOWN

Blasting all variants against Greengenes2 20250505 -> process the blast best hit

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/taxonomy.zip`  

`unzip taxonomy.zip`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/taxonomy_breakdown.py`  

`python2.7 taxonomy_breakdown.py GB_VOL1_20250526_CLEANED_uniq_Greengenes2_20250505_PROCESSED.txt taxonomy.tsv 1 GB_VOL1_20250526_CLEANED_ALL_Phylum_breakdown.txt 188.0`
