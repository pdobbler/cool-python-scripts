### REFERENCE BASED

`wget https://www.biomed.cas.cz/mbu/lbwrf/seed/archive/ITS2_PERMANENT_CLUSTERS_SEEDs.7z`

`7z x ITS2_PERMANENT_CLUSTERS_SEEDs.7z`

`vsearch --uchime_ref GF5_ITS2_PERM_CLUSTERS_SEEDS_ALL_SORTED_RENAMED.fa --db ITS2_PERMANENT_CLUSTERS_SEEDs.fa --nonchimeras chimclean.fasta --chimeras chimeras.fasta`

vsearch v2.21.2_linux_x86_64, 5948.9GB RAM, 1152 cores

https://github.com/torognes/vsearch


Reading file REFERENCE.fa 100%

110417743 nt in 633435 seqs, min 40, max 1415, avg 174

Masking 100%

Counting k-mers 100%

Creating k-mer index 100%

Detecting chimeras 100%

Found 347395 (6.5%) chimeras, 4935564 (92.6%) non-chimeras, and 45776 (0.9%) borderline sequences in 5328735 unique sequences.

Taking abundance information into account, this corresponds to 347395 (6.5%) chimeras, 4935564 (92.6%) non-chimeras, and 45776 (0.9%) borderline sequences in 5328735 total sequences.

### get chimera info table

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/ChimeraKiller/get_chimera_info.py`

`python2.7 get_chimera_info.py GF5_ITS2_PERM_CLUSTERS_SEEDS_ALL_SORTED_RENAMED.fa chimeras.fasta GF5_ITS2_PERM_CLUSTERS_SEEDS_ALL_SORTED_RENAMED_CHIMSTATS.txt`

### check divided parts

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/ChimeraKiller/divide_fasta.py`  

`python2.7 divide_fasta.py ITS2_PERMANENT_CLUSTERS_SEEDs_REFERENCE_5PLUS_SIMPLE.fas`

makeblastdb -in ITS2_PERMANENT_CLUSTERS_SEEDs_REFERENCE_5PLUS_SIMPLE.fas.95.0.seeds -dbtype 'nucl' -out REFERENCE_5PLUS_95.fa



