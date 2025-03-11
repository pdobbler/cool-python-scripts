### CONVERT FASTA

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/EXACT_HIT_MINNING/convert_fasta.py`

`python convert_fasta.py Flexipedes_ITS1.fasta Flexipedes_ITS1_seq_name_file.txt`

### GETTING RES file

`grep -v '>' Flexipedes_ITS1.fasta > Flexipedes_ITS1_seqs.txt`

`dos2unix Flexipedes_ITS1_seqs.txt`

`zgrep --no-group-separator -F -f Flexipedes_ITS1_seqs.txt /mnt/DATA1/RELEASE5/FOR_BLAST_UNITE9_COMPLETE/MAKE_TABLE_DIR/TOTALFINAL/REL5_REANOT_UNITE10/GF5_RAW_TABLE_PROCESSED_UNITE10.txt.gz > Flexipedes_ITS1_RES.txt`

### CONVERT FASTA

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/EXACT_HIT_MINNING/exact_hit_converter.py`

`python2.7 exact_hit_converter.py Flexipedes_ITS1_seq_name_file.txt Flexipedes_ITS1_RES.txt GF5_SAMPLES_METADATA_TOTALFINAL_SAMPLE_ID.txt.gz Flexipedes_ITS1_OUTPUT.txt`

