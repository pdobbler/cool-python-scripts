### GETTING RES file

`grep -v '>' Flexipedes_ITS1.fasta > Flexipedes_ITS1_seqs.txt`

`dos2unix Flexipedes_ITS1_seqs.txt`

`zgrep --no-group-separator -F -f Flexipedes_ITS1_seqs.txt /mnt/DATA1/RELEASE5/FOR_BLAST_UNITE9_COMPLETE/MAKE_TABLE_DIR/TOTALFINAL/REL5_REANOT_UNITE10/GF5_RAW_TABLE_PROCESSED_UNITE10.txt.gz > Flexipedes_ITS1_RES.txt`
