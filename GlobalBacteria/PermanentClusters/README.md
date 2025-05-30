### GREENGENES TAXONOMY BREAKDOWN

Blasting all variants against Greengenes2 20250505 -> process the blast best hit

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/taxonomy.zip`  

`unzip taxonomy.zip`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/taxonomy_breakdown.py`  

`python2.7 taxonomy_breakdown.py GB_VOL1_20250526_CLEANED_uniq_Greengenes2_20250505_PROCESSED.txt taxonomy.tsv 1 GB_VOL1_20250526_CLEANED_ALL_Phylum_breakdown.txt 188.0`

### FILTER NON-SINGLETONS

`awk '/^>/ {keep = ($0 !~ /\|P_1\|V_1\|/)} keep' seqs.fa > filtered_seqs.fa`


### SPLIT FASTA BY TAXONOMY

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/split_by_taxonomy.py`  

`python2.7 split_by_taxonomy.py GB_VOL1_20250526_CLEANED.fa.gz.qualified GB_VOL1_20250526_CLEANED_uniq_Greengenes2_20250505_PROCESSED.txt taxonomy.tsv 1 GB_VOL1_PHYLUM 188.0`

### STATS

```
for file in *.fas; do
    count=$(grep -c '^>' "$file")
    echo -e "$file\t$count" >> counts.txt
done
```

### GROUP SEQUENCES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/group_qualified_sequences_fastest_optimized.py`

`parallel -j $(nproc) python group_qualified_sequences_fastest_optimized.py {} 97.0 false ::: *.fas > output.txt`


