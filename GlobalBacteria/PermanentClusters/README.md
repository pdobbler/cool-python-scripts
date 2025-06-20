### GREENGENES TAXONOMY BREAKDOWN

Blasting all variants against Greengenes2 20250505 -> process the blast best hit

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/taxonomy.zip`  

`unzip taxonomy.zip`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/taxonomy_breakdown.py`  

`python2.7 taxonomy_breakdown.py GB_VOL1_20250526_CLEANED_uniq_Greengenes2_20250505_PROCESSED.txt taxonomy.tsv 1 GB_VOL1_20250526_CLEANED_ALL_Phylum_breakdown.txt 188.0`

### RANK AND SORT

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/generate_SV_ranks_relative_abund.py`

`python2.7 generate_SV_ranks_relative_abund.py GB_VOL1_20250526_CLEANED.fa.gz false`

373 267 631 variants (GB_VOL1_20250526_CLEANED.fa.gz.all)

### FILTER NON-SINGLETONS

`awk '/^>/ {keep = ($0 !~ /\|V_1\|S_1\|/)} keep' GB_VOL1_20250526_CLEANED.fa.gz.all > GB_VOL1_20250526_CLEANED_ranked_multi.fa`

51 603 682 variants (GB_VOL1_20250526_CLEANED_ranked_multi.fa)

`awk '/^>/ {keep = ($0 ~ /\|V_1\|S_1\|/)} keep' GB_VOL1_20250526_CLEANED.fa.gz.all > GB_VOL1_20250526_CLEANED_ranked_single.fa`


### SPLIT FASTA BY TAXONOMY

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/split_by_taxonomy.py`  

`python2.7 split_by_taxonomy.py GB_VOL1_20250526_CLEANED_ranked_multi.fa GB_VOL1_20250526_CLEANED_uniq_Greengenes2_20250505_PROCESSED.txt taxonomy.tsv 1 GB_VOL1_PHYLUM 188.0`

### STATS

```
for file in *.fas; do
    count=$(grep -c '^>' "$file")
    echo -e "$file\t$count" >> counts.txt
done
```

### GROUP SEQUENCES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/group_qualified_sequences_fastest_optimized.py`

```
export TMPDIR=/mnt/DATA1/tmp
parallel --tmpdir /mnt/DATA1/tmp -j $(nproc) python group_qualified_sequences_fastest_optimized.py {} 97.0 false ::: *.fas > output.txt
```

### TRY TO CONTINUE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/split_by_last_seed.py`

`python2.7 split_by_last_seed.py GB_VOL1_PHYLUM_p__Pseudomonadota.fas GB_VOL1_PHYLUM_p__Pseudomonadota.fas.97.0.seed_seqs`

```
mkdir Pseudomonadota_CONTINUE
mv GB_VOL1_PHYLUM_p__Pseudomonadota.fas.done Pseudomonadota_CONTINUE
mv GB_VOL1_PHYLUM_p__Pseudomonadota.fas.undone Pseudomonadota_CONTINUE
mv GB_VOL1_PHYLUM_p__Pseudomonadota.fas.input_seeds_seqs Pseudomonadota_CONTINUE
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/group_qualified_sequences_fastest_optimized_seeds.py`

`python group_qualified_sequences_fastest_optimized_seeds.py GB_VOL1_PHYLUM_p__Pseudomonadota.fas.done 97.0 false`

`python group_qualified_sequences_fastest_optimized_seeds.py GB_VOL1_PHYLUM_p__Pseudomonadota.fas.undone 97.0 false GB_VOL1_PHYLUM_p__Pseudomonadota.fas.input_seeds_seqs`


### MULTIPLE - CONTINUE

```
for file in *.fas
do  
 python2.7 split_by_last_seed.py ${file} ${file}.97.0.seed_seqs
done
```

```
for file in *.done
do
  echo "python group_qualified_sequences_fastest_optimized_seeds.py ${file} 97.0 false"
  echo "python group_qualified_sequences_fastest_optimized_seeds.py ${file%%.done}.undone 97.0 false ${file%%.done}.input_seeds_seqs"
done > continue.sh

mkdir -p /mnt/DATA1/tmp
export TMPDIR=/mnt/DATA1/tmp
cat continue.sh | parallel --tmpdir /mnt/DATA1/tmp > output.txt
```



