Scripts for creating PERMANENT CLUSTERS for GlobalFungi


### PROCESS SEQUENCES TO GET SV RANK


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/generate_SV_ranks_relative_abund.py`


`python2.7 generate_SV_ranks_relative_abund.py GF5_ALL_SAMPLES.fa.its2.gz true`

input formate (FASTA):
```
>GF05024953S|Garcia_2022_AIZ1|814506c13008882b
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTACATATAA...
>GF05020035S|Hannul_2020_ZZ12|299ff196f43d39d7
AAAAAAAAAAAAAAAAAAAAAAAAAAAAGTCAAGAAGGCCGATA...
```

outputs:
- GF5_ALL_SAMPLES.fa.its2.gz.nonqualified
- GF5_ALL_SAMPLES.fa.its2.gz.qualified


### GROUPING QUALIFIED SEQUENCES (IMPORTANT - biopython must be Version: 1.78 or 1.79!!!)
`pip show biopython`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/group_qualified_sequences_fastest_improved.py`

```
python group_qualified_sequences_fast.py REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz.qualified 97.0 false
python group_qualified_sequences_fast.py REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz.qualified 98.5 false

```

outputs:
- GF5_ALL_SAMPLES.fa.its2.gz.qualified.clustered2
- GF5_ALL_SAMPLES.fa.its2.gz.qualified.seeds2


performance:
- group_qualified_sequences_fast.py    - # clusters 633435 total time 11931499.46307 sec (3314.3h)
- group_qualified_sequences_faster.py  - # clusters 633435 total time 4223531.893308 sec (1173.2h)
- group_qualified_sequences_fastest.py - # clusters 633435 total time 2380216.499521 sec (661.2h)

### SEEDS IDENTIFICATION

```
blastn -query GF5_ALL_SAMPLES.fa.its2.gz.qualified.seeds2 -db /mnt/DATA/DATABASES/UNITE9.0/UNITE_9.0_All_euk_ecology -out GF5_ALL_SAMPLES.fa.its2.gz.qualified.seeds2_unite9.0.txt -evalue 1E-5 -outfmt 6 -num_threads 128 -max_target_seqs 10

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr GF5_ALL_SAMPLES.fa.its2.gz.qualified.seeds2_unite9.0.txt | sort -u -k1,1 --merge > GF5_ALL_SAMPLES.fa.its2.gz.qualified.seeds2_unite9.0_best.txt
```


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/PROCESS_BLAST_RESULT.py`


```
python2.7 PROCESS_BLAST_RESULT.py GF5_ALL_SAMPLES.fa.its2.gz.qualified.seeds2 GF5_ALL_SAMPLES.fa.its2.gz.qualified.seeds2_unite9.0_best.txt GF5_ALL_SAMPLES.fa.its2.gz.qualified.seeds2_unite9.0_PROCESSED.txt ITS2
```

### BINNING TO SEEDS


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/split_fasta_by_group_size.py`


```
python2.7 split_fasta_by_group_size.py REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz.nonqualified 1297870
python2.7 split_fasta_by_group_size.py GF5_ALL_SAMPLES.fa.its2.gz.nonqualified 1565100
```

```
# 324467145
# 324467145
```

```
makeblastdb -in REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz.qualified.seeds -dbtype 'nucl' -out BINNING/REL4_ITS2_QUALIFIED_SEEDS
```

```
# for file in *.fas
# do  
#  echo "blastn -query ${file} -db REL5_ITS2_QUALIFIED_SEEDS -out ${file%%.fas}.seeds.txt -evalue 1E-5 -outfmt 6 -num_threads 2 -max_target_seqs 10"
# done > blast_command.sh

# cat blast_command.sh | parallel
```

```
# export LC_ALL=en_US.UTF-8
# export LANG=en_US.UTF-8

# for file in *.txt
# do  
#  echo "sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr ${file} | sort -u -k1,1 --merge > ${file%%.txt}_best.tab"
# done > sort_command.sh

# cat sort_command.sh | parallel
```

*BLAST AND BEST HITS*


```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

for file in *.fas
do
 echo "blastn -query ${file} -db REL5_ITS2_QUALIFIED_SEEDS -outfmt 6 -evalue 1E-5 -num_threads 2 -max_target_seqs 10 | sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge > ${file%%.fas}_best.tab"
done > blast_and_sort_command.sh

mkdir -p /mnt/DATA1/tmp
export TMPDIR=/mnt/DATA1/tmp
cat blast_and_sort_command.sh | parallel --tmpdir /mnt/DATA1/tmp
```

```
# make blast shorter
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/update_blast_results.py`

```
python2.7 update_blast_results.py REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated_nonqualified_seeds_best.txt
```

```
# binning
```


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/PROCESS_BLAST_RESULT.py`


```
python2.7 PROCESS_BLAST_RESULT.py REL4_ITS2_nonqualified.fa REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated_nonqualified_seeds_best.txt.shorter.txt REL4_ITS2_nonqualified_PROCESSED.txt ITS2
```

324467145 sequences loaded correctly - 0 sequnces are empty - Omitted!
FOUND BLAST RESULTS FOR 315998292 SEQUENCES OUT OF 324467145

```
# 97 %
# https://drive.google.com/file/d/1DJLiuxMmQkig8zWC4jg7-A8yvFtSk5xN/view?usp=sharing
# gdrive_download 1DJLiuxMmQkig8zWC4jg7-A8yvFtSk5xN bin_fasta_to_OTU_by_processed_blast.py
https://drive.google.com/file/d/1V7mLOgpU4A0sYsMjw-9oNfy31OLjKBDq/view?usp=sharing
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/bin_fasta_to_CLUSTERS_by_processed_blast_const.py`

```
python2.7 bin_fasta_to_CLUSTERS_by_processed_blast_const.py GF5_ALL_SAMPLES.fa.its2.gz.nonqualified.gz GF5_ALL_SAMPLES_its2_nonqualified_seedsblast_PROCESSED.txt 198.5
```

### CLUSTER STATS


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/count_clusters.py`


```
for file in *.fa
do  
 echo "python2.7 count_clusters.py $file"
done > count_command.sh

cat count_command.sh | parallel
```

Sequences loaded 4411611 clusters found 152306
DONE
Sequences loaded 123981967 clusters found 1
DONE
Sequences loaded 200485178 clusters found 145269
DONE


### GET SHs FOR CLUSTERS

```
blastn -query REL4_ITS2_FINAL_qualified_clustered_seed.fa -db /mnt/DATA/DATABASES/UNITE8.2/UNITE_8.2_All_eukaryotes -out REL4_ITS2_FINAL_qualified_clustered_seed_unite8.2_all_euk_forSH.txt -evalue 1E-5 -outfmt 6 -num_threads 128
```

3.7G Nov 14 08:47 REL4_ITS2_FINAL_qualified_clustered_seed_unite8.2_all_euk_forSH.txt

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/get_SHs_from_complete_blast.py`

```
python2.7 get_SHs_from_complete_blast.py REL4_ITS2_FINAL_qualified_clustered_seed.fa REL4_ITS2_FINAL_qualified_clustered_seed_unite8.2_all_euk_forSH.txt 198.5
```

### GET CLUSTERS FOR SHs

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/get_clusters_for_SHs.py`


```
python2.7 get_clusters_for_SHs.py REL4_ITS2_FINAL_qualified_clustered_seed.fa.shs
```


### GET SAMPLES AND STUDIES FOR CLUSTERS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/get_samples_and_studies_for_clustered_vars.py`

```
python2.7 get_samples_and_studies_for_clustered_vars.py FINAL/REL4_ITS2_FINAL_qualified_clustered_and_binned.fa REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz
```

Sequence vars loaded... 204896789
Sequences processed... total: 2884581662 clustered: 2734213852 non-clustered: 150367810
DONE

98.5
Sequence vars loaded... 250272474
Sequences processed... total: 2884581662 clustered: 2787623835 non-clustered: 96957827
DONE. 250272474


**statistic in selected samples**

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/get_samples_and_studies_for_clustered_vars_with_list.py`

`python2.7 get_samples_and_studies_for_clustered_vars_with_list.py FINAL/GF5_ALL_SAMPLES_its2_nonqualified_clustered_and_binned.fa.gz GF5_ALL_SAMPLES.fa.its2.gz Eligible_samples_ITS2_clusters.txt`

### SUSPICIOUS CLUSTERS - PARTS OF ITS

```
zgrep --no-group-separator -B 1 'TGGGGCTTTGTCACCCGCTCTGTAGGCCCGGCCGGCGCTTGCCGATCAACCAAATTTTTATCCA' REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz > SUSPICIOUS_CLUSTERS/CL018450_variants_complete_all.fa

grep '>' CL018450_variants_complete_all.fa | awk -F'|' '{print $2}' | sort |  uniq -c | wc -l
```

```
titles
GF4S06602b|Sun_2021_PK|ERR4887469.331809
```

```
sed 's/\r//g' CL018450_titles.txt | awk '{print $1" "}' > CL018450_titles_space.txt

dos2unix CL018450_titles_space.txt

zgrep --no-group-separator -A 1 -F -f CL018450_titles_space.txt Sun_2021_PK_qm20_renamed_check_all_extracted.fa.gz > CL018450_Sun_2021_PK_raw.fa
```

### UNITE 9


```
makeblastdb -in UNITEv9_sh_dynamic_all_tax_final.fas -dbtype 'nucl' -out UNITE9_dyn_tax

blastn -query REL4_ITS2_FINAL_qualified_clustered_seed.fa -db UNITE9_dyn_tax -out REL4_ITS2_FINAL_qualified_clustered_seed_unite9_dyn_tax.txt -evalue 1E-5 -outfmt 6 -num_threads 256

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr REL4_ITS2_FINAL_qualified_clustered_seed_unite9_dyn_tax.txt | sort -u -k1,1 --merge > REL4_ITS2_FINAL_qualified_clustered_seed_unite9_dyn_tax_best.txt
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/PROCESS_BLAST_RESULT.py`

```
python2.7 PROCESS_BLAST_RESULT.py REL4_ITS2_FINAL_qualified_clustered_seed.fa REL4_ITS2_FINAL_qualified_clustered_seed_unite9_dyn_tax_best.txt REL4_ITS2_seed_unite9_dyn_tax_PROCESSED.txt ITS2
```

### COUNTS
```
grep '>' GF5_ALL_SAMPLES_ITS2_CLUSTANDBINNED.fa | grep 'BINNED_' | awk -F'|' '{
    split($3, v, "_");
    split($4, s, "_");
    id[$1]++;
    vSum[$1]+=v[2];
    sSum[$1]+=s[2];
}
END {
    for (i in id) {
        printf "%s\t%d\t%d\n", i, vSum[i], sSum[i];
    }
}' > BINNED_output.tsv
```

```
grep '>' GF5_ALL_SAMPLES_ITS2_CLUSTANDBINNED.fa | grep -v 'BINNED_' | awk -F'|' '{
    split($3, v, "_");
    split($4, s, "_");
    id[$1]++;
    vSum[$1]+=v[2];
    sSum[$1]+=s[2];
}
END {
    for (i in id) {
        printf "%s\t%d\t%d\n", i, vSum[i], sSum[i];
    }
}' > CLUSTERED_output.tsv
```


### SAMPLES - BINNED+CL vs NOT

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/get_sample_binned_proportion.py`


```
python2.7 get_sample_binned_proportion.py /mnt/DATA/projects/avetrot/REL4_RAW_and_CLUSTERING_FINAL/RELEASE4_RAW/REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz REL4_ITS2_FINAL_qualified_clustered_and_binned.fa SAMPLES_BINNEDANDCL_VS_NOT.txt
```

###  CLUSTERING 98.5 %

```
for file in *.fas
do  
 echo "blastn -query ${file} -db REL4_ITS2_QUALIFIED_SEEDS_98.5 -out ${file%%.fas}.seeds_98.5.txt -evalue 1E-5 -outfmt 6 -num_threads 4 -max_target_seqs 10"
done > blast_command.sh

cat blast_command.sh | parallel
```

### RESULTS

```
python2.7 count_clusters.py REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz.qualified.clustered2
```

Sequences loaded 4411611 clusters found 228289

```
python2.7 count_clusters.py REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz.nonqualified.binned
```

Sequences loaded 245860863 clusters found 211559



### NON-BINNED PROCESSING

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/split_by_sample_size.py`

`python2.7 split_by_sample_size.py GF5_ALL_SAMPLES.fa.its2.gz.nonqualified.gz.notbinned.gz 4`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/sort_by_SV_ranks.py`

`python2.7 sort_by_SV_ranks.py GF5_ALL_SAMPLES.fa.its2.gz.nonqualified.gz.notbinned.gz.qualified`


### GET TABLE FROM CLUSTERED VARIANTS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/get_table_from_clustered_variants.py`

`python2.7 get_table_from_clustered_variants.py SSU_PERMCLUSTERED_97_AMF_SELECTION.fa SSU_qm20_renamed_PRIMARY_FORBIN_NO_DUPL_AND_AMB.fa.multi.gz`

- For each sample in Table, calculate: the total number of clusters, the total number of sequences, the total number of singletons (within the sample), the total number of doubletons (within the sample):

```
awk 'NR == 1 { for (i = 2; i <= NF; i++) {header[i] = $i;}} 
NR > 1 {
  for (i = 2; i <= NF; i++) {
    if ($i != 0) nonzero[i]++;
    sum[i] += $i;
    if ($i == 1) ones[i]++;
    if ($i == 2) twos[i]++;
  }
}
END {
  for (i = 2; i <= NF; i++) {
    print header[i], "non-zero:", nonzero[i] + 0, "sum:", sum[i] + 0, "ones:", ones[i] + 0, "twos:", twos[i] + 0;
  }
}' SSU_PERMCLUSTERED_97_AMF_SELECTION.fa.TABLE.txt > SSU_PERMCLUSTERED_97_AMF_SELECTION_TABLE_stats.txt
```

- Remove all samples from the table that have fewer than 1000 sequences, and randomly subsample all samples so that they have exactly 1000 sequences:


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/filter_columns.awk`

`chmod +x filter_columns.awk`

`awk -f filter_columns.awk 500 SSU_PERMCLUSTERED_97_AMF_SELECTION.fa.TABLE.txt > SSU_PERMCLUSTERED_97_AMF_SELECTION_TABLE_min500seqs.txt`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/reconstruct_virtual_dataset.py`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/get_subsampled_FASTA.py`

`python2.7 reconstruct_virtual_dataset.py SSU_PERMCLUSTERED_97_AMF_SELECTION_TABLE_min500seqs.txt SSU_PERMCLUSTERED_97_AMF_SELECTION_TABLE_min500seqs.fa`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/get_subsampled_table_from_FASTA.py`

`python2.7 get_subsampled_table_from_FASTA.py SSU_PERMCLUSTERED_97_AMF_SELECTION_TABLE_min500seqs.fa 500 SSU_PERMCLUSTERED_97_AMF_SELECTION_TABLE_SUBSAMPLED_500.txt`


### VARIANTS IN TWO SAMPLES

`awk -F '\t' '$3 ~ /k__Fungi/ && $4 >= 85.0 {print $3}' GF5_ALL_SAMPLES_its2_S2_EukaryomITS_best_PROCESSED.txt | awk -F';' '{print $4}' | sort | uniq -c > fungal_classes_above85sim.txt`

`grep -v "k__Fungi" GF5_ALL_SAMPLES_its2_S2_EukaryomITS_best_PROCESSED.txt | awk -F '\t' '{print $2}' > GF5_ALL_SAMPLES_its2_S2_NONFUNGAL_varaint_names.txt`

`grep "k__Fungi" *_PROCESSED.txt | awk -F '\t' '$3 ~ /k__Fungi/ && $4 < 85.0 {print $2}' > GF5_ALL_SAMPLES_its1_notbinned_S1_VMULTI_below85_varaint_names.txt`  

`cat *_NONFUNGAL_varaint_names.txt *_below85_varaint_names.txt > unidentified_names.txt`

get sequences  

`grep --no-group-separator -A 1 -F -f unidentified_names.txt GF5_ALL_SAMPLES_its1_nonqualified_notbinned_S1_VMULTI.fa > G01_UNIDENTIFIED/S1_unidentified.fa`

`awk -F '\t' '$3 ~ /c__Agaricomycetes/ && $4 >= 85.0 {print $2}' GF5_ALL_SAMPLES_its1_notbinned_S1_VMULTI_best_PROCESSED.txt > c__Agaricomycetes_varaint_names.txt`  
`grep --no-group-separator -A 1 -F -f c__Agaricomycetes_varaint_names.txt GF5_ALL_SAMPLES_its1_nonqualified_notbinned_S1_VMULTI.fa > G02_Agaricomycetes/S1_Agaricomycetes.fa`  

simplified  

```
grep --no-group-separator -A 1 -F -f \
<(awk -F '\t' '$3 ~ /c__Sordariomycetes/ && $4 >= 85.0 {print $2}' GF5_ALL_SAMPLES_its1_notbinned_S1_VMULTI_best_PROCESSED.txt) \
GF5_ALL_SAMPLES_its1_nonqualified_notbinned_S1_VMULTI.fa > S1_Sordariomycetes.fa
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/add_text_to_fasta_title.py`

`python2.7 add_text_to_fasta_title.py S2_Sordariomycetes_varaints.fa.sorted.98.5.clustered S2_Sord`

`python2.7 add_text_to_fasta_title.py S2_Sordariomycetes_varaints.fa.sorted.98.5.seeds S2_Sord`

### GET TABLE FROM CLUSTERED VARIANTS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/get_otutab_for_clusters_REDUCED_FORMATE.py`

`python2.7 get_otutab_for_clusters_REDUCED_FORMATE.py PERMANENT_CLUSTERS/GF5_ALL_SAMPLES.fa.its2.gz FINAL/GF5_ALL_SAMPLES_ITS2_CLUSTANDBINNED.fa.gz GF5_ALL_SAMPLES_ITS2_PERM_OTUTAB_REDUCED_FORMATE.txt`

### S1 PROCESSING

`grep '>' GF5_ALL_SAMPLES_its2_notbinned_S4_S3_without_S2.fa.gz.notbinned | grep -v '|V_1|' > VMULTI.txt`

`grep --no-group-separator -A 1 -F -f VMULTI.txt GF5_ALL_SAMPLES_its2_notbinned_S4_S3_without_S2.fa.gz.notbinned > GF5_ALL_SAMPLES_its2_notbinned_S4_S3_S2_VMULTI.fa`

`blastn -query GF5_ALL_SAMPLES_its2_notbinned_S4_S3_S2_VMULTI.fa -db /mnt/DATA/DATABASES/Eukaryome/EUK_ITS_v1.9_General -outfmt 6 -evalue 1E-5 -num_threads 2 -max_target_seqs 10 | sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge > GF5_ALL_SAMPLES_its2_notbinned_S4_S3_S2_VMULTI_best.tab`

`awk -F'\t' '{print $2}' GF5_ALL_SAMPLES_its1_notbinned_S1_VMULTI_best.tab | awk -F';' '{print $4}' | sort | uniq -c > classes.txt`

### SEEDS PROCESSING

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/rename_permanent_clusters_seeds.py`

`python2.7 rename_permanent_clusters_seeds.py GF5_ITS2_PERM_CLUSTERS_SEEDS_ALL_SORTED.fa 8`

### SEEDS CHIMERA CHECK

`vsearch --uchime_ref GF5_ITS2_PERM_CLUSTERS_SEEDS_ALL_SORTED_RENAMED.fa --db ITS2_PERMANENT_CLUSTERS_SEEDs_S5.fa --nonchimeras chimclean.fasta --chimeras chimeras.fasta`

vsearch v2.21.2_linux_x86_64, 5948.9GB RAM, 1152 cores  
https://github.com/torognes/vsearch  

Reading file ITS2_PERMANENT_CLUSTERS_SEEDs.fa 100%  
110417743 nt in 633435 seqs, min 40, max 1415, avg 174  
Masking 100%  
Counting k-mers 100%  
Creating k-mer index 100%  
Detecting chimeras 100%  
Found 754865 (12.7%) chimeras, 5151188 (86.5%) non-chimeras,  
and 46143 (0.8%) borderline sequences in 5952196 unique sequences.  
Taking abundance information into account, this corresponds to  
754865 (12.7%) chimeras, 5151188 (86.5%) non-chimeras,  
and 46143 (0.8%) borderline sequences in 5952196 total sequences.  







