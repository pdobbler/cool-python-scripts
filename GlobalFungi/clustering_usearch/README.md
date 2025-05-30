### MAKE VARIANTS WITH SIZE INFO

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/make_uniques.py`

`python2.7 make_uniques.py GF5_ALL_SAMPLES.fa.its1.gz` 

### SEPARATE SINGLE AND MULTIVARIANTS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/get_multi_variants.py`

`python2.7 get_multi_variants.py GF5_ALL_SAMPLES.fa.its1.gz.uniq.gz` 

### CLUSTERING

`usearch -cluster_otus GF5_ALL_SAMPLES.fa.its1.gz.uniq.gz.multi -minsize 1 -otus GF5_ALL_SAMPLES.fa.its1_uniq_multi_otus.fa -relabel Otu -uparseout GF5_ALL_SAMPLES_its1_minsize2_uparse.txt`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/PROCESS_UPARSE_v11.0.667_RESULTS.py`

`python2.7 PROCESS_UPARSE_v11.0.667_RESULTS.py GF5_ALL_SAMPLES.fa.its1.gz GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED.fa.gz GF5_ALL_SAMPLES.fa.its1.gz.uniq.gz.multi GF5_ALL_SAMPLES_its1_minsize2_uparse.txt`

Sequences were renamed based on the generated OTUs - total: 1264305616  
Clustered to OTUs     : 1067936446  
Removed chimeras      : 16798002  
Unclustered singletons: 179571168  

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/GET_OTUS_MOST_ABUND_SEQUNCES.py`

`python2.7 GET_OTUS_MOST_ABUND_SEQUNCES.py GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED.fa.gz GF5_ALL_SAMPLES_ITS1_CLUSTERED_MOST_ABUND.fa`

```
mkdir SPLIT
makeblastdb -in GF5_ALL_SAMPLES_ITS1_CLUSTERED_MOST_ABUND.fa -dbtype 'nucl' -out SPLIT/ITS1_CLUSTERS
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/split_fasta_by_group_size.py`

`mv GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED.fa.gz.singletons GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED_singletons.gz`

`python2.7 split_fasta_by_group_size.py GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED_singletons.gz 900000`

```
for file in *.fas
do  
 echo "blastn -query ${file} -db ITS1_CLUSTERS -out ${file%%.fas}.GF5_ALL_SAMPLES_ITS1_CLUSTERS.txt -evalue 1E-5 -outfmt 6 -num_threads 2 -max_target_seqs 10"
done > blast_command.sh

cat blast_command.sh | parallel
```

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

for file in *.txt
do  
 echo "sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr ${file} | sort -u -k1,1 --merge > ${file%%.txt}_best.tab"
done > sorting.sh

cat sorting.sh | parallel
```


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/PROCESS_BLAST_RESULT.py`

`python2.7 PROCESS_BLAST_RESULT.py GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED_singletons.gz GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED_singletons_ITS1_CLUSTERS_best.tab GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED_singletons_PROCESSED.txt ITS1`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/bin_fasta_to_OTU_by_processed_blast.py`

`python2.7 bin_fasta_to_OTU_by_processed_blast.py GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED_singletons.gz GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED_singletons_PROCESSED.txt`

### GET CLASSIC OTUTABLE FROM CLUSTERED FASTA
  
header e.g.: >GB01020442S|An_2019_1acp_Bact|SRR5920425.6512|POS=5|POS=253|OTU00381  
  
`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/GET_OTUTAB_FROM_SEQUENCES.py`
  
`python2.7 GET_OTUTAB_FROM_SEQUENCES.py GB_VOL1_20251903_CLEAN_min12k_max15k_CLUSTERED_AND_BINNED_BACONLY.fa.gz GB_VOL1_20251903_CLEAN_min12k_max15k_CLUSTERED_AND_BINNED_BACONLY_OTUTAB.txt`  


### GET REDUCED OTUTABLE FROM CLUSTERED FASTA

header e.g.: >GF05028373S|Chen12_2022_AKD1|c4327a8b28d158a0|OTU0096999

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/OTU_and_SH_tables/GET_OTUTAB_FROM_SEQUENCES_REDUCED_FORMATE.py`

all clusters and singletons included (no_singletons:false;-:means no selected OTUs list provided)

`python2.7 GET_OTUTAB_FROM_SEQUENCES_REDUCED_FORMATE.py GF5_ALL_SAMPLES_ITS2_CLUSTERED.fa.gz false -`

### CLUSTER NONBINNED SAMPLES INDIVIDUALY

`awk '/^>/ {split($1,a,"|"); sample=a[1]; sub(/^>/,"",sample); print > sample".fasta"; next} {print >> sample".fasta"}' Labouyrie_2023_BBS_minsize2_CLUSTERED_singletons.fa.gz.notbinned`

```
for file in *.fasta; do
    python2.7 make_uniques.py "$file"
done
```

```
for file in *.fasta.uniq; do
    usearch -cluster_otus ${file} -minsize 1 -otus ${file%%.fasta.uniq}_otus.fa -relabel Otu -uparseout ${file%%.fasta.uniq}_uparse.txt
done
```

```
for file in *_otus.fa; do 
    count=$(grep -c "^>" "$file")
    echo -e "${file%%_otus.fa}\t${count}"
done > OTU_counts.txt
```

### SUBSAPLE FASTA BASED ON SAMPLE SIZE (discard samples under the treshold)

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/get_subsampled_FASTA.py`

`python get_subsampled_FASTA.py 12000 Labouyrie_2023_BBS_qm20_renamed_correct_ok_223_283.fa.gz Labouyrie_2023_BBS_12000SEQS/Labouyrie_2023_BBS_qm20_OK_12000seqs.fa discarded_12000.txt`

### SUBSAPLE FASTA BASED ON SAMPLE SIZE RANGE (discard samples under the "discard treshold")

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/get_subsampled_FASTA_FOR_GB.py`

`python get_subsampled_FASTA_FOR_GB.py 12000 GB_VOL1_20251903_CLEAN.fa.gz GB_VOL1_20251903_CLEAN_min12k_max15k.fa GB_VOL1_discarded_12k.txt 15000`

### move OTU name in fasta headers from 6th to 1th position:

`sed -E 's/>([^|]+)\|([^|]+)\|([^|]+)\|([^|]+)\|([^|]+)\|(OTU[0-9]+)/>\6|\1|\2|\3|\4|\5/' input.fasta | gzip > output.fasta.gz`

### Chao1

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/get_Chao1_from_commpressed_formate.py`

`python get_Chao1_from_commpressed_formate.py commpressed_otu_tab.samples_to_otus.nosingle_False.txt Labouyrie_2023_BBS_qm20_OK_12000seqs_Chao1.txt`

### FILTER BACTERIA FROM SILVA 138 OUTPUT

`awk -F'\t' '$5+0>=90 && $4+0>=70' Labouyrie_2023_BBS_qm20_OK_12000seqs_MOST_ABUND_SILVA_138_PROCESSED.txt > filtered_output.txt`

`grep '_Bacteria' filtered_output.txt | grep -v ';Chloroplast;' | awk -F '\t' '{print $1}' | awk -F '|' '{print $1}' > BACTERIAL_OTUs.txt`

`zgrep --no-group-separator -A 1 -F -f BACTERIAL_OTUs.txt Labouyrie_2023_BBS_qm20_OK_12000seqs_CLUSTERED.fa.gz`

### GENERATE OTU TABLE

`(head -n 1 Labouyrie_2023_BBS_qm20_FINAL_10000seqs_OTUTABLE.txt && tail -n +2 Labouyrie_2023_BBS_qm20_FINAL_10000seqs_OTUTABLE.txt | sort -k1,1) > Labouyrie_2023_BBS_qm20_FINAL_10000seqs_OTUTABLE_sorted.txt`

### get OTUs with at least 30 samples

```
awk -F'\t' '
NR == 1 {
    # Store header and initialize column counters
    for (i = 2; i <= NF; i++) colSum[i] = 0;
    header = $0;
    next;
}
{
    nonZero = 0;
    for (i = 2; i <= NF; i++) {
        if ($i + 0 > 0) {
            nonZero++;
            colSum[i]++;
        }
    }
    if (nonZero >= 30) {
        lines[NR] = $0;
    }
}
END {
    # Identify non-zero columns
    keepCols[1] = 1; # always keep OTU column
    for (i = 2; i <= length(colSum) + 1; i++) {
        if (colSum[i] > 0) keepCols[i] = 1;
    }

    # Output header with selected columns
    split(header, hFields, "\t");
    outHeader = hFields[1];
    for (i = 2; i <= length(hFields); i++) {
        if (keepCols[i]) outHeader = outHeader "\t" hFields[i];
    }
    print outHeader;

    # Output filtered lines with selected columns
    for (r in lines) {
        split(lines[r], fields, "\t");
        outLine = fields[1];
        for (i = 2; i <= length(fields); i++) {
            if (keepCols[i]) outLine = outLine "\t" fields[i];
        }
        print outLine;
    }
}' GB_VOL1_20251903_CLEAN_10000seqs_CLUSTERED_AND_BINNED_BACONLY_OTUTAB_min30samples_ELIGIBLE.txt > GB_VOL1_20251903_CLEAN_10000seqs_CLUSTERED_AND_BINNED_BACONLY_OTUTAB_min30samples_ELIGIBLE_FINAL.txt
```

