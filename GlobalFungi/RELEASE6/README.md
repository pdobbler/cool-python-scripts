### Annotate new variants

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/dereplicate_to_md5_gz.py`

`python2.7 dereplicate_to_md5_gz.py itsx_REL6_ITSX_all_parts_ITS1_final.fa.gz itsx_REL6_ITSX_all_parts_ITS1_md5_variants.fa itsx_REL6_ITSX_all_parts_ITS1_derep.map`


580970909 sequences loaded correctly - 0 sequnces are empty - Omitted!
580970909 sequence variants sorted ...
Dereplication is done - 76579647 groups from 580970909 seqs

### ANNOTATE VARIANTS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/split_fasta_by_group_size.py`


```
python2.7 split_fasta_by_group_size.py GB_BOTH_VOL_20260413_97_clustered_SEEDs.fa 350800
```



```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

for file in *.fas
do
 echo "blastn -query ${file} -db greenegenes2_2024_09 -outfmt 6 -evalue 1E-5 -num_threads 2 -max_target_seqs 10 | sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge > ${file%%.fas}_best.tab"
done > blast_and_sort_command.sh

mkdir -p /mnt/DATA1/tmp
export TMPDIR=/mnt/DATA1/tmp
cat blast_and_sort_command.sh | parallel --tmpdir /mnt/DATA1/tmp
```
