### RANK AND SORT

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PRETZEL_clustering/score_variants.py`

`python2.7 score_variants.py GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz`

Processed sequences 2001205289 to variants 514364177 (GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz)

### CLUSTERING

`/home/pretzel/ali`

```
#  /mnt/DATA1/Align/align512
IjaClust ver 3.2
Usage: align <[dir/]input.fasta[.gz|.gzip]> [similarity_threshold] [-o output_dir] [serial] [-m max_threads]
Example: align data/sequences.fasta.gz 97 -o output -m 8
Threshold can be defined with one decimal place. Default value is 97.0
Output file is written to input or provided directory with name based on input file name and threshold
Max threads can be set to 'serial' for single-threaded execution or a positive integer for multi-threading (default is 200 or number of hardware threads)
```

`/mnt/DATA1/Align/align512 GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.gz 97.0 -m 512`
IjaClust ver 3.2
Path: GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.gz zipped Similarity threshold: 97 Max threads: 512
