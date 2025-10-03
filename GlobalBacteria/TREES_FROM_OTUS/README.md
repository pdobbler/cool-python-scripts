### MULTIPLE ALIGNMENT

`mafft --thread 250 --auto otus_0.1_threshold.fas > otus_0.1_threshold_aligned.fasta`  

`mafft --thread 250 --maxiterate 1000 --globalpair otus_0.1_threshold.fas > otus_0.1_threshold_aligned_precise.fasta`

### CONSTRUCT TREE

`conda activate phyml`

```
iqtree2 -s otus_0.1_threshold_aligned_precise.fasta \
        -m MFP \
        -B 1000 \
        -alrt 1000 \
        -T AUTO > output.txt
```

Analyzing sequences: done in 0.0248304 secs using 2.949e+04% CPU  
   1  OTU000001   34.46%    passed     98.07%  
        ...
  33  OTU000045   34.46%    passed     74.16%  
  34  OTU000046   34.46%    failed      0.06%  
  35  OTU000049   34.46%    passed     99.39%  


`mafft --thread 250 --maxiterate 1000 --globalpair otus_0.1_threshold_GOOD_composition.fas > otus_0.1_threshold_GOOD_composition_aligned_precise.fasta`
