### MULTIPLE ALIGNMENT

`mafft --thread 250 --auto otus_0.1_threshold.fas > otus_0.1_threshold_aligned.fasta`  

`mafft --thread 250 --maxiterate 1000 --globalpair otus_0.1_threshold.fas > otus_0.1_threshold_aligned_precise.fasta`

### CONSTRUCT TREE

`conda activate phyml`

```
phyml -i otus_0.1_threshold_aligned_precise.fasta \
      -d nt \
      --sms BEST \
      -b 1000 \
      -o tlr
```
