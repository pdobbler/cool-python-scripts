## Amplicon data

### GET THE DATA

`wget -O AMPLICON_DATA.zip "http://www.biomed.cas.cz/mbu/lbwrf/example/AMPLICON_DATA.zip"`

```
unzip AMPLICON_DATA.zip
cd AMPLICON_DATA/
```

### OTU BASED - CLUSTERING APPROACH

`cd 16S_example_data/`

**PREPARATION OF SAMPLE TABLE FROM FWD and REV tags**
```
head samples.txt
head barcoded_primers_fwd.txt
head barcoded_primers_rev.txt
```

**be sure about formate**
```
dos2unix samples.txt
dos2unix barcoded_primers_fwd.txt
dos2unix barcoded_primers_rev.txt
```

```
awk -F'\t' '{print $2}' samples.txt > fwd_tags.txt
awk -F'\t' '{print $3}' samples.txt > rev_tags.txt
```
