# Amplicon data

## GET THE DATA

`wget -O AMPLICON_DATA.zip "http://www.biomed.cas.cz/mbu/lbwrf/example/AMPLICON_DATA.zip"`

```
unzip AMPLICON_DATA.zip
cd AMPLICON_DATA/
```

## OTU BASED - CLUSTERING APPROACH

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

### CREATE AWK SCRIPT

`nano vlookup.awk`

```
# vlookup script
FNR==NR{
  a[$1]=$col
  next
}
{OFS="\t"; if ($1 in a) {print a[$1], $1} else {print "NA", $1}  }
```

**Ctrl+o & Enter (save)**

**Ctrl+x         (exit)**



**run the script**
```
awk -v col='3' -f vlookup.awk barcoded_primers_fwd.txt fwd_tags.txt > cols_fwd.txt
awk -v col='3' -f vlookup.awk barcoded_primers_rev.txt rev_tags.txt > cols_rev.txt
awk -F'\t' '{print $1}' samples.txt > cols_sample.txt
```

**function "paste" joins the columns of the table**

```
paste -d'\t' cols_sample.txt cols_fwd.txt cols_rev.txt > barcode_table.txt
```

**or**

```
paste <(awk -F'\t' '{print $1}' samples.txt) <(awk -v col='3' -f vlookup.awk barcoded_primers_fwd.txt fwd_tags.txt) <(awk -v col='3' -f vlookup.awk barcoded_primers_rev.txt rev_tags.txt)
```

### DEMULTIPLEXING

**move data to a new folder**

```
mkdir demulti
mv barcode_table.txt demulti/
mv *.fastq demulti/
cd demulti/
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/Amplicons/Course_MB140P94/demultiplexing_fastq.py`
