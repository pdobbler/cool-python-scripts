### GETTING SRA FILES

`fasterq-dump --split-files SRRXXXXXXX`  
`fasterq-dump --split-files --threads 8 SRRXXXXXXX`  

SRRxxxxxxx → NCBI SRA Run (USA)  
ERRxxxxxxx → ENA Run (Europe)  
DRRxxxxxxx → DDBJ Run (Japan)  

### GETTING SPLITTED FILES

```
for file in *.fastq; do
    # Extract the base name without extension
    base_name="${file%.fastq}"
    
    # Run the fasterq-dump command
    echo "Processing: $base_name"
    fasterq-dump --split-files "$base_name"
done
```


### UNZIPPING FILES AND CHECKING

unzipping the file

`unzip Frey_2021_QQ_Bact_seq.zip`

try repairing the ZIP File  (warning [Frey_2021_QQ_Bact_seq.zip]:  2824149584 extra bytes at beginning or within zipfile)

`zip -FF Frey_2021_QQ_Bact_seq.zip --out Frey_2021_QQ_Bact_seq_fixed.zip`

- check zip files content...
```
echo "FILES" > zip_files.txt
for file in *.zip
do
 unzip -l $file >> zip_files.txt
done
```

`/usr/libexec/p7zip/7z x Egbe_2021_JM.7z`

`tar -xvf Cruz_2021_UF_demulti.tar.gz`

`unzip -j archiv.zip -d rozbalene`

- create folders
```
for nazev_souboru in *.zip
do
 IFS="_" read -ra casti <<< "$nazev_souboru"
 slozka="${casti[0]}_${casti[1]}_${casti[2]}"
 mkdir -p "$slozka"
done
```

- unzip to folder without subfolders
```
for nazev_souboru in *.zip
do
 IFS="_" read -ra casti <<< "$nazev_souboru"
 slozka="${casti[0]}_${casti[1]}_${casti[2]}"
 unzip -j "$nazev_souboru" -d "$slozka"
done
```

- unzip to folder without subfolders
`unzip -j Lebre_2023_BGN_seq.zip -d Lebre_2023_BGN`

`mv Lebre_2023_BGN_seq.zip /mnt/DATA/projects/avetrot/RELEASE5/RAW_ZIP_BACKUP/`

- types
`mv Suetsugu_2021_HR /mnt/DATA1/RELEASE5/0_MANUAL_CHECK_NEEDED`
`mv Zhuang_2020_MG /mnt/DATA1/RELEASE5/1_GOOD_SINGLE`
`mv xxx /mnt/DATA1/RELEASE5/2_GOOD_PAIRED`

```
echo "" > info.txt
for d in */ ; do
    echo "$d" >> info.txt
    ls $d | head -4 >> info.txt
done
```

- print folder content

```
echo "" > 2_GOOD_PAIRED_studies.txt
for d in */ ; do
    echo "$d" >> 2_GOOD_PAIRED_studies.txt
done
```

### JOINING - PARALLEL WITH PERCENTAGE OF JOUNED READS IN joining.txt

checking reads length  

`wc -L SRR11880035_2.fastq`

joining  

```
output_file="joining.txt"
> "$output_file"

export output_file

parallel --bar --keep-order '
  sample={1}
  file1="${sample}_1.fastq"
  file2="${sample}_2.fastq"
  line_count=$(fastq-join -v " " -p 15 -m 40 "$file1" "$file2" -o "${sample}_joined" | awk "{ORS=\";\"; print} END {print \"\"}")

  total_reads=$(echo "$line_count" | grep -oP "Total reads: \K[0-9]+")
  total_joined=$(echo "$line_count" | grep -oP "Total joined: \K[0-9]+")

  if [[ -n "$total_reads" && -n "$total_joined" && "$total_reads" -gt 0 ]]; then
    percentage=$(awk "BEGIN {printf \"%.2f\", ($total_joined / $total_reads) * 100}")
  else
    percentage="N/A"
  fi

  echo "$sample ($file1, $file2);$line_count $percentage%" >> "$output_file"
' ::: $(ls *_1.fastq | sed 's/_1.fastq//')

```

### CHECK joining.txt

if joining is <90% (e.g.:70%) try to cut 300 bp to 250 bp...  

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/resize_fastq_length_GZIP.py`

```
for file in *.fastq
do
 sample=${file%%.fastq}
 echo "python2.7 resize_fastq_length_GZIP.py ${file} 250_${sample}.fastq 250"
done > cut_fastq.sh

cat cut_fastq.sh | parallel
```


### QUALITY FILTERING

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/filter_fastq_by_quality_mean_GZIP.py`

```
mkdir QM20

for file in *_joinedjoin
do
 sample=${file%%_joinedjoin}
 echo "python2.7 filter_fastq_by_quality_mean_GZIP.py ${file} QM20/${sample}_qm20.fa 20 True"
done > filter_fastq.sh

cat filter_fastq.sh | parallel
```

### RENAMING

```
mkdir RENAME
cat QM20/*_qm20.fa.gz > RENAME/Study_XXXX_YYYY_qm20.fa.gz
cd RENAME
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/find_and_replace_or_append_text_beginning.py`

`python2.7 find_and_replace_or_append_text_beginning.py Study_XXXX_YYYY_qm20.fa.gz samples.txt Study_XXXX_YYYY_qm20_renamed.fas true`

This should have same counts:  
`grep '>' Study_XXXX_YYYY_qm20_renamed.fas | wc -l`  
`grep '>GB' Study_XXXX_YYYY_qm20_renamed.fas | wc -l`

### RENAMING BY APPENDING NAME

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/add_name_to_titles_FASTA.py`

`python2.7 add_name_to_titles_FASTA.py 1_qm20.fa.gz 'GB01016506S|Delgado_2018_1aae_Bact'`

---

### SEARCH FOR PRIMARY MOTIVE


        Updated sequences: 515F (Parada)–806R (Apprill), forward-barcoded:
        FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
        Original sequences: 515F (Caporaso)–806R (Caporaso), reverse-barcoded:
        FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/search_for_primary_motive_with_reverse.py`

`python2.7 search_for_primary_motive_with_reverse.py Clavel_2021_1GT_qm20_renamed.fa GTGYCAGCMGCCGCGGTAA 2`


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/search_for_primary_motive_with_reverse_cut_primary.py`

`python2.7 search_for_primary_motive_with_reverse_cut_primary.py Clavel_2021_1GT_qm20_renamed.fa GTGYCAGCMGCCGCGGTAA 4`

515F

```
MOTIVE="GTGYCAGCMGCCGCGGTAA"
MISMATCHES=4

# Run the Python script in parallel for all .fas.gz files
ls *_renamed.fas.gz | parallel -j $(nproc) "python2.7 search_for_primary_motive_with_reverse_cut_primary.py {} $MOTIVE $MISMATCHES"
```

519F

```
MOTIVE="CAGCMGCCGCGGTAA"
MISMATCHES=3

# Run the Python script in parallel for all .fas.gz files
ls *_renamed.fas.gz | parallel -j $(nproc) "python2.7 search_for_primary_motive_with_reverse_cut_primary.py {} $MOTIVE $MISMATCHES"
```

### INVESTIGATING FILES

```
#!/bin/bash

# Output file
OUTPUT="ALIGNMENT.fa"
> "$OUTPUT"  # Truncate or create output

for file in *PRIMARY.fa.gz; do
    # Extract prefix before _qm20_renamed.fas.gz.PRIMARY.fa.gz
    prefix="${file%%_qm20_renamed.fas.gz.PRIMARY.fa.gz}"

    # Process 10 sequences quickly and cleanly
    zcat "$file" | awk -v prefix="$prefix" '
        BEGIN { seq_count = 0; printing = 0; }
        /^>/ {
            if (++seq_count > 10) exit;
            header = substr($0, 2);  # strip >
            print ">" prefix "|" header;
            printing = 1;
            next;
        }
        {
            if (printing) print $0;
        }
    ' >> "$OUTPUT"
done
```

RANDOM

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/get_random_fasta_gzip.py`

```
for file in *_qm20_renamed.fas.gz.PRIMARYCUT.fa.gz
do
 sample=${file%%_qm20_renamed.fas.gz.PRIMARYCUT.fa.gz}
 echo "python2.7 get_random_fasta_gzip.py 10 ${file} ${sample}_10rand.fa"
done > random_choose.sh

cat random_choose.sh | parallel
```

### SEARCH FOR SECONDARY MOTIVE

806R - GGACTACHVGGGTWTCTAAT

ATTAGA-ACCC--GTA
mm 4


806R

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/find_and_cut_secondary.py`

```
MOTIVE="ATTAGANACCCNNGTA"
MISMATCHES=4

# Run the Python script in parallel for all .PRIMARYCUT.fa.gz files
ls *.PRIMARYCUT.fa.gz | parallel -j $(nproc) "python2.7 find_and_cut_secondary.py {} $MOTIVE $MISMATCHES"
```

### STATS

parallel "zgrep '>GB' {} | awk -F '|' '{print \$1}' | sort | uniq -c | sed 's/^/{}: /'" ::: *.gz > all_sample_counts.txt

### REMOVE AMBIGUOUS BASES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/FASTA_divide_by_ambiguous_fast.py`

`ls *.gz | parallel -j $(nproc) "python FASTA_divide_by_ambiguous_fast.py {}"`

### FILTER BY LENGTH CUTOFF

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/FASTA_get_lengths_profile_and_divide_by_length.py`

```
MIN=223
MAX=283

ls *_correct.fa.gz | parallel -j $(nproc) "python2.7 FASTA_get_lengths_profile_and_divide_by_length.py {} $MIN $MAX"
```

### MAKE VARIANTS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/FASTA_to_md5_variants_FASTA.py`

`python2.7 FASTA_to_md5_variants_FASTA.py GB_VOL1_20251903_CLEAN.fa.gz`

count singletons

`zgrep ';size=1$' GB_VOL1_20251903_CLEAN.fa.gz_variants.fa.gz | wc -l`


