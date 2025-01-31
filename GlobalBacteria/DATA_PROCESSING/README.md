### GETTING SRA FILES

`fasterq-dump --split-files SRRXXXXXXX`  
`fasterq-dump --split-files --threads 8 SRRXXXXXXX`  

SRRxxxxxxx → NCBI SRA Run (USA)
ERRxxxxxxx → ENA Run (Europe)
DRRxxxxxxx → DDBJ Run (Japan)

### UNZIPPING FILES AND CHECKING

- check zip files content...
echo "FILES" > zip_files.txt
for file in *.zip
do
 unzip -l $file >> zip_files.txt
done

/usr/libexec/p7zip/7z x Egbe_2021_JM.7z
tar -xvf Cruz_2021_UF_demulti.tar.gz

unzip -j archiv.zip -d rozbalene

- create folders
for nazev_souboru in *.zip
do
 IFS="_" read -ra casti <<< "$nazev_souboru"
 slozka="${casti[0]}_${casti[1]}_${casti[2]}"
 mkdir -p "$slozka"
done

- unzip to folder without subfolders
for nazev_souboru in *.zip
do
 IFS="_" read -ra casti <<< "$nazev_souboru"
 slozka="${casti[0]}_${casti[1]}_${casti[2]}"
 unzip -j "$nazev_souboru" -d "$slozka"
done

- unzip to folder without subfolders
unzip -j Lebre_2023_BGN_seq.zip -d Lebre_2023_BGN


mv Lebre_2023_BGN_seq.zip /mnt/DATA/projects/avetrot/RELEASE5/RAW_ZIP_BACKUP/

- types
mv Suetsugu_2021_HR /mnt/DATA1/RELEASE5/0_MANUAL_CHECK_NEEDED
mv Zhuang_2020_MG /mnt/DATA1/RELEASE5/1_GOOD_SINGLE
mv xxx /mnt/DATA1/RELEASE5/2_GOOD_PAIRED

echo "" > info.txt
for d in */ ; do
    echo "$d" >> info.txt
    ls $d | head -4 >> info.txt
done


- print folder content
echo "" > 2_GOOD_PAIRED_studies.txt
for d in */ ; do
    echo "$d" >> 2_GOOD_PAIRED_studies.txt
done

### JOINING - PARALLEL WITH PERCENTAGE OF JOUNED READS IN joining.txt

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

