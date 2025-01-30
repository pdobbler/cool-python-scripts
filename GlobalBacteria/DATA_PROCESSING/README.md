### GETTING SRA FILES

fasterq-dump --split-files SRRXXXXXXX
fasterq-dump --split-files --threads 8 SRRXXXXXXX

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

### JOINING

```
output_file="joining.txt" > "$output_file"
for file in *_1.fastq
do
  sample=${file%%_1.fastq}
  echo -n "$sample;" >> "$output_file"
  line_count=$(fastq-join -v " " -p 15 -m 40 ${sample}_1.fastq  ${sample}_2.fastq -o ${sample}_joined | awk '{ORS=";"; print} END {print ""}')
  echo $line_count >> "$output_file"
done
```

```
output_file="joining.txt"
> "$output_file"

export output_file

parallel --bar --keep-order '
  sample={1}
  echo -n "$sample;" >> "$output_file"
  line_count=$(fastq-join -v " " -p 15 -m 40 ${sample}_1.fastq ${sample}_2.fastq -o ${sample}_joined | awk "{ORS=\";\"; print} END {print \"\"}")
  echo $line_count >> "$output_file"
' ::: $(ls *_1.fastq | sed 's/_1.fastq//')
```



