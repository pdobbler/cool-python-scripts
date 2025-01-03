
### get SH community matrix

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/OTU_and_SH_tables/TABLE_BY_SAMPLES_AND_SHS.py`


- get complete table (WARNING - could be very large):

`python2.7 TABLE_BY_SAMPLES_AND_SHS.py GF5_RAW_TABLE_TAB_SH.txt.gz GF5_RAW_TABLE_SAMPLES.txt.gz - -`

- get all SH based on samples:

`python2.7 TABLE_BY_SAMPLES_AND_SHS.py GF5_RAW_TABLE_TAB_SH.txt.gz GF5_RAW_TABLE_SAMPLES.txt.gz test_samples.txt -`

- get all samples based on SHs:

`python2.7 TABLE_BY_SAMPLES_AND_SHS.py GF5_RAW_TABLE_TAB_SH.txt.gz GF5_RAW_TABLE_SAMPLES.txt.gz - test_shs.txt`

- get selection of samples and SH:

`python2.7 TABLE_BY_SAMPLES_AND_SHS.py GF5_RAW_TABLE_TAB_SH.txt.gz GF5_RAW_TABLE_SAMPLES.txt.gz test_samples.txt test_shs.txt`


### GET OTUTABLE FROM CLUSTERED FASTA

header e.g.: >GF05023557S|Fernan_2020_JZ12|e7c017be559f77ec|OTU0027655

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/OTU_and_SH_tables/GET_OTUTAB_FROM_SEQUENCES.py`

`python2.7 GET_OTUTAB_FROM_SEQUENCES.py GF5_ALL_SAMPLES_ITS2_CLUSTERED_ECTO.fa GF5_ALL_SAMPLES_ITS2_CLUSTERED_ECTO_OTUTABLE.txt`


### GET REDUCED OTUTABLE FROM CLUSTERED FASTA

header e.g.: >GF05028373S|Chen12_2022_AKD1|c4327a8b28d158a0|OTU0096999

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/OTU_and_SH_tables/GET_OTUTAB_FROM_SEQUENCES_REDUCED_FORMATE.py`

all clusters and singletons included (no_singletons:false;-:means no selected OTUs list provided)

`python2.7 GET_OTUTAB_FROM_SEQUENCES_REDUCED_FORMATE.py GF5_ALL_SAMPLES_ITS2_CLUSTERED.fa.gz false -`

### AWK script for selecting columns with at least 40 non-zero values

e.g.  
| Samples     | SH0737578.10FU | SH0759515.10FU | SH1049840.10FU | SH1019594.10FU | SH0880862.10FU |
|-------------|---------------:|---------------:|---------------:|---------------:|---------------:|
| GF01002699S |             0  |             0  |             0  |             0  |             8  |
| GF05008849S |             0  |             0  |             0  |             0  |             0  |
| GF05000787S |            20  |             0  |             0  |            78  |             0  |
| GF04000029S |            20  |             4  |             0  |             0  |            88  |
| GF04002266S |             0  |             6  |             0  |             0  |            88  |
| GF01016797S |             0  |             0  |            40  |             4  |             8  |
| GF05000787V |            15  |             9  |             0  |             0  |             8  |
| GF04010997S |             0  |             9  |            40  |             4  |             0  |
| GF03009554S |             0  |             0  |             0  |             4  |             0  |

```
awk -F'\t' '
{
    # For each column, we check for non-zero values
    if (NR > 1) {
        for (i = 2; i <= NF; i++) {
            if ($i != "" && $i != "0") count[i]++;
        }
    } else {
        # Uložíme hlavičku
        for (i = 1; i <= NF; i++) header[i] = $i;
    }
}
END {
    # We print the corresponding rows
    while ((getline line < FILENAME) > 0) {
        split(line, fields, FS);
        first = 1;
        for (i = 1; i <= NF; i++) {
            if (i == 1 || count[i] > 40) {
                printf "%s%s", (first ? "" : OFS), fields[i];
                first = 0;
            }
        }
        printf "\n";
    }
}
' OFS='\t' GF5_SH_IN_SAMPLES.relative.txt > GF5_SH_IN_SAMPLES_40moreSamples.relative.txt
```


