### JGI FUNGAL GENOMES

**blastp**

`diamond blastp -d /mnt/DATA/DATABASES/FUNGAL_PROTEINS_JGI/JGI_FUNGAL_PROTEINS_ANNOTATED_20240403 -q genecalling_fgs.faa -e 1E-5 -o genecalling_JGI_FUN_20240403_best.txt -f 6 -p 256 -b12 -c1`

**PROCESS THE RESULTS**

GET FUNGAL LIST - abreviation vs. organism

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/ANNOTATION/abr_list_2024.txt`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/ANNOTATION/replace_fungal_annot_by_taxname.py`

`python2.7 replace_fungal_annot_by_taxname.py genecalling_JGI_FUN_20240403_best.txt abr_list_2024.txt genecalling_JGI_FUN_best_reformate.txt D`


**UPDATE THE TABLE**


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/ANNOTATION/abr_list_2024.txt`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/ANNOTATION/replace_fungal_annot_by_taxname.py`

`python2.7 replace_fungal_annot_by_taxname.py genecalling_JGI_FUN_20240403_best.txt abr_list_2024.txt genecalling_JGI_FUN_best_reformate_update.txt U`

GET TAXONOMY FOR A NEW HITS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/ANNOTATION/taxonomy_2024.txt`

UPDATE TABLE AND TAX TREE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/ANNOTATION/update_taxonomy_by_fungal_hits.py`

`python2.7 update_taxonomy_by_fungal_hits.py MG_MAPTAB_NORMPERSAMPLE_GENES_BESTTAX_CAZy_KOG_KO.tab 18 22 ga0566199_proteins_nolinebreaksJGI_FUNGAL_PROT_20240403_best_reformate_update.txt taxonomy_2024.txt TAX_tree.tab`


**KOG TREE**

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/ANNOTATION/KOG_TAB_2021_03_03.txt`






