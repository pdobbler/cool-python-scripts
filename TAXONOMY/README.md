## GET TAXONOMY FROM NCBI

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/TAXONOMY/retrieve_taxonomy_from_ncbi_optimal.py`

**Based on NCBI accessions**

`python2.7 retrieve_taxonomy_from_ncbi_optimal.py missing_acc.txt A missing_acc_taxonomy.txt`

**Based on taxon names**

`python2.7 retrieve_taxonomy_from_ncbi_optimal.py missing_names.txt N missing_names_taxonomy.txt`

## PREPARE NVBInr AND TAXONOMY TABLE

*get NCBInr database:: (diamond can also use the blast.db files. if needed can be used instead of building a new one)*
`wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz`


*get prot.accession2taxid.FULL file*
`wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz`

*get taxdump from NCBI*
`wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`

*build diamond database with taxonomy information*

`diamond makedb --in nr.gz --db NCBInr_diamond_datedownoload --taxonmap prot.accession2taxid.FULL.gz --taxonnodes nodes.dmp --threads 200`

*activate/create environment, here is ncbi-taxonomist*
`conda activate ncbi-taxonomist`

*get taxonkit https://bioinf.shenwei.me/taxonkit/usage/#lineage*
`conda install -c bioconda taxonkit`

*export correct taxdump file, --db flag not always work*
`export TAXONKIT_DB='/mnt/DATA/DATABASES/NCBI_nr_dmnd_09_2022/taxdump_new_2024/'`

*get taxid list from 'acc2taxID' file*
`zcat prot.accession2taxid.FULL.gz |cut -f2 | sort | uniq > all_taxid_052024.txt`

*get lineages from taxid_list*
`cat all_taxid_052024.txt | taxonkit lineage > ncbi_052024_lineages_update.txt`

*formating flat tsv table*
`cat ncbi_052024_lineages_update.txt | taxonkit reformat --output-ambiguous-result -r Unassigned -F -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" > taxonomy_ncbi_052024_update.txt`

*remove extra 2nd columns, put TaxID last between []*
`awk -F'\t' '{print $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t" "[" $1 "]"}' taxonomy_ncbi_052024_update.txt > taxonomy_ncbi_052024__update.txt`

*needs to remove first and last lines and colunm names*
`cat  taxonomy_ncbi_052024__update.txt | sed '1,1d' | sed '$d' | sed '1i domain\tphylum\tclass\torder\tfamily\tgenus\torganism\ttax_key'> taxonomy_ncbi_052024__update__done.txt`

*fix cases where the genera from species name have brakets []?*
OBS.. Square brackets ([ ]) around a genus indicates that the name awaits appropriate action by the research community to be transferred to another genus.
