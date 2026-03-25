
### get primers taxonomy coverage 

`wget https://www.biomed.cas.cz/mbu/lbwrf/seed/archive/SILVA138_SEED2_20230818.zip`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PRIMERS_COVERAGE/test_primer_pair_coverage.py`

`python2.7 test_primer_pair_coverage.py silva_138_1_fixed_taxa_MithoChloro_FINAL_20230818.fasta 27F/1492R AGAGTTTGATCMTGGCTCAG TACGGYTACCTTGTTAYGACTT`

### get coverage table

`mv PRIMERS_COVERAGE`

`mv *_taxonomy_coverage.txt PRIMERS_COVERAGE/`

`python3 merge_taxonomy.py`
