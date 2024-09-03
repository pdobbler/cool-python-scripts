### EXTRACT ITS

**`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/WHOLE_ITS/test_dereplicated.fa`**

`/home/kdanielmorais/bioinformatics/tools/ITSx_1.0.11/ITSx -i test_dereplicated.fa --cpu 2 --only_full T -t F -o itsx_test_dereplicated`


### EXTRACT COMPLETE ITS

**`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/WHOLE_ITS/get_whole_its_region.py`**


`python2.7 get_whole_its_region.py test_dereplicated.fa itsx_test_dereplicated.positions.txt test_dereplicated_complete_ITS.fa`

