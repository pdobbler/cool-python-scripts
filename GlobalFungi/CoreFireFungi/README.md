### PROCESS SRA


`mkdir SRAs`  
`find ./ -type f -name "*.sra" -exec mv {} SRAs/ \;`  
  
`conda env list`  
`conda activate /mnt1/florian/.conda/envs/sra_env`  

```
for f in *.sra; do
    fasterq-dump --split-files "$f"
done
```
