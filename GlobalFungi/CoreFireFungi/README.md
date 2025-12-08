### PROCESS SRA


`mkdir SRAs`  
`find ./ -type f -name "*.sra" -exec mv {} SRAs/ \;`  
```
for f in *.sra; do
    fasterq-dump --split-files "$f"
done
```
