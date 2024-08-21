### 1. Download all the files to your folder
```
Your_folder
├── cluster.json
├── config.yaml
├── sample2json.py
├── linear_sciRNA_18bp_barcode.txt.gz
├── Snakefile
├── fq         # move your raw fastq files to this folder
│   ├── DNA_example_R1_001.fastq.gz
│   └── DNA_example_R2_001.fastq.gz
└── script
    ├── barcode_hash_v2.py
    ├── fq_barcode_correction.py  
    ├── raw_fq_update.py
    └── sample2json.py
```

### 2.create samples.json file

`python3 sample2json.py --fastq_dir fq`

### 3. run the pipeline

`snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common,scavenger -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &`
