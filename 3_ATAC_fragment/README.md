## Note: In the latest scHiCAR protocol, read 1 corresponds to the ATAC-seq fragment.
### 1. Download all the files to your folder
```
Your_folder
├── Snakefile
├── cluster.json
├── config.yaml #modify this file based on your file paths
├── sample2json.py
├── fq  # This folder contains processed R1 FASTQ files (e.g., *_cutME_L001_R1_001.fastq.gz)
│   ├── DNA_example_cutME_L001_R1_001.fastq.gz
└── script
    ├── scHiCAR_R2_parse.py
```

### 2.Create samples.json file

`python3 sample2json.py --fastq_dir fq` or `python3 sample2json.py --fastq_dir ../2_DNA_preprocess/05_cutME_fq`

### 3. Run snakemake pipeline (customize -p as needed based on your HPC environment)

`snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &`
