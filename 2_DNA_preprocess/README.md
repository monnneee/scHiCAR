### 1. Download all the files to your folder
```
Your_folder
├── ME_index
├── Snakefile
├── cluster.json
├── sample2json.py
├── scHiCAR_DNA_18bp_barcode.txt.gz
├── fq  # move your raw fastq files to this folder
│   ├── DNA_example_R1_001.fastq.gz
│   └── DNA_example_R2_001.fastq.gz
└── script
    ├── barcode_hash_v2_ME.py
    ├── fq_barcode_correction_R1_ME.py
    ├── raw_fq_update.py
```

### 2.Create samples.json file

`python3 sample2json.py --fastq_dir fq`

### 3. Run snakemake pipeline (customize -p as needed based on your HPC environment)
Before running Snakemake, please make sure all required Python packages used in the `.py` files under the `script` folder are installed.

`snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &`
