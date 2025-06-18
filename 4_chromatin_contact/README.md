### 1. Download all the files to your folder
```
Your_folder
├── Snakefile
├── cluster.json
├── config.yaml
├── sample2json.py
├── fq    # This folder contains preprocessed FASTQ files
│   ├── DNA_example_cutME_L001_R1_001.fastq.gz
│   └── DNA_example_cutME_L001_R2_001.fastq.gz
└── script
    ├── extract_barcode_based_on_knee.R
    ├── filter_pairs.py
    ├── read_summary.R
```

### 2.Create samples.json file

`python3 sample2json.py --fastq_dir fq`

### 3. Run snakemake pipeline (customize -p as needed based on your HPC environment)
Before running Snakemake, please make sure all required R/Python packages used in the `.R`/`.py` file under the `script` folder are installed.

`snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &`

The output `05_filtered/*.dedup.filtered.pairs.gz` files can be used in downstream pseudo-bulk or single-cell analysis.
