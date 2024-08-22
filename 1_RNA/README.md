### 1. Download all the files to your folder
```
Your_folder
├── cluster.json
├── config.yaml
├── sample2json.py
├── linear_sciRNA_18bp_barcode.txt.gz
├── Snakefile
├── fq  # move your raw fastq files to this folder
│   ├── RNA_example_R1_001.fastq.gz
│   └── RNA_example_R2_001.fastq.gz
└── script
    ├── barcode_hash_v2.py
    ├── fq_barcode_correction_R1.py
```

### 2.create samples.json file

`python3 sample2json.py --fastq_dir fq`

### 3. run the pipeline

`snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common,scavenger -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &`

### 4. align reads to genome and generate filtered matrix (`barcodes.tsv`, `features.tsv`, and `matrix.mtx`) for use in standard scRNA-seq downstream analysis (such as [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial)).
