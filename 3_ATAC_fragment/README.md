## Note: In the latest scHiCAR protocol, read 1 corresponds to the ATAC-seq fragment.
### 1. Download all the files to your folder
```
Your_folder
├── Snakefile
├── cluster.json
├── config.yaml #modify this file based on your file paths
├── sample2json.py
├── fq  # This folder contains preprocessed R1 FASTQ files (e.g., *_cutME_L001_R1_001.fastq.gz)
│   ├── DNA_example_cutME_L001_R1_001.fastq.gz
└── script
    ├── extract_barcode_based_on_knee.R
    ├── scHiCAR_R2_parse.py
```

### 2.Create samples.json file

`python3 sample2json.py --fastq_dir fq`

### 3. Run snakemake pipeline (customize -p as needed based on your HPC environment)
Before running Snakemake, please make sure all required R/Python packages used in the `.R`/`.py` file under the `script` folder are installed.

`snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &`

### Output files:
- The **02_fragment/*.tsv.gz** files can be used in standard scATAC-seq downstream analysis (such as with [ArchR]([https://www.archrproject.com/bookdown/creating-arrow-files.html](https://www.archrproject.com/articles/Articles/tutorial.html)), starting with **reformatFragmentFiles** to prepare for **createArrowFiles**).
- The **02_fragment/*.cut_rank** files list high-quality barcodes used to filter the final fragment files.
