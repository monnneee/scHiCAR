### 1. Download all the files to your folder
```
Your_folder
├── ME_index
├── Snakefile
├── cluster.json
├── sample2json.py
├── scHiCAR_RNA_18bp_barcode.txt.gz
├── fq  # move your raw fastq files to this folder
│   ├── RNA_example_R1_001.fastq.gz
│   └── RNA_example_R2_001.fastq.gz
└── script
    ├── barcode_hash_v2.py
    ├── fq_barcode_correction_R1.py
```

### 2.Create samples.json file

`python3 sample2json.py --fastq_dir fq`

### 3. Run snakemake pipeline (customize -p as needed based on your HPC environment)
Note: Before running Snakemake, please make sure all required Python packages used in the *.py files under the [script](https://github.com/monnneee/scHiCAR/tree/v2/1_RNA_preprocess/script) folder are installed.

`snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &`

### 4. Align reads to the genome and generate a filtered matrix folder that includes the files `barcodes.tsv`, `features.tsv`, and `matrix.mtx`
```
gunzip -c sciHiCAR_RNA_18bp_barcode.txt.gz > sciHiCAR_RNA_18bp_barcode.txt

STAR --runMode alignReads \
--genomeDir PATH_TO_STAR_INDEX_folder \
--runThreadN 12 \
--outFileNamePrefix RNA_example \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
--soloType CB_UMI_Simple \
--soloFeatures GeneFull \
--soloCBwhitelist sciHiCAR_RNA_18bp_barcode.txt \
--soloCBstart 1 \
--outSAMmapqUnique 255 \
--soloCBlen 18 \
--soloUMIstart 19 \
--soloUMIlen 16 \
--soloCBmatchWLtype Exact \
--soloUMIdedup 1MM_CR \
--soloStrand Forward \
--soloUMIfiltering - \
--readFilesIn 03_corrected_fq/RNA_example_all_L001_R2_001.fastq.gz 03_corrected_fq/RNA_example_all_L001_R1_001.fastq.gz \
--readFilesCommand zcat \
--genomeSAindexNbases 2 \
--soloBarcodeReadLength 0 \
--soloCellFilter EmptyDrops_CR
--limitBAMsortRAM 200000000000 > log 2>&1 &
```
The STAR output **GeneFull/filtered** folder can be used in standard scRNA-seq downstream analysis (such as with [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial)).
