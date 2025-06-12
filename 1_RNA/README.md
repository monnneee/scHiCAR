### 1. Download all the files to your folder
```
Your_folder
├── cluster.json
├── sample2json.py
├── sciHiCAR_RNA_18bp_barcode.txt.gz
├── Snakefile
├── ME_index
├── fq  # move your raw fastq files to this folder
│   ├── RNA_example_R1_001.fastq.gz
│   └── RNA_example_R2_001.fastq.gz
└── script
    ├── barcode_hash_v2.py
    ├── fq_barcode_correction_R1.py
```

### 2.Create samples.json file

`python3 sample2json.py --fastq_dir fq`

### 3. Run snakemake pipeline （customize -p as needed based on your HPC environment）

`snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &`

### 4. Align reads to the genome and generate a filtered matrix folder that includes the files `barcodes.tsv`, `features.tsv`, and `matrix.mtx`
```
gunzip -c sciHiCAR_RNA_18bp_barcode.txt.gz > sciHiCAR_RNA_18bp_barcode.txt

STAR --runMode alignReads \
--genomeDir ./GRCm38_STAR_2.7.6a \
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
--soloUMIlen 12 \
--soloCBmatchWLtype Exact \
--soloUMIdedup 1MM_CR \
--soloStrand Forward \
--soloUMIfiltering - \
--readFilesIn snakemake_output/03_corrected_fq/RNA_example_L001_R2_001.fastq.gz snakemake_output/03_corrected_fq/RNA_example_L001_R1_001.fastq.gz \
--readFilesCommand zcat \
--genomeSAindexNbases 2 \
--soloBarcodeReadLength 0 \
--soloCellFilter EmptyDrops_CR > Result/Log/sciHiCAR-RNA-2-230124-230303_L001_STAR.log 2>&1
```
The STAR output **GeneFull/filtered** folder can be used in standard scRNA-seq downstream analysis (such as with [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial)).
