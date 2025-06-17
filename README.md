# scHiCAR
## Pipelines for scHiCAR data processing

### Tools requirement

- **Python**: python3.7 or later
- **snakemake**:  `pip install snakemake==5.13.0`
- **cutadapt**: `pip install cutadapt==3.3`
- **STAR**: [v2.7.5c](https://github.com/alexdobin/STAR/releases/tag/2.7.5c)
- **BWA**: [v0.7.17](https://github.com/lh3/bwa/releases/tag/v0.7.17)
- **SAMTOOLS**: [v1.12](https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2)

### 1. Preproces raw FASTQ files of the RNA library with Snakemake and align sequences to the genome ([code](https://github.com/monnneee/scHiCAR/tree/v2/1_RNA_preprocess))

#### a. Snakemake procedures:
- Trim specific sequences at the 5′ end of Read 1.  
- Extract RNA barcodes from Read 1 and append them to its 5′ end.  
- Remove the template-switching oligo (TSO) from the 5′ end of Read 2.  
- Remove poly(A) tails and adaptor sequences from the 3′ end of Read 2.  
- Split FASTQ files into two subsets based on primer strategy: oligo-dT vs. random hexamer.  
- Extract and count all barcodes from the dataset.  
- Correct barcodes that have only one mismatch relative to the whitelist.  
- Compress and merge filtered FASTQ files from the two primer strategies.

The resulting files (`03_corrected_fq/*_all_*.fastq.gz`) are ready for alignment using STAR.

#### b. Generate filtered gene expression matrices (`barcodes.tsv`, `features.tsv`, and `matrix.mtx`) with STAR.

### 2. Preproces raw FASTQ files of the DNA library with Snakemake ([code](https://github.com/monnneee/scHiCAR/tree/v2/2_DNA_preprocess))

#### Snakemake procedures:
- Trim specific sequences at the 5′ end of both Read 1 and Read 2.  
- Extract barcodes from the read sequences and append them to the read names following the `@` symbol.  
- Generate a list of all extracted barcodes and count their occurrences.  
- Compare extracted barcodes against a provided whitelist.  
- Correct barcodes that have only one mismatch relative to the whitelist.  
- Compress the filtered FASTQ files.  
- Remove ME (mosaic end) sequences from the reads.  

The resulting files (`05_cutME_fq/*_cutME_*.fastq.gz`) are ready for generating ATAC fragment files and chromatin contact pair files.

### 3. Generate ATAC fragment files with Snakemake (`*.tsv.gz`)([code](https://github.com/monnneee/scHiCAR/tree/v2/3_ATAC_fragment))

#### Snakemake procedures:
- Align R1 reads to the reference genome using Snaptools with BWA, and sort BAM files by read name.
- Convert BAM files to fragment-level BED format.
- Extract high-quality cell barcodes based on the knee point of the barcode rank curve.

### 4. Generate chromatin contact pair files with Snakemake (`*.dedup.pairs.gz`)([code](https://github.com/monnneee/scHiCAR/tree/v2/4_chromatin_contact))

### 5. Downsteam pseudo-bulk / single-cell analysis ([code](https://github.com/monnneee/scHiCAR/tree/v2/))
