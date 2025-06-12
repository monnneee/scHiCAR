# scHiCAR
## Pipelines for scHiCAR data processing

### Tools requirement

- **Python**: python3.7 or later
- **snakemake**:  `pip install snakemake==5.13.0`
- **cutadapt**: `pip install cutadapt==3.3`
- **STAR**: [v2.7.10a](https://github.com/alexdobin/STAR/releases/tag/2.7.10a)

### 1. Process RNA FASTQ files with Snakemake and align sequences to the genome ([1_RNA](https://github.com/monnneee/scHiCAR/tree/v2/1_RNA))

    1.1. Snakemake procedures:
    - Trim specific sequences at the 5′ end of Read 1.  
    - Extract RNA barcodes from Read 1 and append them to its 5′ end.  
    - Remove the template-switching oligo (TSO) from the 5′ end of Read 2.  
    - Remove poly(A) tails and adaptor sequences from the 3′ end of Read 2.  
    - Split FASTQ files into two subsets based on priming strategy: oligo-dT vs. random hexamer.  
    - Extract and count all barcodes from the dataset.  
    - Compare extracted barcodes against a provided whitelist and correct those with a single mismatch.  
    - Compress and merge filtered FASTQ files from the two primer strategies.

    The resulting files (`03_corrected_fq/*_all_*.fastq.gz`) are ready for alignment using STAR.

    1.2. Generate filtered gene expression matrices (`barcodes.tsv`, `features.tsv`, and `matrix.mtx`) using STAR, suitable for downstream scRNA-seq clustering analysis.

### 2. Process raw FASTQ files of DNA library with Snakemake ([2_DNA](https://github.com/monnneee/scHiCAR/tree/main/2_DNA))
a. Extract ***DNA barcodes*** from the read sequence and add them to the read name. Remove the adaptors from the read sequence. If a read sequence does not contain any DNA barcodes, remove the entire read.

b. Generate a filtered list of DNA barcodes by removing background noise

### 3. Generate pseudo-bulk FASTQ files of DNA library ([3_create_pseudo-bulk_fastq](https://github.com/monnneee/scHiCAR/tree/main/3_create_pseudo-bulk_fastq))
Match the ***DNA barcodes*** corresponding to the cells in each cluster or cell type identified from the RNA library. Extract reads with read names containing these DNA barcodes from the processed FASTQ files of the DNA library, and generate a new FASTQ files for each cluster or cell type.

The downstream scripts for pseudo-bulk FASTQ files used in scHiCAR paper refer to the [nf-core/hicar](https://github.com/jianhong/hicar/tree/dev2rc) pipeline, the details can be found in [paper_scripts.md](https://github.com/monnneee/scHiCAR/blob/main/3_create_pseudo-bulk_fastq/paper_scripts.md)
