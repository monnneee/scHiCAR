# scHiCAR
## Pipelines for scHiCAR data processing

### Tools requirement

- **Python**: python3.9 or above
- **snakemake**:  `pip install snakemake==5.13.0`
- **cutadapt**: `pip install cutadapt==3.3`
- **STAR**: [v2.7.5c](https://github.com/alexdobin/STAR/releases/tag/2.7.5c)

### 1. Process FASTQ files of RNA library with Snakemake ([README.md](https://github.com/monnneee/scHiCAR/tree/main/1_RNA))
a. Extract ***RNA barcodes*** from the read sequence and add them to the beginning of read 1 (*_R1_001.fastq). Remove the adaptors from the read sequence. If a read sequence does not contain any RNA barcodes, remove the entire read.

b. genrate filtered matrix (`barcodes.tsv`, `features.tsv`, and `matrix.mtx`) for use in standard scRNA-seq downstream analysis.

### 2. Process FASTQ files of DNA library with Snakemake ([README.md](https://github.com/monnneee/scHiCAR/tree/main/2_DNA))
a. Extract ***DNA barcodes*** from the read sequence and add them to the read name. Remove the adaptors from the read sequence. If a read sequence does not contain any DNA barcodes, remove the entire read.

b. Generate a filtered list of DNA barcodes by removing background noise

### 3. Generate pseudo-bulk FASTQ files of DNA library ([README.md](https://github.com/monnneee/scHiCAR/tree/main/3_create_pseudo-bulk_fastq))
Match the ***DNA barcodes*** corresponding to the cells in each cluster or cell type identified from the RNA library. Extract reads with read names containing these DNA barcodes from the processed FASTQ files of the DNA library, and generate a new FASTQ files for each cluster or cell type.
