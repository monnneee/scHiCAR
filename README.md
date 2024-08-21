# scHiCAR
## Pipelines for scHiCAR data processing

### Requirements

To run this pipeline, you need to install the following software:

- **Python3**
- **snakemake**:  `pip install snakemake==5.13.0`
- **cutadapt**: `pip install cutadapt==3.3`
- **STAR**: https://github.com/alexdobin/STAR/releases/tag/2.7.5c
  
### 1. Process FASTQ files of RNA library with Snakemake (requires downloading contents of the '1_RNA' folder)
a. Extract ***RNA barcodes*** from the read sequence and add them to the beginning of read 1 (*_R1_001.fastq). Remove the adaptors from the read sequence. If a read sequence does not contain any RNA barcodes, remove the entire read.

b. genrate filtered matrix (`barcodes.tsv`, `features.tsv`, and `matrix.mtx`) with `STARsolo`

### 2. Process FASTQ files of DNA library with Snakemake (requires downloading contents of the '2_DNA' folder)
Extract ***DNA barcodes*** from the read sequence and add them to the read name. Remove the adaptors from the read sequence. If a read sequence does not contain any DNA barcodes, remove the entire read.
### 3. Generate pseudo-bulk FASTQ files of DNA library
Match the ***DNA barcodes*** corresponding to the cells in each cluster or cell type identified from the RNA library. Extract reads with read names containing these DNA barcodes from the processed FASTQ files of the DNA library, and generate a new FASTQ files for each cluster or cell type.
