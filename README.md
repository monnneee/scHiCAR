# scHiCAR
Pipelines for scHiCAR data processing
### 1. Process FASTQ files of RNA library with Snakemake
Extract RNA barcodes from the read sequence and add them to the read name. Remove the RNA barcodes from the read sequence. If a read sequence does not contain any RNA barcodes, remove the entire read.
### 2. Process FASTQ files of DNA library with Snakemake
Extract RNA barcodes from the read sequence and add them to the read name. Remove the RNA barcodes from the read sequence. If a read sequence does not contain any RNA barcodes, remove the entire read.
### 3. Generate pseudo-bulk FASTQ files of DNA library
Match the DNA barcodes corresponding to the cells in each cluster or cell type identified from the RNA library. Extract reads with read names containing these DNA barcodes from the processed FASTQ files of the DNA library, and generate a new FASTQ files for each cluster or cell type.
