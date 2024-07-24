### 1. export DNA barcodes of each cell in the same cluster/cell type
Rscript dna_barcode.R
### 2. extract DNA reads based on barcodes
zcat DNA_fastq/sample_R1.fastq.gz | awk -F ' ' '{if(NR%4==1){print $1}}' > DNA_fastq.readName #extract all read names from read1 fastq file of DNA library
ls *_dna_barcode 
for i in {}
do
python3 match.py -l DNA_fastq.readName
-s ${i}_dna_barcode  -o ${i}.readName

