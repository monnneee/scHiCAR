### 1.extract DNA barcodes for each sample and matched with RNA barcodes
for i in {sample1,sample2,...,sample20} #modified based on your sample names
do
sort -k1b,1 ${i}_dna_barcode|join -j 1 - ATAC-RNA_barcode.dict|awk '{print"'$i'_"$2"\t'$i'_"$1}' OFS='\t' > ${i}_RNA_ATAC.barcode # the 1st column is RNA barcode and 2nd column is matched DNA barcode
done
cat sample*_RNA_ATAC.barcode > total_RNA_ATAC.barcode # merge total samples together
### 2. export DNA barcodes of each cell from the same cluster/cell type for each sample
Rscript dna_barcode.R
for i in {}
cluster_list[[i]],"/total_dna_barcode
### 3. extract DNA reads based on barcodes
zcat DNA_fastq/sample_R1.fastq.gz | awk -F ' ' '{if(NR%4==1){print $1}}' > DNA_fastq.readName #extract all read names from read1 fastq file of DNA library
ls *_dna_barcode 
for i in {}
do
python3 match.py -l DNA_fastq.readName
-s ${i}_dna_barcode  -o ${i}.readName

