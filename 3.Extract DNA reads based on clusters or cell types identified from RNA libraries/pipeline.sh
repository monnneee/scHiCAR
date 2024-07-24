### 1.extract DNA barcodes for each sample and matched with RNA barcodes
for i in {sample1,sample2,...,sample20} #modified based on your sample names
do
sort -k1b,1 ${i}_dna_barcode|join -j 1 - ATAC-RNA_barcode.dict|awk '{print"'$i'_"$2"\t'$i'_"$1}' OFS='\t' > ${i}_RNA_ATAC.barcode # the 1st column is RNA barcode and 2nd column is matched DNA barcode
done
cat sample*_RNA_ATAC.barcode > total_RNA_ATAC.barcode # merge the total samples together that are used in the 'dna_barcode.R' script

### 2. export DNA barcodes of each cell from the same cluster/cell type for each sample
Rscript dna_barcode.R #export DNA barcodes for each cluster/cell type that includes mixed samples
for i in {sample1,sample2,...,sample20} #all cluster names or  cell types are consistent with the "cluster_list" in "dna_barcode.R"
do
awk '{if($1=="'$i'"){print$2}}' Astro/Astro_barcode.txt > Astro/${i}_Astro_barcode.txt
cluster$i/total_dna_barcode
done

### 3. extract DNA reads based on barcodes
zcat DNA_fastq/sample_R1.fastq.gz | awk -F " " '{if(NR%4==1){print $1}}' > DNA_fastq.readName #extract all read names from read1 fastq file of DNA library
ls *_dna_barcode 
for i in {}
do
python3 match.py -l DNA_fastq.readName
-s ${i}_dna_barcode  -o ${i}.readName

