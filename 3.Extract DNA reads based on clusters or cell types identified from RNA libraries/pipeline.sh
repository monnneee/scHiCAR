### 1.extract DNA barcodes for each sample and matched with RNA barcodes
for i in {sample1,sample2,...,sample20} #modified based on your sample names
do
sort -k1b,1 ${i}_dna_barcode|join -j 1 - ATAC-RNA_barcode.dict|awk '{print"'$i'_"$2"\t'$i'_"$1}' OFS='\t' > ${i}_RNA_ATAC.barcode # the 1st column is RNA barcode and 2nd column is matched DNA barcode
done
cat sample*_RNA_ATAC.barcode > total_RNA_ATAC.barcode # merge the total samples together that are used in the 'dna_barcode.R' script

### 2. export DNA barcodes of each cell from the same cluster/cell type for each sample
Rscript dna_barcode.R #outputs DNA barcodes for each cluster that includes mixed samples
for i in {1..25} #all clusters or  cell types are consistent with the "cluster_list" in "dna_barcode.R"
do
perl -i -p -e "s/_/\t/g" cluster${i}/total_dna_barcode # split DNA barcode lines: 1st column is sampleID and 2nd column is DNA barcodes
for j in {1..20}  #all samples
do
awk '{if($1=="sample'$j'"){print$2}}' cluster${i}/total_dna_barcode > cluster${i}/sample${j}_barcode.txt
done
done #outputs DNA barcodes for each samples in each cluster

### 3. extract DNA reads based on barcodes for each samples in each cluster
mkdir final_readName
for m in {1..20}  #all samples
do
zcat DNA_fastq/sample${m}_R1.fastq.gz | awk -F " " '{if(NR%4==1){print $1}}' > sample${m}_DNA_readName #extract all read names from read1 fastq file
for n in {1..25} #all clusters
do
python3 match.py -l sample${m}_DNA_readName -s cluster${n}/sample${m}_barcode.txt -o final_readName/cluster${n}_sample${m}_readName
done
done

### 4. generate fastq file for each cluster or cell type
for i in {}
do
python3 match.py -l DNA_fastq.readName
-s ${i}_dna_barcode  -o ${i}.readName

