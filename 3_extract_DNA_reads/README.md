````
```
This pipeline was used for extracting DNA reads from same clusters or cell types that identidied by RNA library
```
````

### 1. extract DNA barcodes for each sample and matched with RNA barcodes
```bash
for i in {1..20}  #all samples 
  
do
  
sort -k1b,1 sample${i}_dna_barcode|join -j 1 - ATAC-RNA_barcode.dict|awk '{print"sample'$i'_"$2"\tsample'$i'_"$1}' OFS='\t' > sample${i}_RNA_ATAC.barcode # the 1st column is RNA barcode and 2nd column is matched DNA barcode  
  
done  
  
cat sample*_RNA_ATAC.barcode > total_RNA_ATAC.barcode # merge the total samples together that are used in the 'dna_barcode.R' script  
```
  
### 2. export DNA barcodes of each cell from the same cluster/cell type for each sample
```bash  
Rscript dna_barcode.R # outputs DNA barcodes for each cluster that includes mixed samples  
  
for i in {1..25} # all clusters or  cell types are consistent with the "cluster_list" in "dna_barcode.R"  
  
do
  
perl -i -p -e "s/_/\t/g" cluster${i}/total_dna_barcode # split DNA barcode lines: 1st column is sampleID and 2nd column is DNA barcodes  
  
for j in {1..20}  # all samples  
  
do
  
awk '{if($1=="sample'$j'"){print$2}}' cluster${i}/total_dna_barcode > cluster${i}/sample${j}_barcode.txt  
  
done  
  
done #outputs DNA barcodes for each samples in each cluster  
```
  
### 3. extract DNA reads based on barcodes for each samples in each cluster
```bash
mkdir final_readName  
  
for m in {1..20} # all samples  
  
do
  
zcat DNA_fastq/sample${m}_R1.fastq.gz | awk -F " " '{if(NR%4==1){print $1}}' > sample${m}_DNA_readName # extract all read names from read1 fastq file 

for n in {1..25} # all clusters  
  
do
  
python3 match.py -l sample${m}_DNA_readName -s cluster${n}/sample${m}_barcode.txt -o final_readName/cluster${n}_sample${m}_readName  
  
done  
  
done  
```
  
### 4. generate fastq file for each cluster or cell type
```bash
for i in {1..25} # all clusters  
  
do
  
for j in {1..20}  # all samples  
  
do
  
seqtk subseq DNA_fastq/sample${i}_R1.fastq.gz final_readName/cluster${i}_sample${j}_readName >> cluster${i}_R1.fastq  
  
seqtk subseq DNA_fastq/sample${i}_R2.fastq.gz final_readName/cluster${i}_sample${j}_readName >> cluster${i}_R2.fastq  
  
done  
  
pigz -p 8 cluster${i}_R1.fastq  
  
pigz -p 8 cluster${i}_R2.fastq  
  
done  
```
  
### 5. The above cluster${i}_R*.fastq as pseudo-bulk fastq files were used for running nf-core/hicar (https://github.com/nf-core/hicar)
```bash
nextflow pull jianhong/hicar -r dev2rc #dev2rc is the newest version  
  
nextflow run jianhong/hicar -profile singularity --genome mm10 -r dev2rc --input samplesheet.csv --skip_fastqc --skip_cutadapt --outdir result --skip_interactions --skip_tads --skip_diff_analysis --skip_peak_qc --skip_igv --skip_trackhub --skip_circos --pairtools_parse_version parse2 -resume
```  
