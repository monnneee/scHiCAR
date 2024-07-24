```R
### 1. add DNA barcodes to seurat object
library(Seurat)
rna_seurat_object<-readRDS("rna_seurat.rds")
dna_barcode<-read.table("./total_RNA_DNA_barcode.txt",sep='\t',header=F,row.names=1) #1st column is RNA barcode and 2nd column is DNA barcode
dna_bd<-dna_barcode$V2
names(dna_bd)<-rownames(dna_barcode)
rna_seurat_object$temp<-names(rna_seurat_object$orig.ident)
rna_filtered_object<-subset(rna_seurat_object,temp %in% rownames(dna_barcode)) #filter out cells without matched DNA barcodes
rna_filtered_object$temp<-NULL
rna_filtered_object <- AddMetaData(object = rna_filtered_object, metadata = dna_bd,col.name = 'dna_barcode')
```

### 2. extract DNA barcodes after identified clusters or annotated cell types:
```R
cluster<-names(table(rna_filtered_object$cluster)) #get total cluster name
cluster_list<-as.list(cluster)
for (i in 1:length(cluster)){
mkdir cluster_list[[i]]
bd<-subset(rna_filtered_object,cluster %in% cluster_list[[i]])$dna_barcode
names(bd)<-NULL
write.table(bd,paste(cluster_list[[i]],"/total_dna_barcode",sep=""),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
} #each cluster has an output file named "total_dna_barcode" that includes mixed samples
```
