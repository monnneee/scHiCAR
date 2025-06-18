### 1. Match each RNA barcode with its corresponding DNA barcode

This step is based on the RNA object generated from the filtered expression matrix using Seurat, along with the predefined 6-bp RNA–DNA barcode pairing scheme used during library preparation.
```
library(Seurat)
data<- Read10X(data.dir = "1_RNA_preprocess/GeneFull/filtered")
rna<-CreateSeuratObject(counts = data, project = "xw346", min.cells = 3, min.features = 200)
df<-as.data.frame(rna$orig.ident)
colnames(df)<-"RNAlibrary"
df$DNAbarcode1<-substr(gsub("^(.*)_","",rownames(df)),1,6)
df$DNAbarcode2<-substr(gsub("^(.*)_","",rownames(df)),7,12)
df$DNAbarcode3<-substr(gsub("^(.*)_","",rownames(df)),13,18)
```

#### 1.1 replace the first 6 bp
```
barcode_map <- c( #predefined 6-bp RNA–DNA barcode pairing scheme used during library preparation.
  "TCATCC" = "TACCCG",
  "AGTCAA" = "GAGTTT",
        ...
  "ATGTCA" = "CGATTT",
  "CCCTAC" = "CTGGAG"
)
df$DNAbarcode1 <- sapply(df$DNAbarcode1, function(barcode) {
  for (old_prefix in names(barcode_map)) {
    if (startsWith(barcode, old_prefix)) {
      return(sub(paste0("^", old_prefix), barcode_map[old_prefix], barcode))
    }
  }
  return(barcode)
})
```

#### 1.2 reverse complement the last 6bp
```
reverse_complement <- function(seq) {
  complement <- setNames(c("T", "G", "C", "A"), c("A", "C", "G", "T"))
  comp_seq <- unname(complement[unlist(strsplit(seq, ""))])
  rev_seq <- paste(rev(comp_seq), collapse = "")
  return(rev_seq)
}
df$DNAbarcode3<-unname(sapply(df$DNAbarcode3, reverse_complement))
df$DNAbarcode<-paste(df$DNAbarcode1,df$DNAbarcode2,df$DNAbarcode3,sep="")
df$DNAbarcode1<-NULL
df$DNAbarcode2<-NULL
df$DNAbarcode3<-NULL
```

#### 1.3 filter cells based on the high-quliaty barcode list generate from [3_ATAC_fragment](https://github.com/monnneee/scHiCAR/edit/v2/3_ATAC_fragment) and [4_chromatin_contact](https://github.com/monnneee/scHiCAR/edit/v2/4_chromatin_contact)
```
atac_cut_rank<-read.table("3_ATAC_fragment/02_fragment/*.barcode.cut_rank")$V1
pairs_cut_rank<-read.table("4_chromatin_contact/03_dedup/*.barcode.cut_rank")$V1
valid_DNA_barcodes <- intersect(atac_cut_rank, pairs_cut_rank)
df2 <- df[df$DNAbarcode %in% valid_DNA_barcodes, ]
rna$temp<-names(rna$orig.ident)
rna_filter<-subset(rna, temp %in% rownames(df2))
rna_filter$temp<-NULL
rna_filter_DNAbarcode<-df2$DNAbarcode
names(rna_filter_DNAbarcode)<-rownames(df2)
rna_filter<-AddMetaData(object = rna_filter,metadata = rna_filter_DNAbarcode,col.name = 'DNAbarcode')
```
