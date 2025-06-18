# 1. Match each RNA barcode with its corresponding DNA barcode

This step is based on the RNA object generated from the filtered expression matrix using Seurat, along with the predefined 6-bp RNA–DNA barcode pairing scheme used during library preparation.

#### 1.1 Load RNA expression matrix
```
library(Seurat)
data<- Read10X(data.dir = "1_RNA_preprocess/GeneFull/filtered")
rna<-CreateSeuratObject(counts = data, project = "example", min.cells = 3, min.features = 200)
```

#### 1.2 Extract RNA barcode and split into 3 DNA parts
```
df <- data.frame(RNAlibrary = rna$orig.ident)
df$DNAbarcode1<-substr(rownames(df),1,6)
df$DNAbarcode2<-substr(rownames(df),7,12)
df$DNAbarcode3<-substr(rownames(df),13,18)
```

#### 1.3 Replace the first 6bp using predefined map
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

#### 1.4 Reverse complement the last 6bp
```
reverse_complement <- function(seq) {
  complement <- setNames(c("T", "G", "C", "A"), c("A", "C", "G", "T"))
  comp_seq <- unname(complement[unlist(strsplit(seq, ""))])
  rev_seq <- paste(rev(comp_seq), collapse = "")
  return(rev_seq)
}
df$DNAbarcode3<-unname(sapply(df$DNAbarcode3, reverse_complement))
df$DNAbarcode<-paste0(df$DNAbarcode1,df$DNAbarcode2,df$DNAbarcode3) # combine full DNA barcode
df$DNAbarcode1<-NULL
df$DNAbarcode2<-NULL
df$DNAbarcode3<-NULL
```

#### 1.5 filter cells based on the high-quliaty barcode list generate from [3_ATAC_fragment](https://github.com/monnneee/scHiCAR/edit/v2/3_ATAC_fragment) and [4_chromatin_contact](https://github.com/monnneee/scHiCAR/edit/v2/4_chromatin_contact)
```
atac_cut_rank<-read.table("3_ATAC_fragment/02_fragment/*.barcode.cut_rank")$V1
pairs_cut_rank<-read.table("4_chromatin_contact/03_dedup/*.barcode.cut_rank")$V1
valid_DNA_barcodes <- intersect(atac_cut_rank, pairs_cut_rank)
df2 <- df[df$DNAbarcode %in% valid_DNA_barcodes, ]
```

#### 1.6 Subset RNA object and add DNA barcode information to meta.data
```
rna_filter <- subset(rna, cells = rownames(df2))
rna_filter <- AddMetaData(rna_filter, metadata = setNames(df2$DNAbarcode, rownames(df2)), col.name = "DNAbarcode")
```

`write.table(rna_filter@meta.data,"example_metadata.txt",quote=F,sep='\t',col.names=T,row.names=T)` outputs the cell type and DNA barcode information after cell clustering and annotation with Seurat. The metadata table will be used for downstream pseudo-bulk and single-cell analysis of the DNA library.

# 2. Pseudo-bulk analysis

#### 2.1 Generate pseudo-bulk ATAC fragment files for each cell type
```
for i in {celltype1,celltype2,celltype3,...,celltypeN}
do
python3 split_ATAC_fragment.py -l 3_ATAC_fragment/03_filtered/*.filtered.tsv -s ${i}.DNAbarcode -o ${i}.ATAC.fragment.tsv
done
```

#### 2.2 Call open chromatin peaks for each cell types
macs2: [install](https://github.com/macs3-project/MACS/wiki/Install-macs2)
```
genome_size=hs
for i in {celltype1,celltype2,celltype3,...,celltypeN}
do
macs2 callpeak --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits --qval 0.01 --gsize $genome_size --format BED --name ${i} --treatment ${i}.ATAC.fragment.tsv
done
```

#### 2.2 Generate pvalue BIGWIG file for open chromatin visualization
slopBed: [install](https://github.com/arq5x/bedtools2/releases/tag/v2.31.0)

bedClip and bedGraphToBigWig: [install](https://github.com/ENCODE-DCC/kentUtils)

```
chrom_size=The/PATH/of/chrom/sizes/file
for i in {celltype1,celltype2,celltype3,...,celltypeN}
do
sval=$(wc -l ${i}.ATAC.fragment.tsv | awk '{printf "%f", $1/1000000}') #counts the number of tags per million in the (compressed) BED file
macs2 bdgcmp -t ${i}_treat_pileup.bdg -c ${i}_control_lambda.bdg --o-prefix ${i} -m ppois -S $sval
slopBed -i ${i}_ppois.bdg -g $chrom_size -b 0 | bedClip stdin $chrom_size ${i}.pval.signal.bedgraph
grep -E 'chrX|chrY' ${i}.pval.signal.bedgraph > temp
grep -E -v 'chrX|chrY' ${i}.pval.signal.bedgraph|sort -k1.4n -k 2,2n|cat - temp1 > ${i}.pval.signal.srt.bedgraph
bedGraphToBigWig ${i}.pval.signal.srt.bedgraph $chrom_size ${i}_sig.pval.signal.bigwig
rm -f ${i}.pval.signal.bedgraph ${i}.pval.signal.srt.bedgraph
rm temp
done

