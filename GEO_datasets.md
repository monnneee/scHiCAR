# scHiCAR Processed Data

We have uploaded the processed scHiCAR data to [GSE305889](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305889).

The processed data from the scHiCAR datasets used in our [Nature Biotechnology paper](https://www.nature.com/articles/s41587-026-03013-7) are described below.

## 1. cellline_scHiCAR

The human cell line scHiCAR dataset, [GSE267117](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267117), corresponds to **Fig. 1** of our Nature Biotechnology paper.

The following processed files are available:

* **`.pairsam` file**: Read 2 represents the ATAC-seq reads, while paired Read 1 and Read 2 represent the two ends of each chromatin contact.
* **RNA matrix**: Gene expression matrix.

> **Note:** For the cell line datasets, the RNA and DNA data were generated from **two different batches**. The two libraries were evaluated separately to assess the quality and advantages of the RNA and DNA libraries of scHiCAR, and were not matched at the individual cell level.

## 2. mouse_brain_scHiCAR_2

[GSE267126](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267126)

This dataset corresponds to **Figs. 3–5** of our Nature Biotechnology paper.

For each cell type, the following processed files are available:

* **`.pairsam` file**
* **RNA matrix**
* **Metadata files**

### Extracting ATAC Fragments and Chromatin Contact Pairs

For both the human cell line and mouse brain datasets described above, ATAC fragments and chromatin contact pairs can be extracted from the `.pairsam` files by following **Steps 3–6** in our processing pipeline:

[scHiCAR processing pipeline – Step 3: Extract Read 2 from PETs to call open chromatin peaks](https://github.com/DiaoLab/scHiCAR/blob/main/3_create_pseudo-bulk_fastq/paper_scripts.md#3-extract-read2-from-pets-to-call-open-chromatin-peaks-with-q-value-cutoff-001)

## 3. mouse_brain_scHiCAR_1

[GSE305439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305439)

This dataset generated using our **latest protocol**, corresponds to **Fig. 2** of our Nature Biotechnology paper.

Unlike the earlier protocol, the latest protocol captures **ATAC-seq information in Read 1**.

The following processed files are available:

* **ATAC fragment files**
* **Chromatin contact pair files**
* **RNA matrix files**
* **Metadata files** 

## 4. mouse_skeletal_muscle_scHiCAR

[GSE304674](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304674)
This dataset was also generated using our **latest protocol** and corresponds to **Fig. 6** of our Nature Biotechnology paper.

The following processed files are available:

* **ATAC fragment files**
* **Chromatin contact pair files**
* **RNA matrix files**
* **Metadata files**

## Summary of scHiCAR Datasets

| Dataset                       | GEO Accession                                                             | Figure    | Protocol         | Processed Files                                     |
| ----------------------------- | ------------------------------------------------------------------------- | --------- | ---------------- | --------------------------------------------------- |
| Cell line scHiCAR                         | [GSE267117](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267117) | Fig. 1    | Earlier protocol | `.pairsam`, RNA matrix, metadata                    |
| Mouse brain scHiCAR (1.62 million cells)  | [GSE267126](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267126) | Figs. 3-5 | Earlier protocol | `.pairsam`, RNA matrix, metadata                    |
| Mouse brain scHiCAR (5,313 cells)         | [GSE305439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305439) | Fig. 2    | Latest protocol  | ATAC fragments, contact pairs, RNA matrix, metadata |
| Mouse skeletal muscle scHiCAR             | [GSE304674](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304674) | Fig. 6    | Latest protocol  | ATAC fragments, contact pairs, RNA matrix, metadata |

