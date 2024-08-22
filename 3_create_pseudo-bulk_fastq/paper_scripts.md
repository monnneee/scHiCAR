### The following scripts use Astrocytes as an example to demonstrate how we performed downstream analysis using pseudo-bulk FASTQ files in our scHiCAR paper. 
### These scripts are adapted from the nf-core/hicar pipeline (https://github.com/jianhong/hicar/tree/dev2rc), with modifications: while the original nf-core pipeline only uses R2 reads from long-range (>10kb) read pairs for open chromatin peaks calling, here, we use all R2 reads for peak calling.

```
bwa mem -SP -t 12 $BWA_INDEX Astro_R1.fastq.gz Astro_R2.fastq.gz | samtools view -bhS --threads 12 -o Astrocyte_REP1_T1.bam -
```
