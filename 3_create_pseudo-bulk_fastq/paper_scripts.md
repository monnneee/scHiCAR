### The following pipeline use astrocyte as an example cell type to show how we performed downstream analysis using pseudo-bulk FASTQ files in our paper <mark>(add link)</mark>. 
### This pipeline is adapted from the [nf-core/hicar](https://github.com/jianhong/hicar/tree/dev2rc) pipeline, with modifications: while the original nf-core pipeline only uses R2 reads from long-range (>10kb) read pairs for open chromatin peaks calling, here, we use all R2 reads for peak calling.

To run this pipeline, you need to install the following software:
- **BWA**: [install](https://github.com/lh3/bwa)
- **samtools**: [install](http://www.htslib.org/download/)
- **pairtools**: [install](https://pairtools.readthedocs.io/en/latest/installation.html)
- **macs2**: [install](https://github.com/macs3-project/MACS/wiki/Install-macs2)
- **slopBed**: [install](https://github.com/arq5x/bedtools2/releases/tag/v2.31.0)
- **bedClip and bedGraphToBigWig**: [install](https://github.com/ENCODE-DCC/kentUtils)
- **pairix**: [install](https://github.com/4dn-dcic/pairix)
- **cooler**: [install](https://cooler.readthedocs.io/en/latest/quickstart.html)

#### 1. Align reads to genome and index the sorted BAM File
```
bwa mem -SP -t 12 $BWA_INDEX Astro_R1.fastq.gz Astro_R2.fastq.gz | samtools view -bhS --threads 12 -o Astrocyte.bam -
samtools sort -m 4915M -@ 6 -o Astrocyte.srt.bam -T Astrocyte.srt Astrocyte.bam
samtools index -@ 6 Astrocyte.srt.bam
samtools stats --threads 1 Astrocyte.srt.bam > Astrocyte.stats # collects statistics (optional)
samtools flagstat --threads 1 Astrocyte.srt.bam > Astrocyte.flagstat #counts the number of alignments for each FLAG type (optional)
samtools idxstats --threads 1 Astrocyte.srt.bam > Astrocyte.idxstats #reports alignment summary statistics
```
#### 2. Parse alignments and generate paired end tags (PETs)
```
pairtools parse2 \
        -c genome.fa.sizes \
        --min-mapq 10 --max-insert-size 2000 --max-inter-align-gap 50 --report-position outer \
        --add-pair-index --no-flip --drop-seq --expand --max-expansion-depth 6 \
        --output-stats Astrocyte.pairsam.stat \
        -o Astrocyte.pairsam.gz \
        Astrocyte.bam
```
#### 3. Extract read2 from PETs to call open chromatin peaks with q-value cutoff (0.01)
```
zcat Astrocyte.pairsam.gz | awk 'BEGIN {OFS="\t"};  /^[^#]/ { if ($7 == "+") {$5 = $5 + 4} else if ($7 == "-") {$5 = $5 - 5};  print $4, $5, $5+1, $1, "1", $7}' | \
awk '{gsub(":.*","",$4);print $0}' OFS='\t'|grep -v '!'| \
sort -k1,1 -k2,2n --parallel=20 -T $TMP_DIR |uniq > Astrocyte.R2.ATAC.bed

macs2 callpeak --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits --qval 0.01 --gsize ${genome_size} --format BED --name Astrocyte --treatment Astrocyte.R2.ATAC.bed
```
#### 4. Generate pvalue BIGWIG file for open chromatin visualization
```
sval=$(wc -l Astrocyte.R2.ATAC.bed | awk '{printf "%f", $1/1000000}')

macs2 bdgcmp -t Astrocyte_treat_pileup.bdg -c Astrocyte_control_lambda.bdg --o-prefix Astrocyte -m ppois -S $sval

slopBed -i Astrocyte_ppois.bdg -g genome.fa.sizes -b 0 | bedClip stdin genome.fa.sizes Astrocyte.pval.signal.bedgraph

sort -k1,1 -k2,2n Astrocyte.pval.signal.bedgraph|grep -v chrM > Astrocyte.pval.signal.srt.bedgraph

bedGraphToBigWig Astrocyte.pval.signal.srt.bedgraph genome.fa.sizes Astrocyte_sig.pval.signal.bigwig
```
#### 5. Select uniquely mapped PETs and assign CviQI restriction fragments to PETs
```
pairtools select \
        "(pair_type=='UU') or (pair_type=='UR') or (pair_type=='RU')" \ #unique-unique, unique-rescued, rescued-unique
        -o Astrocyte.selected.pairs.gz \
        --output-rest Astrocyte.unselected.pairs.gz \
        Astrocyte.pairsam.gz

cooler digest -o genome_CviQI.bed genome.fa.sizes genome.fa CviQI

pairtools restrict \
        -f genome_CviQI.bed \
        -o Astrocyte.restrict.pairs.gz \
        Astrocyte.selected.pairs.gz
```
#### 6. Remove PETs with mapped to the same digestion fragment，flip, and deduplicate the rest PETs
```
pairtools select \
        "(COLS[-6]==COLS[-3]) and (chrom1==chrom2)" \
        -o Astrocyte.selected.pairs.gz \
        --output-rest Astrocyte.unselected.pairs.gz \
        Astrocyte.restrict.pairs.gz

pairtools flip \
        -c genome.fa.sizes \
        --nproc-in 2 --nproc-out 2 \
        -o Astrocyte.flip.gz \
        Astrocyte.unselected.pairs.gz
 
pairtools sort \
        --tmpdir ./ \
        --nproc 12 \
        --memory 100G \
        -o Astrocyte.sorted.pairs.gz \
        Astrocyte.flip.gz

pairtools dedup \
        --max-mismatch 1 --method max \
        -o Astrocyte.dedup.pairs.gz \
        --output-stats Astrocyte.dedup.pairs.stat \
        Astrocyte.sorted.pairs.gz

pairix Astrocyte.dedup.pairs.gz #generate index file *.px2
```
### 7. Aggregate PETs into contact matrix in the cooler format (10kb resolution)
```
cooler cload \
        pairix --max-split 2 \
        --nproc 12 \
        genome.fa.sizes:10000 \
        Astrocyte.dedup.pairs.gz \
        Astrocyte.10000.cool

cooler balance --cis-only p 12 Astrocyte.10000.cool #balance matrix

cooler zoomify \ #Generate a multi-resolution cooler file by coarsening.
        --balance -r 10000N \
        -n 12 \
        -o Astrocyte.10000.mcool \
        Astrocyte.10000.cool
```
The *.cool file was used to call interaction peak with [Peakachu](https://github.com/tariks/peakachu). The *.mcool was used to visualization with [Higlass](https://github.com/higlass/higlass)
