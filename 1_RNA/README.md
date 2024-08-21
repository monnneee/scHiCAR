### 1. move your raw fastq files to 'fq' folder

Make sure that the files have the suffix  `_R1_001.fastq.gz`  or `_R2_001.fastq.gz`.

### 2.create samples.json file

`python3 sample2json.py --fastq_dir fq`

### 3. run the pipeline

`snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common,scavenger -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &`
