### 1. Download all the files from this folder


### 2. move your raw fastq files to 'fq' folder

Make sure that the files have a suffix like `_R1_001.fastq.gz`  or `_R2_001.fastq.gz`.

### 3.create samples.json file

`python3 sample2json.py --fastq_dir fq`

### 4. run the pipeline

`snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common,scavenger -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &`
