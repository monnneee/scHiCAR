# 1.create samples.json file

python3 sample2json.py --fastq_dir fq/

# 2. run the pipeline

snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common,scavenger -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &
