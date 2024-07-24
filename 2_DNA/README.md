python3 sample2json.py --fastq_dir fq # ./fq contains your original *fastq.gz  
  
snakemake --latency-wait 100 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common,scavenger -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &  
