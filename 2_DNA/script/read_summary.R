## generate the statistis for each samples
library(magrittr)
raw = snakemake@input[[1]]
dedup = snakemake@input[[2]]
out = snakemake@output[[1]]

sample_name = stringr::str_split(raw, "\\/", simplify = T) %>% as.vector() %>% .[2]
sample_name = stringr::str_remove(sample_name,'.raw.pairsam.stat')
all_pairs = readr::read_tsv(raw, col_names = F)
all_read = all_pairs$X2[all_pairs$X1=='total']

dep_pairs = readr::read_tsv(dedup, col_names = F)

dep_reads = dep_pairs$X2[dep_pairs$X1=='total_dups']
non_duplicated = dep_pairs$X2[dep_pairs$X1=='total_nodups']
trans = dep_pairs$X2[dep_pairs$X1=='trans']
cis = dep_pairs$X2[dep_pairs$X1=='cis']
longrange =  dep_pairs$X2[dep_pairs$X1=='cis_20kb+']


table = tibble::tribble(
    ~ group, ~ reads,
    'total', all_read,
    'duplicate',dep_reads,
    'non_duplicated',non_duplicated,
    'trans',trans,
    'cis',cis,
    'longRang',longrange
)
table$sample = sample_name

readr::write_tsv(table, path = out)
