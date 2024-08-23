## this program will extract paired information of cells whose number are larger than knee
library(GNET2)
# generate empty vector based chrom.size
print ('usage= extract_barcode_based_on_knee.R rank_file knee_ratio  cut_rank_file')
argv <- commandArgs(TRUE);
rank_file <- argv[1];
knee_ratio <- as.numeric(argv[2]);
cut_rank_file <- argv[3];

barcode_rank <- read.table(rank_file, quote="\"", comment.char="");
knee_point <- kneepointDetection(barcode_rank$V1[(1:length(barcode_rank$V1))%%100==1])*100;
new_knee_point <- knee_point*knee_ratio;
cut_rank <- barcode_rank[1:new_knee_point,];
write.table(cut_rank, file = cut_rank_file, quote = FALSE, row.names = FALSE, col.names = FALSE,append= FALSE);

