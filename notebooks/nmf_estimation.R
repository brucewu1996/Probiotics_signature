library(NMF)
library(ggplot2)

lacto_df = read.table('/home/bruce1996/repo/Microbiome_health_indicator/tutorial/data/lacto_abundance_matrix.txt',sep = '\t',header = T,row.names = 1)
bifido_df = read.table('/home/bruce1996/repo/Microbiome_health_indicator/tutorial/data/bifido_abundance_matrix.txt',sep = '\t',header = T,row.names = 1)

ranks <- 2:7
bifido_mat <- as.matrix(bifido_df)
i0 <- which(colSums(bifido_mat) == 0)
i_na <- which(colSums(is.na(bifido_mat)) > 0)
nmf_input = bifido_mat[, -c(col_0, col_na)] + 10 ** -8
bifido_estim.coad <- nmf(nmf_input,ranks,nrun = 5,.opt='v')
bifido_p = plot(bifido_estim.coad) + ggtitle('Clustering evaluation of Bifidobacterium')
ggsave('/home/bruce1996/repo/Microbiome_health_indicator/tutorial/nmf_evaluation/bifido_nmf_evaluation.png',dpi =300,width = 6,height = 4)

lacto_mat <- as.matrix(lacto_df)
col_0 <- which(colSums(lacto_mat) == 0)
col_na <- which(colSums(is.na(lacto_mat)) > 0)
nmf_input = lacto_mat[, -c(col_0, col_na)] + 10 ** -8
lacto_estim.coad <- nmf(nmf_input,ranks,'lee',nrun = 5,.opt='v')
lacto_p = plot(lacto_estim.coad) + ggtitle('Clustering evaluation of Lactobacillus subtype')
ggsave('/home/bruce1996/repo/Microbiome_health_indicator/tutorial/nmf_evaluation/lacto_nmf_evaluation.png',dpi =300,width = 6,height = 4)


