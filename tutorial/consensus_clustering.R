library(ConsensusClusterPlus)

consensus_matrix <- function(consensus_res,consensus_number){
  m = consensus_res[[consensus_number]]$consensusMatrix
  samples = names(res[[4]]$consensusClass)
  df = as.data.frame(m)
  rownames(df) = samples
  colnames(df) = samples
  return(df)
}

probiotic_consensus <- function(exp_m,n_target,consensus_path,output_path,max_k =10){
  # exp_m : matrix; row is feature, col is sample
  # n_target : integer; number of consensus clustering
  # max_k : integer; max number of cluster number
  # consensus_path : output path of consensus clustering result
  dir.create(consensus_path, showWarnings = FALSE)
  
  results = ConsensusClusterPlus(exp_m,maxK=max_k,reps=50,pItem=0.8,pFeature=1,
                                 title=consensus_path,clusterAlg="pam",distance="euclidean",seed=1262118388.71279,plot='png')
  label = data.frame(cluster = results[[n_target]]$consensusClass)
  write.table(label,output_path,sep='\t',quote = F)
}

consensus_p = "/home/bruce1996/repo/Microbiome_health_indicator/tutorial/consensus_evaluation"
output_p = "/home/bruce1996/repo/Microbiome_health_indicator/tutorial/consensus_clustering/consensus_label.txt"

exp_m = read.table("/home/bruce1996/repo/Microbiome_health_indicator/tutorial/sig_matrix/sig_proprotion_matrix.txt",
                   header = T,row.names = 1,sep = '\t',encoding = "UTF-8")
exp_m = as.matrix(exp_m)
res = ConsensusClusterPlus(exp_m,maxK=10,reps=50,pItem=0.8,pFeature=1,
                           clusterAlg="hc",title=consensus_p,
                           distance = 'euclidean',
                           seed=1262118388.71279,
                           plot="png")
n_cluster = 4
df = as.data.frame(res[[n_cluster]]$consensusClass)
colnames(df) = c('cluster')
cm = consensus_matrix(res,n_cluster)
write.table(cm,file =output_p,sep = '\t',quote = F)






