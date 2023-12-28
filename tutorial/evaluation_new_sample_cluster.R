#new 
library(readxl)
library(ConsensusClusterPlus)

get_consensus_matrix <- function(consensus_clustering_result,consensus_number){
  # Extract the consensus matrix from ConsensusClusteringPlus output.
  # 
  # Args:
  #   consensus_clustering_result (list): The ConsensusClusteringPlus output.
  # consensus_number (int): The number of cluster.
  m = consensus_clustering_result[[consensus_number]]$consensusMatrix
  samples = names(consensus_clustering_result[[consensus_number]]$consensusClass)
  df = as.matrix(m)
  rownames(df) = samples
  colnames(df) = samples
  return(df)
}

determine_new_sample_cluster <- function(exp_m, metadata, cluster_colnames='cluster', new_sample_label='None') {
  # Evaluation new sample cluster based on previous clustering result.
  # 
  # Args:
  #   exp_m (matrix): The origin matrix for ConsensusClusteringPlus.
  # metadata (data.frame): The data.frame including the patient metadata.
  # cluster_colnames (str, optional): The colname of clustering result. Defaults to 'cluster'.
  # new_sample_label (str, optional): The element in clustering result column for new / non-clustered samples. Defaults to 'None'.
  
  cluster_label <- unique(metadata[[cluster_colnames]])  # Unique cluster labels are cluster1, 2, 3, 4 + none
  cluster_label <- setdiff(cluster_label, new_sample_label)
  consensus_res = ConsensusClusterPlus(exp_m,maxK=10,reps=50,pItem=0.8,pFeature=1,
                             clusterAlg="hc",
                             distance = 'euclidean',
                             seed=1262118388.71279,
                             plot="png")
  consensus_matrix <- get_consensus_matrix(consensus_res,length(cluster_label))
  #subset the sample list
  new_samples = rownames(metadata)[metadata$cluster == new_sample_label]
  #Create a blank consensus similarity matrix
  cluster_consensus_sim <- matrix(0, nrow = length(new_samples), ncol = length(cluster_label))
  rownames(cluster_consensus_sim) <- new_samples
  colnames(cluster_consensus_sim) <- cluster_label
  
  for (sample in new_samples) {
    for (cluster in cluster_label) {
      cluster_sample_list <- rownames(metadata)[metadata[[cluster_colnames]] == cluster]
      # Calculate the average consensus for new sample to cluster_sample_list (samples in specific cluster)
      cluster_consensus_sim[sample, cluster] <- mean(consensus_matrix[sample, cluster_sample_list])
    }
  }
  return(cluster_consensus_sim)
}
#example
repo_dir = "/home/bruce1996/repo/Microbiome_health_indicator/"
exp_m = read.table(paste0(repo_dir,"tutorial/sig_matrix/sig_proprotion_matrix.txt"),
                   header = T,row.names = 1,sep = '\t',encoding = "UTF-8")
tmp = read_excel(paste0(repo_dir,"tutorial/data/TPMIC_Diagnosis_297_KCF_0704.xlsx"))
metadata = as.data.frame(tmp)
rownames(metadata) = metadata$ID
metadata = metadata[colnames(exp_m),]
exp_m = as.matrix(exp_m)
##########
metadata$cluster = sample(c('a','b','c','d'),dim(metadata)[1],replace=TRUE)
cluster_colnames = 'cluster'
new_sample_label = 'd'
res = determine_new_sample_cluster(exp_m = exp_m ,metadata = metadata,new_sample_label = 'd')

consensus_res = ConsensusClusterPlus(exp_m,maxK=10,reps=50,pItem=0.8,pFeature=1,
                                     clusterAlg="hc",
                                     distance = 'euclidean',
                                     seed=1262118388.71279,
                                     plot="png")
