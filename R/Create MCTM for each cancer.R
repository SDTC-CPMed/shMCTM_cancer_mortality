# Create MCTM for each cancer 
library(reticulate)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(nichenetr)
library(tidyverse)
library(stringr)
library(readxl)
library(plyr)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(RColorBrewer)

# 
wd = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-Personal/Github/SCBD'
input.dir=paste0(wd,"/Data")
outdir= paste0(wd,"/Output/MCTM")
if (dir.exists(outdir)==F){dir.create(outdir)
  print("outdir created")}

#########################
# Run Nichenet ####
#Prepare Nichenet files
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds")) #ligand_target_matrix
ligand_target_matrix[1:5,1:5] 
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
head(weighted_networks$lr_sig)

DEG_method = "MAST"
cancer_name = c("Lung","Liver","Colon","Ovary","Breast")

for (cancer in cancer_name){
  #1.import DEGs
  DEGs = read.table(file = paste0(wd,"/Data/DEGs_merged_",cancer,".txt"), header= T, sep="\t", stringsAsFactor = T)
  DEGs = DEGs[, apply(DEGs, 2, function(y) any(!is.na(y)))]
  head(DEGs)

  #2.import background genes
  BGs = read.table(file = paste0(wd,"/Data/background_genes_",cancer,".csv"), header= T, sep=",", stringsAsFactor = T)
  head(BGs)
  BGs = BGs[,2:dim(BGs)[2]]

  # 3.setup for nichenet
  outputdir <- paste0(outdir,"/",cancer)
  if (dir.exists(outputdir)==F){dir.create(outputdir)
    print("outputdir created")}

  all_expressed_genes <- BGs
  colnames(all_expressed_genes)
  dim(all_expressed_genes)

  DEGs_human <- DEGs
  colnames(DEGs_human)
  DEGs_human[1:5,1:7]

  # 4. run nichnet
  library(nichenetr) 
  library(Matrix)
  library(dplyr)
  
  all_ligand_activity_MCTM <- data.frame(matrix(ncol=6,nrow=0))
  colnames(all_ligand_activity_MCTM) <- c("test_ligand","auroc","aupr","pearson", "Sender", "Target")
  
  best_upstream_ligands_all_interaction <- data.frame(matrix(ncol=3,nrow=0))
  colnames(best_upstream_ligands_all_interaction) <- c("test_ligand", "Sender", "Target")
  
  for(i in 1:dim(DEGs_human)[2]){ 
    
    for(j in 1:dim(DEGs_human)[2]){

      sender_cluster = colnames(DEGs_human)[i]
      receiver_cluster = colnames(DEGs_human)[j]   
      
      if (receiver_cluster %in% colnames(DEGs_human) ==F){next}
      if (receiver_cluster %in% colnames(all_expressed_genes) ==F){next}
      
      expressed_genes_sender = as.vector(DEGs_human[,i])
      expressed_genes_receiver = as.vector(all_expressed_genes[,receiver_cluster])
      
      # gene set of interest - set DEG as gene set of interest ####
      geneset_oi <- as.vector(DEGs_human[, j] %>% .[. %in% rownames(ligand_target_matrix)] )
      
      # set background expressed genes  ####
      background_expressed_genes = as.vector(expressed_genes_receiver %>%
                                               .[. %in% rownames(ligand_target_matrix)])
      
      # Define a set of potential ligands ####
      ligands = lr_network %>% pull(from) %>% unique()        
      expressed_ligands = intersect(ligands,expressed_genes_sender)  
      
      if (length(expressed_ligands)==0) {
        next
      }
      
      receptors = lr_network %>% pull(to) %>% unique()
      expressed_receptors = intersect(receptors,expressed_genes_receiver)
      if (length(expressed_receptors)==0) {
        next
      }
      
      lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)      #选取ligand 和receptor都要有表达的
      head(lr_network_expressed)
      
      potential_ligands = lr_network_expressed %>% pull(from) %>% unique()  #只有lr_network在当前数据中有表达的ligands才考虑作为potentially active ligands for the NicheNet analysis
      head(potential_ligands)
      
      if (length(potential_ligands)[1]==0){
        next}
      
      ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                                    background_expressed_genes = background_expressed_genes, 
                                                    ligand_target_matrix = ligand_target_matrix, 
                                                    potential_ligands = potential_ligands)
      
      ligand_activities %>% arrange(-pearson)
      best_upstream_ligands_interaction <- ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
      
      best_upstream_ligands_interaction <- cbind(best_upstream_ligands_interaction, Sender=sender_cluster, Target=receiver_cluster)
      best_upstream_ligands_all_interaction <- rbind(best_upstream_ligands_all_interaction, best_upstream_ligands_interaction)
      
      ligand_activities <- cbind(ligand_activities,Sender=sender_cluster, Target=receiver_cluster)
      all_ligand_activity_MCTM <- rbind(all_ligand_activity_MCTM,ligand_activities)
      
      ####Find DSs
      best_upstream_ligands = ligand_activities %>% filter(pearson > 0) %>% pull(test_ligand) %>% unique()
      active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix) %>% bind_rows() %>% drop_na()
      if (dim(active_ligand_target_links_df)[1]==0){
        next
      }
      
      active_ligand_target_links_df = cbind(active_ligand_target_links_df,sender=sender_cluster,receiver=receiver_cluster)
      active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix)
      
      dir.create(paste0(outputdir,"/DS"))
      write.csv(active_ligand_target_links_df,paste0(outputdir,"/DS/active_ligand_target_links_df_",sender_cluster,"_",receiver_cluster,".csv"),row.names = T, col.names = T, sep = ',',quote = F)
      write.csv(active_ligand_target_links,paste0(outputdir,"/DS/active_ligand_target_links_",sender_cluster,"_",receiver_cluster,".csv"),row.names = T, col.names = T, sep = ',',quote = F)
    }
  }

  write.table(all_ligand_activity_MCTM,
              file=paste0(outputdir,"/all_ligand_activity.csv"),
              sep=",", quote=F, col.names=T,row.names=F)
  
  # positive PCC only####
  #  get PCC positive ligands
  PCC_positive_ligand  = all_ligand_activity_MCTM %>% filter(pearson > 0)
  length(unique(PCC_positive_ligand$test_ligand))
  dim(PCC_positive_ligand)
  write.table(PCC_positive_ligand,
              file=paste0(outputdir,"/PCC_positive_ligand.txt"),
              sep="\t", quote=F, col.names=T,row.names=F)
  PCC_positive_ligand = read.table(file =paste0(outputdir,"/PCC_positive_ligand.txt"), header= T, sep="\t",stringsAsFactor = F)
  
  # calculate the most frequent ligands from all clusters
  PCC_positive_ligand_freq_ranked <- data.frame(table(PCC_positive_ligand$test_ligand))
  PCC_positive_ligand_freq_ranked <- PCC_positive_ligand_freq_ranked %>% arrange(-Freq)

  write.table(PCC_positive_ligand_freq_ranked,
              file=paste0(outputdir,"/PCC_positive_ligand_freq_ranked.txt"),
              sep="\t", quote=F, col.names=T,row.names=F)

  # make cluster interaction matrix - for cytoscape
  PCC_positive_cluster_intera_ye <- data.frame(table(PCC_positive_ligand$Sender,PCC_positive_ligand$Target))
  colnames(PCC_positive_cluster_intera_ye) = c("Sender", "Target", "No.Interactions")

  write.table(PCC_positive_cluster_intera_ye,
              file=paste0(outputdir,"/cluster_interaction_PCC_positive.txt"),
              sep="\t", quote=F, col.names=T,row.names=F)
  dim(PCC_positive_cluster_intera_ye)

  positive_ligands_counts = matrix(NA,1,4)
  rownames(positive_ligands_counts) = cancer
  colnames(positive_ligands_counts) = c("N_positive_ligands","N_all_ligands","N_positive_interactions","N_all_interactions")
  positive_ligands_counts[1,1] = length(unique(PCC_positive_ligand$test_ligand))
  positive_ligands_counts[1,2] = length(unique(all_ligand_activity_MCTM$test_ligand))
  positive_ligands_counts[1,3] = dim(PCC_positive_ligand)[1]
  positive_ligands_counts[1,4] = dim(all_ligand_activity_MCTM)[1]
  write.table(positive_ligands_counts,
              file=paste0(outputdir,"/positive_ligands_counts.txt"),
              sep="\t", quote=F, col.names=T,row.names=T)

  # Centrality ####
  library(igraph)
  library(CINNA)
  
  all_ligand_activity_MCTM <- PCC_positive_ligand
  ligands_interaction <- all_ligand_activity_MCTM # cell-cell interactions + ligand information
  ligands_interaction <- as.matrix(ligands_interaction)
  head(ligands_interaction)
  
  # Create graph
  g <- graph_from_edgelist(el = ligands_interaction[,6:7], directed = T)
  
  out <- matrix(NA, nrow = 4, ncol = length(V(g)$name))
  
  rownames(out) <- c("node_degree_all", "node_degree_in", "node_degree_out", "closeness")
  colnames(out) <- paste("Cluster_", V(g)$name, sep="")
  
  # Calculate centrality degree
  out[1,] <- centr_degree(g, mode = "all")$res # all
  out[2,] <- centr_degree(g, mode = "in")$res # in degree
  out[3,] <- centr_degree(g, mode = "out")$res # out degree
  
  # Calculate closeness
  out[4,] <- closeness(g, mode = "all")
  
  # Calculate centralities
  #pr_cent<-proper_centralities(g)
  #calculate_centralities(g, include = pr_cent[1:20])  %>% pca_centralities(scale.unit = TRUE)
  proper_centralities(g)
  centrality_matrix <- calculate_centralities(g, include = c("eigenvector centralities", "Kleinberg's hub centrality scores", 
                                                             "Laplacian Centrality", "Leverage Centrality"))
  p = pca_centralities( centrality_matrix  ) #figure out the order of most important centrality types based on your graph structure
  ggsave(plot=p,  path = outputdir, dpi = 300, width = 25, height = 20, units = "cm" ,filename="PCA_centrality_matrix.pdf")
  
  tsne_centralities( centrality_matrix, dims = 2, perplexity = 1, scale=TRUE) #Another method for distinguishing which centrality measure has more information or in another words has more costs is using (t-SNE) t-Distributed Stochastic Neighbor Embedding analysis(Van Der Maaten 2014). 
  p = visualize_heatmap( centrality_matrix , scale = TRUE  )
  ggsave(plot=p,  path = outputdir, dpi = 300, width = 25, height = 20, units = "cm" ,filename="centrality_matrix.pdf")
  
  
  centrality_matrix <- matrix(unlist(centrality_matrix), ncol = length(V(g)$name), byrow = T)
  colnames(centrality_matrix) <- paste("Cluster_", V(g)$name, sep="")
  rownames(centrality_matrix) <- c("eigenvector centralities", "Kleinberg's hub centrality scores", 
                                   "Laplacian Centrality", "Leverage Centrality")
  
  out <- rbind(out, centrality_matrix)
  
  write.table(out, file=paste0(outputdir,"/Cell-cell_centrality_summary_PCC_positive.txt"), sep="\t", col.names = NA, row.names = T)
  
  print(paste0(cancer, "NicheNet done!"))
}











