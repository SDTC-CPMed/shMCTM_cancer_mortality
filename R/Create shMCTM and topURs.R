# Create shMCTM 
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

wd = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-Personal/Github/SCBD'

input.dir=paste0(wd,"/Data")
outdir= paste0(wd,"/Output/MCTM/shMCTM")
if (dir.exists(outdir)==F){dir.create(outdir)
  print("outdir created")}
cancer_name = c('Lung','Liver','Ovary','Colon','Breast')
#########################
# find shared UR####
ligand_list = list()
for (cancer in cancer_name){
  ligands <- read.csv(paste0(wd,"/Output/MCTM/",cancer,"/PCC_positive_ligand.txt"), sep="\t",header= T)
  ligand_list = c(ligand_list,list(ligands$test_ligand %>% unique ) )
}
names(ligand_list) = cancer_name
ligand_list = c(ligand_list['Breast'],ligand_list['Colon'],ligand_list['Liver'],ligand_list['Lung'],ligand_list['Ovary'])

shared_ligands <- intersect(unlist(ligand_list[1]), unlist(ligand_list[2])) %>% 
  intersect(.,unlist(ligand_list[3])) %>% 
  intersect(.,unlist(ligand_list[4])) %>% 
  intersect(.,unlist(ligand_list[5]))
length(shared_ligands)                     
write.table(shared_ligands, file=paste0(outdir,"/shared_UR.txt"), sep="\t", col.names = F, row.names = F)

# venn plot of all URs
library(VennDiagram)
library(gplots)
library(RColorBrewer)
cancer_color <- c("#ABDDA4","#FFFFBF","#FDAE61","#D7191C","#2B83BA")

venn_list <- ligand_list
venn.plot <- venn.diagram(x=venn_list, 
                          col="black",fill=cancer_color , # #paletteer_d("ggthemes::Classic_Green_Orange_6",5)
                          filename = NULL,
                          cex=1,cat.cex = 1,cat.fontface = "bold",alpha=0.8,
                          # fontfamily=3,
                          margin = 0.05)
pdf(paste0(outdir,"/Venn_all_URs.pdf"))
grid.draw(venn.plot)
dev.off()

# find shared DS ####
ligand_target_merged <- read.csv(file = paste0(wd,"/Data/ligand_target_merged_withFC.csv"),header = T,sep = ",") 

venn_list = list(Ovary=ligand_target_merged[ligand_target_merged$cancer=="Ovary","target"]  %>% unique,
                 Breast=ligand_target_merged[ligand_target_merged$cancer=="Breast","target"]  %>% unique,
                 Colon=ligand_target_merged[ligand_target_merged$cancer=="Colon","target"]  %>% unique,
                 Liver=ligand_target_merged[ligand_target_merged$cancer=="Liver","target"]  %>% unique,
                 Lung=ligand_target_merged[ligand_target_merged$cancer=="Lung","target"]  %>% unique)

shared_targets <- intersect(unlist(venn_list[1]), unlist(venn_list[2])) %>% 
  intersect(.,unlist(venn_list[3])) %>% 
  intersect(.,unlist(venn_list[4])) %>% 
  intersect(.,unlist(venn_list[5]))
length(shared_targets)                     
write.table(shared_targets, file=paste0(outdir,"/shared_DS.txt"), sep="\t", col.names = F, row.names = F)

venn.plot <- venn.diagram(x=venn_list, 
                          col="black",fill=cancer_color , 
                          filename = NULL,
                          cex=1,cat.cex = 1,cat.fontface = "bold",alpha=0.8,
                          # fontfamily=3,
                          margin = 0.05)
pdf(paste0(outdir,"/Venn_all_DSs.pdf"))
grid.draw(venn.plot)
dev.off()

# Identify shMCTM ####
DEG_merge_all = read.csv(paste0(wd,"/Data/DEGs_of_all_celltype_cancer_FC.csv"),sep = ';',header = T,row.names = 1)
DEG_merge_all = select(DEG_merge_all,-p_val,-pct.1,	-pct.2,	-p_val_adj)

colnames(ligand_target_merged)[7] = 'avg_log2FC.target'
ligand_target_merged_FC = merge(ligand_target_merged,DEG_merge_all,by.x=c('ligand','cancer'),by.y=c('gene','cancer'),all.x=T)
head(ligand_target_merged_FC)
ligand_target_merged_FC = ligand_target_merged_FC %>% select(-target_merge,)
colnames(ligand_target_merged_FC)[7] = 'avg_log2FC.ligand'
colnames(ligand_target_merged_FC)[8] = 'sender'
dim(ligand_target_merged_FC)
ligand_target_merged_FC = ligand_target_merged_FC[is.na(ligand_target_merged_FC$avg_log2FC.ligand )==F,]

shared_ligands = read.csv(paste0(outdir,"/shared_UR.txt"), sep="\t",header=F) %>% .[,'V1'] 
shared_targets = read.csv(paste0(outdir,"/shared_DS.txt"), sep="\t",header=F) %>% .[,'V1'] 

##shUR - Shared URs that has the same expression change in > 4 cancers in the same cell type ####
ligand_target_merged_FC.shared = ligand_target_merged_FC[ligand_target_merged_FC$ligand %in% shared_ligands & ligand_target_merged_FC$target %in% shared_targets,]

ligand_target_merged_FC.shared.UR = ligand_target_merged_FC.shared[,c('cancer','sender', 'ligand','avg_log2FC.ligand')] %>% unique
ligand_target_merged_FC.shared.UR = ligand_target_merged_FC.shared.UR[is.na(ligand_target_merged_FC.shared.UR$avg_log2FC.ligand) == F,]
ligand_target_merged_FC.shared.UR = ligand_target_merged_FC.shared.UR %>% mutate(direction = ifelse(avg_log2FC.ligand>0,1,
                                                                                                    ifelse(avg_log2FC.ligand<0,-1,NA)) )
shUR = data.frame()
for (i in unique(ligand_target_merged_FC.shared.UR$sender)){
  # i='Fibroblast'
  a = ligand_target_merged_FC.shared.UR[ligand_target_merged_FC.shared.UR$sender==i,] %>% select(cancer ,ligand, direction)
  a.wide = spread(a,key=c(ligand), value=direction) %>% as.data.frame()
  rownames(a.wide) = a.wide$cancer
  a.wide = a.wide %>% select(-cancer)
  a.wide = a.wide[rownames(a.wide) != 'NA',]
  a.wide = t(a.wide) %>% as.data.frame()
  a.wide$sum = apply(a.wide,1,function(x){return(sum(x,na.rm=T))})
  a.wide = a.wide[abs(a.wide$sum) >= 4,]
  
  if(dim(a.wide)[1]==0){next}
  
  a.wide$gene = rownames(a.wide)
  a.wide$cell = i
  a.wide$type = 'UR'
  shUR = rbind(shUR,a.wide)
}
head(shUR)
unique(shUR$gene)

a = unique(shUR$gene)
b = read.csv('/Users/yelin.zhao/Library/CloudStorage/OneDrive-Personal/Projects/cancer/datafile/outputs_5cancers/new_outputs/MCDM/scUR_celltype.csv',sep=',')
b = b$URs_same_fc %>% unique
intersect(a,b)
c = a[a %in% b == F]

##shDS - Shared DSs that has the same expression change in > 4 cancers in the same cell type ####
ligand_target_merged_FC.shared.DS = ligand_target_merged_FC.shared[,c('cancer','receiver', 'target','avg_log2FC.target')] %>% unique
ligand_target_merged_FC.shared.DS = ligand_target_merged_FC.shared.DS[is.na(ligand_target_merged_FC.shared.DS$avg_log2FC.target) == F,]
ligand_target_merged_FC.shared.DS = ligand_target_merged_FC.shared.DS %>% mutate(direction = ifelse(avg_log2FC.target>0,1,
                                                                                                    ifelse(avg_log2FC.target<0,-1,NA)) )
shDS = data.frame()
for (i in unique(ligand_target_merged_FC.shared.DS$receiver)){
  # i='Fibroblast'
  a = ligand_target_merged_FC.shared.DS[ligand_target_merged_FC.shared.DS$receiver==i,] %>% select(cancer ,target, direction)
  a.wide = spread(a,key=c(target), value=direction) %>% as.data.frame()
  rownames(a.wide) = a.wide$cancer
  a.wide = a.wide %>% select(-cancer)
  a.wide = a.wide[rownames(a.wide) != 'NA',]
  a.wide = t(a.wide) %>% as.data.frame()
  
  if(i=="Pericyte"){a.wide$Lung = NA
  a.wide = a.wide[,c('Breast','Colon','Liver','Lung','Ovary')]}
  
  a.wide$sum = apply(a.wide,1,function(x){return(sum(x,na.rm=T))})
  a.wide = a.wide[abs(a.wide$sum) >= 4,]
  
  if(dim(a.wide)[1]==0){next}
  
  a.wide$gene = rownames(a.wide)
  a.wide$cell = i
  a.wide$type = 'DS'
  shDS = rbind(shDS,a.wide)
}
head(shDS)
unique(shDS$gene) %>% length

shMCTM = rbind(shUR,shDS)
shMCTM$cell_gene = paste0(shMCTM$cell,'_',shMCTM$gene)
unique(shMCTM$gene) %>% length

write.table(shMCTM, file=paste0(outdir,"/shMCTM_genes.txt"), sep="\t", col.names = T, row.names = F)

## shMCTM interaction ####
ligand_target_merged_FC.shared = ligand_target_merged_FC %>% select(-cancer,-avg_log2FC.ligand ,-avg_log2FC.target) %>% unique
ligand_target_merged_FC.shared$sender_ligand = paste0(ligand_target_merged_FC.shared$sender,'_',ligand_target_merged_FC.shared$ligand)
ligand_target_merged_FC.shared$receiver_target = paste0(ligand_target_merged_FC.shared$receiver,'_',ligand_target_merged_FC.shared$target)

shMCTM_interaction = ligand_target_merged_FC.shared[ligand_target_merged_FC.shared$sender_ligand %in% shMCTM[shMCTM$type=='UR','cell_gene'],]
shMCTM_interaction = shMCTM_interaction[shMCTM_interaction$receiver_target %in% shMCTM[shMCTM$type=='DS','cell_gene'],]

head(shMCTM_interaction)
write.table(shMCTM_interaction, file=paste0(outdir,"/shMCTM_interaction.txt"), sep="\t", col.names = T, row.names = F)

#########################
# Identify top shURs ####
library(reshape2)
# df1 = read.csv(paste0(outdir,'/shMCTM_comb2_shMCTM_secondtry.csv'),row.names = 1) %>% as.data.frame %>% unique
shMCTM_interaction = read.csv2(paste0(outdir,"/shMCTM_interaction.txt"), sep="\t")

df = as.data.frame(shMCTM_interaction[,c('ligand',    'receiver', 'target',   'sender')]) %>% unique
df.wide = dcast(df, receiver ~ ligand, value.var = "target", fun.aggregate = length) 
rownames(df.wide) = df.wide$receiver
df.wide = df.wide %>% select(-receiver)
df.wide.t = t(df.wide)  
df.wide.t.scaled =scale(df.wide.t) #scaling is necessary to get right distance metrics

# decide how many clusters we will cut####
#Deciding based on multiple metrics
library(NbClust)
res <- NbClust(df.wide.t, diss=NULL, distance = "euclidean",  method = "complete", index = "all")

# check elbow
library("factoextra")
P = fviz_nbclust(df.wide.t, hcut, hc_method="complete", method = "wss") +
  geom_vline(xintercept = 2, linetype = 2)+
  labs(subtitle = "Elbow method")
pdf(file=paste0(outdir,"/Elbow_clusterof_UR_selection_scaled.pdf"), width = 5, height = 5)
print(p)
dev.off()

data = df.wide.t
# make heatmap ####
p = Heatmap(data,
            # split
            row_split = 2,
            # column_split = 2,  #split by dendrograms
            # row_km = 2,
            # column_km = 2 , #split rows and colums
            # cluster_column_slices = F,
            # split = data.frame(cyl = mtcars$cyl, am = mtcars$am), 
            # column_split = anno_col$Group,
            #设置title,legend,font
            rect_gp = gpar(col = "white"),
            name = "No.interactions", #title of legend
            column_title = "Downstream cell types", 
            column_title_side = "bottom",  
            column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
            column_names_gp = gpar(fontsize = 9),
            row_title = "shURs",
            row_title_side ="right", 
            row_title_gp = gpar(fontsize = 12) ,
            row_names_gp = gpar(fontsize = 9), # Text size for row names
            cluster_columns = T,
            cluster_rows = T,
            border = T,
            col = colorRamp2(c(min(data),max(data)), c("white","red")),
)
# pdf(file=paste0(outdir,"heatmap_UR_selection_nonscaled.pdf"), width = 5, height = 8)
pdf(file=paste0(outdir,"/heatmap_UR_selection_scaled.pdf"), width = 5, height = 8)
print(p)
dev.off()

draw(p)
ro = row_order(p)[[2]]
draw(p)
row_order(p)
rownames(data)[ro]
top_shURs = rownames(data)[ro]
write.csv(top_shURs,file = paste0(outdir,'/top_shURs.csv'))




