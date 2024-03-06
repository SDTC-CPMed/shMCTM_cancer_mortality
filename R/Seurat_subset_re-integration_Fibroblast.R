#analysis based on each subset/Fibroblast
rm(list = ls())
library(tidyverse)
library(nichenetr)
library(Seurat) # Seurat V4
library(MAST)
library(patchwork)
library(ggplot2)
library(SingleR)
library(SummarizedExperiment)
library(reshape2)

######################################################################################## 
#1 pilot - look all cell types####
seed = 42
nPCs = 30
res= 0.8
reduction= "umap"
genes_to_check = c('PTPRC','CD45', ## immune cell marker ('PTPRC'='CD45')
                   'PTPRC','CD3D', 'CD3E','CD4', 'CD8A','FOXP3','KLRD1', ## Tcells
                   'CD8A','GZMK','CD8B',  ## CD8+T cells
                   # 'PDCD1','LAG3','HAVCR2','TOX',# CD8+T exhausted (only high with these markers) and Proliferating
                   # 'MKI67','TOP2A','STMN2', #CD8A+ T proliferating -except above, also highly express these genes
                   # 'JUNB', 'FOS', 'ATF3', 'HSPA1A','DNAJB1', 'DNAJB1',#CD8A+T exhausted IEG
                   'GZMB', # shared by CTL and NK 
                   'CD4','IL7R', ##CD4+T
                   'FOXP3','CD4','IL2RA','IKZF2','TNFRSF4',  ## Treg
                   'CD40LG','IL7R', ## T memory
                   'CD40LG','IL7R','STAT4','CD3G', ## T helper
                   'CD16','CD16A','CD16B','CCR5','CCR6','CD27','TNFRSF7','TNFSF7', #T gamma delta
                   'GNLY','KLRD1','KLRC1','KLRF1','GZMB','NKG7','NCAM1','HAVR2C','NCR1','FGFBP2','FCGR3A','PRF1',  ##NK(CD56=NCAM1) (but no cd3d,cd4,cd8a)
                   'ZNF683','CD8A','CD8B','CD3','NCAM1','CD56', ##NKT  
                   'CD19', 'CD79A','CD79B','MS4A1' , # B cells
                   'IGHG1', 'MZB1', 'SDC1','IGHG4','CD38',  # plasma 
                   'CD68', 'LYZ', 'AIF1', #myeloid
                   'VCAN','FCN1','S100A12', ## monocyte
                   'CD68', 'CD163', 'CD14',  'CD86', 'CCL22','S100A4','CD207','CCL17', ## DC(belong to monocyte)
                   'CD40','CD80', 'HLA-DOB','DOB', # DC-activated
                   'CD1C','CLEC9A', ## mDCs
                   'LILRA4','IL3RA','TCF4','TCL1A','CLEC4C', ## pDCs
                   'CD68','CD163', 'LAMP2','LGALS3', 'MRC1','MSR1' ,'S100A8','CD14','CD11B','APOE','C1QA','C1QB','ITGMX','CD11C','ITGAM','CD11B', ## Macrophage (belong to monocyte)
                   'FCGR3B', ## Neutrophil
                   'CD33','KIT','VIM','MS4A2','TPSAB1', 'CPA3',  ##Mast cells
                   'EPCAM', 'KRT19','KRT7','KRT8','KRT18', 'PROM1',  ## epi or tumor
                   'FGF7', 'MME','COL1A1','ACTA2','PECAM1', 'VWF' ,'PROX1','PDGFRA', ## Fibroblasts,Endothelial
                   'MME','CD10','ACTA2', 'COL1A1', 'FN1', 'BGN', 'DCN',     ##stromal_fibroblasts('MME'='CD10')
                   'PECAM1',"CD31", 'ENG','VWF', 'CD36',  ##stromal_endo('PECAM1'="CD31")
                   'MCAM','RGS5' #Pericytes
                   # 'ASGR1', # hepatocyte
                   # 'FXYD2',  # cholangiocyte
                   # 'CLDN18', 'SFTPA1', 'SFTPA2', 'SFTPC', #alveolar cell 
                   # 'S100B', 'PLP1' #enteric glia- enteric nervous system
)

##. Visualization for major cell type####
for (i in c("Myeloids","Lymphocytes","Stromal","Epithelial")){
  load(paste0(outdir_subset,"/all_integrated_",i,".RData"))
  
  DefaultAssay(all.integrated) <- "integrated"

  p1 <- DimPlot(all.integrated, reduction = reduction, group.by = "group")
  p2 <- DimPlot(all.integrated, reduction = reduction, group.by = "Celltype_final",label = TRUE)
  # p3 <- DimPlot(all.integrated, reduction = reduction, split.by = "group") #To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
  p3 <- DimPlot(all.integrated, reduction = reduction, group.by = "cancer_name")
  p4 <- DimPlot(all.integrated, reduction = reduction, split.by = "group" )
  p =  (p1 + p2 + p3)
  ggsave(plot=p,  path = outdir_subset_PCA, dpi = 100, width = 16, height = 6, units = "in" ,filename=paste0("/DimPlot_celltypefinal_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))

  FeaturePlot(all.integrated, features = c("MDK","SPP1","PLAU",'HOXAS'),split.by = "group", max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  p <- FeaturePlot(all.integrated, features = c("MDK","SPP1","PLAU",'IL6','LIF','CCL5','CSF3','ACTA2','FAP','CTGF'),split.by = "group", max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  ggsave(plot=p,  path = outdir_subset_PCA, device = "pdf", width = 12, height = 45 ,filename=paste0("Feature_Plot","_",i,"_bygroup_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  p <- FeaturePlot(all.integrated, features = c("MDK","SPP1","PLAU",'IL6','LIF','CCL5','CSF3','ACTA2','FAP','CTGF'),split.by = "cancer_name", max.cutoff = 3, label = TRUE)
  ggsave(plot=p,  path = outdir_subset_PCA, device = "pdf",width = 25, height = 45 ,filename=paste0("Feature_Plot","_",i,"_bycancer_",reduction,"_seed",seed,"_res",res,".pdf"))
  

  table(all.integrated$seurat_clusters)
  table(all.integrated$seurat_clusters,all.integrated$group)

  DefaultAssay(all.integrated) <- "RNA"
  # p <- FeaturePlot(all.integrated, features = unique(genes_to_check), max.cutoff = 3, label = TRUE,
  #                  cols = c("grey", "red"))
  # ggsave(plot=p,  path = outdir_subset_PCA, device = "png",dpi = 150, width = 18, height = 45, units = "in" ,filename=paste0("Feature_Plot","_",i,"_seed",seed,"_res",res,".png"))

  p <- DotPlot(all.integrated, features = unique(genes_to_check),
               assay='RNA' , ##can use this to check orignial or normalized expression
               group.by = 'Celltype_final' ##this can change to other group method like 'seurat_clusters'
  )  + scale_color_viridis_c()  + coord_flip()
  ggsave(plot=p,  path = outdir_subset_PCA, dpi = 300, width = 30, height = 45, units = "cm" ,filename=paste0("check_markers_celltypefinal_",i,"_seed",seed,"_res",res,".pdf"))

  print(paste0(i, " done!"))
}

######################################################################################## 
#2. re-sub merging for Fibroblast####
##2.1. subset celltype from each cancers ####
Fibroblast_sublist = list()

outdir = "/"
file_list= list.files(path=outdir, pattern="all_integrated_Celltype_final_")
cancer_name= sapply(strsplit(file_list,".RData"),'[[',1)
cancer_name= sapply(strsplit(cancer_name,"_"),'[[',5)

for (i in 1:length(file_list)) {
  load(paste0(outdir,"/",file_list[i]))
  print(i)
  DefaultAssay(all.integrated) <- "RNA"
  seuratObj = all.integrated
  seuratObj@meta.data %>% head()
  
  seuratObj@meta.data$Celltype_final %>% table()
  seuratObj@meta.data$cancer_name = cancer_name[i]
  
  Fibroblast_sub <- subset(x = seuratObj, 
                        subset = Celltype_final == "Fibroblast")
  
  Fibroblast_sublist = c(Fibroblast_sublist,Fibroblast_sub)
  
  print(paste0(cancer_name[i], " done!"))
  rm(all.integrated,seuratObj)
}

save(Fibroblast_sublist,file = paste0(outdir_subset_final,"/Fibroblast_sublist.RData"))

Fibroblast_sublist[[1]]$Celltype_final %>% table
Fibroblast_sublist[[1]]@meta.data %>% head()


##2.2. merge the same celltype across cancers - seurat integration####
library(SeuratDisk)
library(MySeuratWrappers)
library(SeuratWrappers)
###2.2.1 select integration features####
Fibroblast_sublist <- lapply(X = Fibroblast_sublist, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Fibroblast_sublist)
save(features,file = paste0(outdir_subset_final,"/features_Fibroblast.RData"))

###2.2.2 integration  #### 
load(paste0(outdir_subset_final,"/features_Fibroblast.RData"))
anchors <- FindIntegrationAnchors(object.list = Fibroblast_sublist, anchor.features = features)
all.integrated <- IntegrateData(anchorset = anchors)
save(all.integrated,file = paste0(outdir_subset_final,"/all_integrated_Fibroblast_beforePCA.RData"))
rm(features)

###2.2.3 Clustering after integration #### 
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
for (i in c("Fibroblast")){
  i="Fibroblast"
  load(paste0(outdir_subset_final,"/all_integrated_",i,"_beforePCA.RData"))
  all.integrated <- subset(x = all.integrated, 
                           subset = (cancer_name != "nopericyteLung" & cancer_name !="withpericyteLung"))
  all.integrated@meta.data <- all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  DefaultAssay(all.integrated) <- "integrated"
  seed = 42
  nPCs = 10
  res= 0.1
  reduction= "umap"
  
  all.integrated <- ScaleData(all.integrated, verbose = T)
  all.integrated <- RunPCA(all.integrated,seed.use = seed, verbose = T)
  p <- ElbowPlot(all.integrated,ndims = 30)
  ggsave(plot=p,  path = outdir_subset_final, dpi = 300, width = 25, height = 20, units = "cm" ,filename=paste0(i,"_PCA-ElbowPlot.pdf"))
  all.integrated <- RunUMAP(all.integrated, seed.use = seed, reduction = "pca", dims = seq(nPCs))
  all.integrated <- RunTSNE(all.integrated, seed.use = seed, reduction = "pca", dims=seq(nPCs))
  all.integrated <- FindNeighbors(all.integrated, reduction = "pca", dims = seq(nPCs))
  all.integrated <- FindClusters(all.integrated, random.seed = seed, resolution = res)
  table(all.integrated$seurat_clusters)
  
  # DefaultAssay(all.integrated) <- "RNA"
  # all.integrated <- NormalizeData(all.integrated)
  # all.integrated <- ScaleData(all.integrated, verbose = T)
  all.integrated@assays$RNA@counts[1:70,1:2]
  all.integrated@assays$RNA@data[1:70,1:2]
  all.integrated@assays$RNA@scale.data[1:70,1:2]
  
  all.integrated@assays$integrated@counts[1:70,1:2]
  all.integrated@assays$integrated@data[1:70,1:2]
  all.integrated@assays$integrated@scale.data[1:70,1:2]
  
  save(all.integrated, file = paste0(outdir_subset_final,"/all_integrated_",i,".RData"))
}

#visualization
for (i in c("Fibroblast")){
  i="Fibroblast"
  load(paste0(outdir_subset_final,"/all_integrated_",i,".RData"))
  all.integrated@meta.data <- all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  DefaultAssay(all.integrated) <- "integrated"
  seed = 42
  nPCs = 10
  res= 0.1
  reduction= "umap"
  
  #DimPlot
  p1 <- DimPlot(all.integrated, reduction = reduction, group.by = "group")
  p2 <- DimPlot(all.integrated, reduction = reduction, label = TRUE)
  p3 <- DimPlot(all.integrated, reduction = reduction, group.by = "cancer_name")
  p =  (p1 + p2 + p3)
  ggsave(plot=p,  path = outdir_subset_final, width = 9, height = 3,filename=paste0("/DimPlot_celltypefinal_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  p <- DimPlot(all.integrated, reduction = reduction, split.by = "group",label=T) #To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
  ggsave(plot=p,  path = outdir_subset_final, width = 6, height = 3,filename=paste0("/DimPlot_celltypefinal_bygroup_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  p <- DimPlot(all.integrated, reduction = reduction, split.by = "cancer_name",label=T) #To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
  ggsave(plot=p,  path = outdir_subset_final, width = 9, height = 2,filename=paste0("/DimPlot_celltypefinal_bycancer_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  #Dotplot
  DefaultAssay(all.integrated) <- "RNA"
  markercluster0 = c("COL1A1" ,"FN1" ,"MMP11","CTHRC1" ,  "COL1A2","COL3A1,SPARC","COL5A2","POSTN")
  top8UR = c("FN1","COL1A1","PLAU" ,"CLEC11A","MDK")
  DS_pw = c('NFKB1', 'VIM', 'TGFB1', 'COL1A1', 'COL1A2', 'COL3A1', 'MMP1', 'MMP10', 'MMP9')
  UR_pw = c('CYR61','RARRES2' )
  CAF_markers = c('ACTA2','IL6','IL10','LIF','FBLN1','PDGFRA','POSTN','KRT19','CD74')
  genes_to_check = c(CAF_markers,markercluster0,top8UR) %>% unique

  all.integrated$new_id = paste0(all.integrated$seurat_clusters,'_',all.integrated$group)
  p=DotPlot(object = all.integrated, features = CAF_markers, 
            assay="RNA",group.by='new_id') + guides(color = guide_colorbar(title = 'Scaled Average Expression'))+
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)+
    theme(axis.text.x=element_text(angle=90,hjust=1))
  ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 7, height = 7 ,filename=paste0("DotPlot_CAF_markers","_",i,"_bygroup_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  p=DotPlot(object = all.integrated, features = CAF_markers, 
            assay="RNA",cols = c("white", "red"))+RotatedAxis()  + 
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) +
    guides(color = guide_colorbar(title = 'Scaled Average Expression'))+ 
    theme(axis.text.x=element_text(angle=90,hjust=1)) 
  ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 7, height = 7 ,filename=paste0("DotPlot_CAF_markers","_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
 
  #FeaturePlot
  DefaultAssay(all.integrated) <- "RNA"
  p = FeaturePlot(all.integrated, features = c('IL6','ACTA2','FN1', 'COL1A1', 'MMP11', 'CTHRC1', 'COL1A2','POSTN','CD74', 'HLA-DRA','KRT19'),
                  # split.by = "group", 
                  ncol = 6,
                  max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 20, height = 8 ,filename=paste0("Feature_Plot","_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))

  p = FeaturePlot(all.integrated, features = c('ACTA2','FN1','IL6', 'HSPA1A', 'KRT7','CD74'),
                  # split.by = "group", 
                  ncol = 3,
                  max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 16, height = 8 ,filename=paste0("Feature_Plot","_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  p = FeaturePlot(all.integrated, features = c('IL6','ACTA2','FN1', 'COL1A1', 'MMP11', 'CTHRC1', 'COL1A2','POSTN','CD74', 'HLA-DRA','KRT19'),
                  split.by = "group",
                  ncol = 6,
                  max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 8, height = 22 ,filename=paste0("Feature_Plot","_",i,"_bygroup_",reduction,"_seed",seed,"_res",res,".pdf"))
  
}

###2.2.4 cell count####
all.integrated$cancer_group = paste0(all.integrated$cancer_name,'_',all.integrated$group)
count = table(all.integrated$seurat_clusters,all.integrated$cancer_group) %>% as.matrix() %>% as.data.frame.matrix()
write.csv(count,file = paste0(outdir_subset_final,"/Cluster_cellcount.csv"),sep=",")
percent = apply(count,2,function(x){x/sum(x)})
percent.long = melt(percent, value.name = "Percentage")
colnames(percent.long) = c("Cluster","Cancer_group","Percentage")
percent.long$Cluster = factor(percent.long$Cluster )
percent.long$Cancer = lapply(as.character(percent.long$Cancer_group), function(x){strsplit(x,"_") %>% unlist %>% .[1]}) %>% unlist
percent.long$Group = lapply(as.character(percent.long$Cancer_group), function(x){strsplit(x,"_") %>% unlist %>% .[2]}) %>% unlist

head(percent.long)

p = ggplot(percent.long, aes(x = Group, y = Percentage,fill= Cluster)) + 
  geom_bar(stat = "identity", color = "black") + scale_color_viridis_c() +
  facet_grid(~Cancer,scales='fixed')
pdf(file=paste0(outdir_subset_final,"/Cluster_percentage.pdf"), width = 5, height = 5)
print(p)
dev.off()

##2.3 find markers for clusters ####
for (i in c("Fibroblast")){ #findallmarkers
  load(paste0(outdir_subset_final,"/all_integrated_",i,".RData"))
  print(i)
  DefaultAssay(all.integrated) <- "RNA"
  
  all.integrated@meta.data = all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  Idents(all.integrated) <- all.integrated$seurat_clusters
  head(Idents(all.integrated))
  
  all.markers <- FindAllMarkers(all.integrated, only.pos = TRUE, 
                                min.pct = 0.1, logfc.threshold = 0.25)
  write.table(all.markers,
              file=paste0(outdir_subset_final,"/ClusterMarkers/",i,"_total_marker_genes_",nPCs,"PC.txt"),
              sep="\t",quote = F,row.names = F)
  
  all.markers = read.csv2(paste0(outdir_subset_final,"/ClusterMarkers/",i,"_total_marker_genes_",nPCs,"PC.txt"),
                          sep="\t")
  
  all.markers.ranked = all.markers[order(all.markers$cluster,-as.numeric(all.markers$avg_log2FC)),]
  
  # top genes per cluster
  marker.sig <- all.markers %>% 
    mutate(Ratio = round(as.numeric(pct.1)/as.numeric(pct.2),3)) %>%
    filter(as.numeric(all.markers$p_val_adj) <= 0.05)  # 
  
  for(cluster_id in unique(marker.sig$cluster)){
    cl4.genes <- marker.sig %>%
      filter(cluster == cluster_id) %>%
      arrange(desc(avg_log2FC))
    cl4.genes <- cl4.genes[1:min(nrow(cl4.genes),10),"gene"]
  }
  
  top <- marker.sig %>% group_by(cluster) %>% 
    top_n(n = 30, wt = avg_log2FC)
  
  #top-marker- dotplot
  pdf(paste0(outdir_subset_final,"/ClusterMarkers/MarkerGene-DotPlot_all_cluster_",i,"_","res",res,"_PC",nPCs,"PC.pdf"),width = 20,height = 5)
  DotPlot(all.integrated, features = unique(top$gene))+
    # scale_color_viridis_c()+
    RotatedAxis() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
  dev.off()
  
  #top-marker- heatmap
  DefaultAssay(all.integrated)  <- 'integrated'
  pdf(paste0(outdir_subset_final,"/ClusterMarkers/MarkerGene-Heatmap_all_cluster_",i,"_","res",res,"_PC",nPCs,"PC.pdf"),width= 10, height= 10 )
  DoHeatmap(all.integrated, features = top$gene,size = 2) +
    theme(legend.position = "none", 
          axis.text.y = element_text(size = 7))
  dev.off()
}

###2.3.0 dotPlot all scMCTM genes ####
load(paste0(outdir_subset_final,"/all_integrated_Fibroblast.RData"))
all.integrated@meta.data <- all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
# DefaultAssay(all.integrated) <- "integrated"
DefaultAssay(all.integrated) <- "RNA"
table(all.integrated$seurat_clusters)
  
# shMCTM genes that highly expressed in fibroblast
allscUR.fib = read.csv(paste0(wd, "/cancer/datafile/outputs_5cancers/new_outputs/MCDM/CT_high_shUR_shDS/Findallmarkers_celltype/marker.high.shUR.csv"))
allscUR.fib = allscUR.fib[allscUR.fib$type =="shUR" & allscUR.fib$celltype == "Fibroblast","gene"]
allscDS.fib = read.csv(paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/MCDM/CT_high_shUR_shDS/Findallmarkers_celltype/marker.high.shDS.csv"))
allscDS.fib = allscDS.fib[allscDS.fib$type =="shDS" & allscDS.fib$celltype == "Fibroblast","gene"]
  
#order by cluster marker
all.markers = read.csv2(paste0(outdir_subset_final,"/ClusterMarkers/Fibroblast_total_marker_genes_",nPCs,"PC.txt"),
                          sep="\t")
all.markers.ranked = all.markers[order(all.markers$cluster,-as.numeric(all.markers$avg_log2FC)),]
  
allscUR.ranked = all.markers.ranked$gene[all.markers.ranked$gene %in% allscUR.fib] %>% unique
allscDS.ranked = all.markers.ranked$gene[all.markers.ranked$gene %in% allscDS.fib] %>% unique

#plot - order based on cluster marker order
p= DotPlot(object = all.integrated, features = allscUR.ranked, assay="RNA")  +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
  data_for_plot = p$data # get data of average expression and percentage expression
ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 8, height = 7 ,filename=paste0("DotPlot","_fibroblast_allscUR.ordered","_seed",seed,"_res",res,".pdf"))
  
all.integrated$new_id = paste0(all.integrated$seurat_clusters,'_',all.integrated$group)
p= DotPlot(object = all.integrated, features = allscUR.ranked, assay="RNA",group.by = "new_id")  +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
  data_for_plot = p$data # get data of average expression and percentage expression
ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 8, height = 7 ,filename=paste0("DotPlot","_fibroblast_allscUR_bygroup.ordered","_seed",seed,"_res",res,".pdf"))
  
p= DotPlot(object = all.integrated, features = allscDS.ranked, assay="RNA")  +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
data_for_plot = p$data # get data of average expression and percentage expression
ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 16, height = 7 ,filename=paste0("DotPlot","_fibroblast_allscDS.ordered","_seed",seed,"_res",res,".pdf"))
  
p= DotPlot(object = all.integrated, features = allscDS.ranked, assay="RNA",group.by = "new_id")  +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
data_for_plot = p$data # get data of average expression and percentage expression
ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 16, height = 7 ,filename=paste0("DotPlot","_fibroblast_allscDS_bygroup.ordered","_seed",seed,"_res",res,".pdf"))
  

### 2.3.1. Check whether top 10 cluster0 marker genes were DEGs against control ####
#order by cluster marker
all.markers = read.csv2(paste0(outdir_subset_final,"/ClusterMarkers/Fibroblast_total_marker_genes_",nPCs,"PC.txt"),
                        sep="\t")
all.markers.ranked = all.markers[order(all.markers$cluster,-as.numeric(all.markers$avg_log2FC)),]
all.markers.ranked.C0 = all.markers.ranked[all.markers.ranked$cluster ==0,"gene"][1:10]

#DEGs tumor vs. control - all Fibroblast
DEGs_of_all_celltype_cancer = read.csv2(paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/DEGs/merge_DEGs_of_all_celltype_cancer.csv"),sep=",",row.names = 1)
DEGs_of_all_celltype_cancer.fib = DEGs_of_all_celltype_cancer[DEGs_of_all_celltype_cancer$cell == "Fibroblast" & DEGs_of_all_celltype_cancer$X %in% all.markers.ranked.C0,]
table(DEGs_of_all_celltype_cancer.fib$X,DEGs_of_all_celltype_cancer.fib$cancer)

write.table(DEGs_of_all_celltype_cancer.fib,
            file=paste0(outdir_subset_final,"/ClusterMarkers/top10C0_DEG_TvsN_Fib.txt"),
            sep="\t",quote = F,row.names = F)

#DEGs tumor vs. control - in C0
load(paste0(outdir_subset_final,"/all_integrated_Fibroblast.RData"))
DefaultAssay(all.integrated) <- "RNA"

all.integrated@meta.data = all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
Idents(all.integrated) <- all.integrated$seurat_clusters
head(Idents(all.integrated))

DEG_method <- "MAST"  #wilcox or MAST和DESeq2
DEG = data.frame()
for (i in 0:(length(unique(all.integrated$seurat_clusters))-1)) {
    print(i)
  # for (j in 1:5){
  #   print(unique(temp$cancer_name)[j])
    Idents(all.integrated) <- "seurat_clusters"
    temp <- subset(all.integrated, idents = i)
    # temp <- subset(temp, cancer_name == unique(temp$cancer_name)[j])
    Idents(temp) <- "group"
    
    if(dim(table(temp$Celltype_final,temp$group))[2] ==1){next} #skip if cell type only exists in one group
    
    iDEGs <- FindMarkers(temp, ident.1 = "T", ident.2 = "N", verbose = FALSE,
                         min.pct = 0.1, #only select genes that expressed in more than 10% of cells 
                         test.use = DEG_method #"wilcox"(default), can change to wilcox or MAST和DESeq2
    )
    head(iDEGs, n = 15)
    iDEGs <- iDEGs %>% filter(p_val_adj <= 0.05)
    iDEGs$cluster = i
    DEG = rbind(DEG,iDEGs)
    # }
}
write.csv(DEG,paste0(outdir_subset_final,"/ClusterMarkers/DEG_TvsN_allclusters","_",DEG_method,".csv"),row.names = T, col.names = T, sep = ',',quote = F)

DEGs = read.csv2(paste0(outdir_subset_final,"/ClusterMarkers/DEG_TvsN_allclusters","_",DEG_method,".csv"),sep=",")
DEGs_of_all_celltype_cancer.fib = DEGs[DEGs$X %in% all.markers.ranked.C0,]

write.table(DEGs_of_all_celltype_cancer.fib,
            file=paste0(outdir_subset_final,"/ClusterMarkers/top10C0_DEG_TvsN_C0.txt"),
            sep="\t",quote = F,row.names = F)

### 2.3.2. GO,KEGG for cluster0 marker ####
  library(org.Hs.eg.db)
  library(org.Hm.eg.db)  
  library(clusterProfiler)  
  library(dplyr)  
  library(ggplot2)  
  library(DOSE)
  library(topGO)
  library(pathview)
  #Write function to create all GO, KEGG output####
  GO_Kegg_Enrichment_Vis <- function(x,datasource,nshow=10,fromtype="SYMBOL",totype=c("ENSEMBL","ENTREZID"),orgdb="org.Hs.eg.db", out.dir){
    KEGG_results <- matrix(ncol=9)
    colnames(KEGG_results) <- c('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID','Count')
    ##translate####
    test = bitr(x,  
                fromType=fromtype,  
                toType=totype,   
                OrgDb=orgdb) 
    head(test)
    
    ##KEGG enrichment####
    kk <- enrichKEGG(gene = test$ENTREZID,
                     organism = 'hsa', #human ='hsa', mouse= 'mmu' #https://www.genome.jp/kegg/catalog/org_list.html 
                     pvalueCutoff = 1)
    head(kk,2)
  
    #show results
    write.csv(summary(kk),paste0(out.dir,"/KEGG-enrich_",datasource,".csv"),row.names =FALSE)
    write.csv(kk,file = paste0(out.dir,"/KEGG-enrich_all_",datasource,".csv"),row.names =T)
    KEGG_results <- kk
    barplot
    p = dotplot(kk,orderBy = "x",title="Enrichment KEGG_dot",showCategory=nshow)
    ggsave(plot=p,path = paste0(out.dir), filename = paste0("KEGG-enrich_dotplot_",datasource,".pdf"),dpi = 300, width = 15, height = 5, units = "cm")
    p = barplot(kk,orderBy = "x",title="Enrichment KEGG_dot",showCategory=nshow)
    ggsave(plot=p,path = paste0(out.dir), filename = paste0("KEGG-enrich_barplot_",datasource,".pdf"),dpi = 300, width = 15, height = 5, units = "cm")

    ##GO enrichment ####
    ggo <- groupGO(gene = test$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC",level = 3,readable = TRUE)
    ego_ALL <- enrichGO(gene = test$ENTREZID,   
                        keyType = "ENTREZID",  
                        OrgDb=org.Hs.eg.db,  
                        ont = "ALL",    
                        pvalueCutoff = 0.01,  
                        pAdjustMethod = "fdr",  
                        minGSSize = 1,   
                        maxGSSize = 500,   
                        qvalueCutoff = 0.05,   
                        readable = TRUE)   
    dim(ego_ALL)
    head(ego_ALL@result)
    ego_ALL@result$ONTOLOGY %>% table
    ego_ALL <- setReadable(ego_ALL, OrgDb = org.Hs.eg.db)
    write.csv(summary(ego_ALL),paste0(out.dir,"/GO-all_enrich_",datasource,".csv"),row.names =FALSE)

    p=barplot(ego_ALL,showCategory=nshow,
              font.size=5,
              split="ONTOLOGY",
              color ="pvalue",title="EnrichmentGO_ALL")+
      facet_grid(ONTOLOGY~.,scale="free")+
      scale_y_discrete(labels=function(x) str_wrap(x, width=60))#条状图，按p从小到大排，绘制前10个Term
    ggsave(plot=p,path = paste0(out.dir), filename = paste0("GO-all_barplot_",datasource,".pdf"),dpi = 300, width = 30, height = 20, units = "cm")
  }
  
  out.dir = paste0(outdir_subset_final,'/GOKEGG')
  for (cluster in unique(all.markers.ranked$cluster)){
    if (cluster == 3){
      gene_of_interest = all.markers.ranked[(all.markers.ranked$cluster == cluster & all.markers.ranked$avg_log2FC > 0.5),"gene"]
    }else
    {gene_of_interest = all.markers.ranked[(all.markers.ranked$cluster == cluster & all.markers.ranked$avg_log2FC > 1.25),"gene"]}
    
    GO_Kegg_Enrichment_Vis(gene_of_interest,datasource=paste0("Fib_markers_cluster",cluster),nshow=5,out.dir = out.dir)
    print(paste0("cluster done: ", cluster))
  }
  







