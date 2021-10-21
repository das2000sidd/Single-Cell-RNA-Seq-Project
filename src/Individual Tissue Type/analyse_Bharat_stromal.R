setwd("/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells")


library(Seurat)

SC89_raw=Read10X('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/CellRanger_run_output/All_Cell_Ranger_Output/SC89_bronchus_stromal/outs/raw_feature_bc_matrix/')
SC86_raw=Read10X('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/CellRanger_run_output/All_Cell_Ranger_Output/SC86_parenchyma_stromal/outs/raw_feature_bc_matrix/')
SC108_raw=Read10X('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/CellRanger_run_output/All_Cell_Ranger_Output/SC108_lung_stromal/outs/raw_feature_bc_matrix/')
SC143_raw=Read10X('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/CellRanger_run_output/All_Cell_Ranger_Output/SC143_lung_stromal/outs/raw_feature_bc_matrix/')
SC151_raw=Read10X('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/CellRanger_run_output/All_Cell_Ranger_Output/SC151_lung_stromal/outs/raw_feature_bc_matrix/')
SC147_raw=Read10X('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/CellRanger_run_output/All_Cell_Ranger_Output/SC147_lung_stromal/outs/raw_feature_bc_matrix/')


SC86=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC86_seurat_min_3_cells_min_200_genes.rds')
SC89=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC89_seurat_min_3_cells_min_200_genes.rds')
SC108=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC108_seurat_min_3_cells_min_200_genes.rds')
SC143=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC143_seurat_min_3_cells_min_200_genes.rds')
SC151=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC151_seurat_min_3_cells_min_200_genes.rds')
SC147=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC147_seurat_min_3_cells_min_200_genes.rds')

## SC86 is Donor1_3
## SC89 is Donor1_6
## SC108 is PB1_2
## SC143 is TX_Donor_3
## SC147 is TX_Recipient_3
## SC151 is PB2_3

SC86_doublet=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Donor1_3_doublets.csv",header = F,stringsAsFactors = F)
SC89_doublet=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Donor1_6_doublets.csv",header = F,stringsAsFactors = F)
SC108_doublet=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/PB1_2_doublets.csv",header = F,stringsAsFactors = F)
SC143_doublet=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/TX_Donor_3_doublets.csv",header = F,stringsAsFactors = F)
SC147_doublet=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/TX_Recipient_3_doublets.csv",header = F,stringsAsFactors = F)
SC151_doublet=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/PB2_3_doublets.csv",header = F,stringsAsFactors = F)



library(scDblFinder)


SC86.sce <- as.SingleCellExperiment(SC86)
SC89.sce <- as.SingleCellExperiment(SC89)
SC108.sce <- as.SingleCellExperiment(SC108)
SC143.sce <- as.SingleCellExperiment(SC143)
SC147.sce <- as.SingleCellExperiment(SC147)
SC151.sce <- as.SingleCellExperiment(SC151)



SC86.sce=scDblFinder(SC86.sce,verbose = TRUE)
SC89.sce=scDblFinder(SC89.sce,verbose = TRUE)
SC108.sce=scDblFinder(SC108.sce,verbose = TRUE)
SC143.sce=scDblFinder(SC143.sce,verbose = TRUE)
SC147.sce=scDblFinder(SC147.sce,verbose = TRUE)
SC151.sce=scDblFinder(SC151.sce,verbose = TRUE)


#SC86.seurat.with.doublet=as.Seurat(SC86.sce)

SC86.doublet.info=data.frame(colnames(SC86.sce),SC86.sce$scDblFinder.class)
SC89.doublet.info=data.frame(colnames(SC89.sce),SC89.sce$scDblFinder.class)
SC108.doublet.info=data.frame(colnames(SC108.sce),SC108.sce$scDblFinder.class)
SC143.doublet.info=data.frame(colnames(SC143.sce),SC143.sce$scDblFinder.class)
SC147.doublet.info=data.frame(colnames(SC147.sce),SC147.sce$scDblFinder.class)
SC151.doublet.info=data.frame(colnames(SC151.sce),SC151.sce$scDblFinder.class)

#SC86.doublet=subset(SC86.doublet.info,SC86.doublet.info$SC86.sce.scDblFinder.class=="doublet")
#SC89.doublet=subset(SC89.doublet.info,SC89.doublet.info$SC89.sce.scDblFinder.class=="doublet")
#SC108.doublet=subset(SC108.doublet.info,SC108.doublet.info$SC108.sce.scDblFinder.class=="doublet")
#SC143.doublet=subset(SC143.doublet.info,SC143.doublet.info$SC143.sce.scDblFinder.class=="doublet")
#SC147.doublet=subset(SC147.doublet.info,SC147.doublet.info$SC147.sce.scDblFinder.class=="doublet")
#SC151.doublet=subset(SC151.doublet.info,SC151.doublet.info$SC151.sce.scDblFinder.class=="doublet")


rownames(SC86.doublet.info)=SC86.doublet.info$colnames.SC86.sce.
rownames(SC89.doublet.info)=SC89.doublet.info$colnames.SC89.sce.
rownames(SC108.doublet.info)=SC108.doublet.info$colnames.SC108.sce.
rownames(SC143.doublet.info)=SC143.doublet.info$colnames.SC143.sce.
rownames(SC147.doublet.info)=SC147.doublet.info$colnames.SC147.sce.
rownames(SC151.doublet.info)=SC151.doublet.info$colnames.SC151.sce.



View(SC86@meta.data)

SC86=AddMetaData(SC86,SC86.doublet.info,col.name = c("colnames.SC86.sce.","Doublet_Identification_Call"))
SC89=AddMetaData(SC89,SC89.doublet.info,col.name = c("colnames.SC89.sce.","Doublet_Identification_Call"))
SC108=AddMetaData(SC108,SC108.doublet.info,col.name = c("colnames.SC108.sce.","Doublet_Identification_Call"))
SC143=AddMetaData(SC143,SC143.doublet.info,col.name = c("colnames.SC143.sce.","Doublet_Identification_Call"))
SC147=AddMetaData(SC147,SC147.doublet.info,col.name = c("colnames.SC147.sce.","Doublet_Identification_Call"))
SC151=AddMetaData(SC151,SC151.doublet.info,col.name = c("colnames.SC151.sce.","Doublet_Identification_Call"))


head(SC86@meta.data)
head(SC89@meta.data)


colnames(SC86@meta.data)[4]="Cell_id"
colnames(SC89@meta.data)[4]="Cell_id"
colnames(SC108@meta.data)[4]="Cell_id"
colnames(SC143@meta.data)[4]="Cell_id"
colnames(SC147@meta.data)[4]="Cell_id"
colnames(SC151@meta.data)[4]="Cell_id"


SC86@meta.data$Cell_Id=rownames(SC86@meta.data)
meta.SC86=SC86@meta.data
test_cell_id=0
for(index in 1:nrow(meta.SC86)){
  if(meta.SC86[index,4]==meta.SC86[index,6]){
    test_cell_id=test_cell_id+1
  }
}



SC86@meta.data$Sample_id="SC86"
SC89@meta.data$Sample_id="SC89"
SC108@meta.data$Sample_id="SC108"
SC143@meta.data$Sample_id="SC143"
SC147@meta.data$Sample_id="SC147"
SC151@meta.data$Sample_id="SC151"


SC86@meta.data$Phenotype="Control(Bharat)"
SC89@meta.data$Phenotype="Control(Bharat)"
SC108@meta.data$Phenotype="COVID19"
SC143@meta.data$Phenotype="Transplant Donor"
SC147@meta.data$Phenotype="transplant Recipient"
SC151@meta.data$Phenotype="COVID19"



bharat_stromal=merge(SC86, y = c(SC89,SC108,SC143,SC147,SC151), add.cell.ids = c("SC86", "SC89","SC108","SC143","SC147","SC151"), project = "Stromal")
rm(all_epi_bharat)
rm(SC86,SC89,SC108,SC143,SC147,SC151)
View(bharat_stromal@meta.data)
x=bharat_stromal@meta.data
x_incmpl=x[!complete.cases(x),]

habermann=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/All_Habermann_etal_samples_Seurat_object_after_various_QC_ONLY_USE_THIS.rds')


head(habermann@meta.data)


bharat_stromal <- NormalizeData(bharat_stromal, normalization.method = "LogNormalize", scale.factor = 10000,verbose = TRUE)
bharat_stromal <- FindVariableFeatures(bharat_stromal, selection.method = "vst", nfeatures = 2000,verbose = TRUE)
variable.genes=VariableFeatures(bharat_stromal,selection.method = "vst",verbose=TRUE)
bharat_stromal <- ScaleData(bharat_stromal,verbose = TRUE,features = variable.genes)
bharat_stromal <- RunPCA(bharat_stromal, features = variable.genes,verbose = TRUE,ndims.print = 1:50)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_stromal_Cells_Elbowplot_no_doublet_Removal.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

ElbowPlot(bharat_stromal)

dev.off()
vars.per.PC=bharat_stromal@reductions$pca@stdev*100/sum(bharat_stromal@reductions$pca@stdev)
vars.per.PC*100/sum(vars.per.PC)
sum(vars.per.PC[1:40])

bharat_stromal <- RunUMAP(bharat_stromal, dims=1:40,verbose = TRUE,reduction = "pca")


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_UMAP_to_see_doublet_effect.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
DimPlot(bharat_stromal, group.by = "Doublet_Identification_Call",raster = FALSE,reduction = "umap")
dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_UMAP_phenotype_stromal_Cells.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
DimPlot(bharat_stromal, group.by = "Phenotype",raster = FALSE,reduction = "umap")
dev.off()




bharat_stromal_no_doublet=subset(bharat_stromal, subset = Doublet_Identification_Call!="doublet")



bharat_stromal_no_doublet <- NormalizeData(bharat_stromal_no_doublet, normalization.method = "LogNormalize", scale.factor = 10000,verbose = TRUE)
bharat_stromal_no_doublet <- FindVariableFeatures(bharat_stromal_no_doublet, selection.method = "vst", nfeatures = 2000,verbose = TRUE)
variable.genes=VariableFeatures(bharat_stromal_no_doublet,selection.method = "vst",verbose=TRUE)
bharat_stromal_no_doublet <- ScaleData(bharat_stromal_no_doublet,verbose = TRUE,features = variable.genes)
bharat_stromal_no_doublet <- RunPCA(bharat_stromal_no_doublet, features = variable.genes,verbose = TRUE,ndims.print = 1:50)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_stromal_Cells_Elbowplot_doublet_removed.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

ElbowPlot(bharat_stromal_no_doublet)

dev.off()
vars.per.PC=bharat_stromal_no_doublet@reductions$pca@stdev*100/sum(bharat_stromal_no_doublet@reductions$pca@stdev)
vars.per.PC*100/sum(vars.per.PC)
sum(vars.per.PC[1:40])

bharat_stromal_no_doublet <- RunUMAP(bharat_stromal_no_doublet, dims=1:40,verbose = TRUE,reduction = "pca")


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_stromal_UMAP_doublet_removed.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
DimPlot(bharat_stromal_no_doublet, group.by = "Doublet_Identification_Call",raster = FALSE,reduction = "umap")
dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_stromal_UMAP_doublet_removed_phenotype.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
DimPlot(bharat_stromal_no_doublet, group.by = "Phenotype",raster = FALSE,reduction = "umap")
dev.off()

saveRDS(bharat_stromal_no_doublet,file='Bharat_Stromal_No_doublet_for_clustering.rds')

bharat_stromal_no_doublet=FindNeighbors(bharat_stromal_no_doublet,verbose = TRUE,dims = 1:40)

bharat_stromal_no_doublet <- FindClusters(bharat_stromal_no_doublet, resolution = 0.5,algorithm = 4,method = "igraph",verbose = TRUE) ## Leiden algorithm

table(Idents(bharat_stromal_no_doublet))



saveRDS(bharat_stromal_no_doublet,file="Bharat_Stromal_cells_leiden_clustered_for_further_analysis.rds")



bharat_stromal_no_doublet=BuildClusterTree(bharat_stromal_no_doublet,verbose = TRUE)
PlotClusterTree(bharat_stromal_no_doublet)


library(cowplot)
p2 <- DimPlot(bharat_stromal_no_doublet, reduction = "umap", label = TRUE, repel = TRUE,label.size = 9)
p3 = DimPlot(bharat_stromal_no_doublet, group.by = "Phenotype",raster = FALSE,reduction = "umap")



pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_stromal_UMAP_doublet_removed_leiden_clustered.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(p2)

dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_stromal_UMAP_doublet_removed_leiden_clustered_possibly_pheno_on_cluster.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(p3)

dev.off()


library(ggplot2)

goi=c("FBLN1", "MMP2", "LAMC3", "COX4I2", "KRT8", "MRC1", "SCARA5", "PDGFRA", "COL1A1", "HAS1", "ACTA2", "MYH11", "MKI67", "HBM")
heatmap.object=DoHeatmap(bharat_stromal_no_doublet, features = goi,group.by = "seurat_clusters",slot = 'data')
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_stromal_UMAP_doublet_removed_leiden_clustered_marker_gene_by_cluster_heatmap_blue_white_red.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)

dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_stromal_UMAP_doublet_removed_leiden_clustered_marker_gene_by_cluster_violin_plot.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

VlnPlot(bharat_stromal_no_doublet, features = goi,group.by = "seurat_clusters",flip=TRUE,log=TRUE)

dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_stromal_UMAP_doublet_removed_leiden_clustered_dendogram.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

PlotClusterTree(bharat_stromal_no_doublet)

dev.off()


Idents(bharat_stromal_no_doublet)

heatmap.object=DoHeatmap(bharat_stromal_no_doublet, features = goi,group.by = "seurat_clusters",slot = 'data',group.colors = c("blue","red"))  #+ scale_fill_gradient(colors = c("blue", "white", "red"))


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_stromal_UMAP_doublet_removed_leiden_clustered_marker_gene_by_cluster_heatmap.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)


dev.off()



#### Use the expression pattern of NM to combine clusters


heatmap.object=DoHeatmap(bharat_stromal_no_doublet, features = "HBM",group.by = "seurat_clusters",slot = 'data',group.colors = c("blue","red"))  #+ scale_fill_gradient(colors = c("blue", "white", "red"))

VlnPlot(bharat_stromal_no_doublet, features = "HBM",group.by = "seurat_clusters",flip=TRUE,log=FALSE)
VlnPlot(bharat_stromal_no_doublet, features = "COL1A1",group.by = "seurat_clusters",flip=TRUE,log=FALSE)
VlnPlot(bharat_stromal_no_doublet, features = "FBLN1",group.by = "seurat_clusters",flip=TRUE,log=FALSE)
VlnPlot(bharat_stromal_no_doublet, features = "MMP2",group.by = "seurat_clusters",flip=TRUE,log=FALSE)
VlnPlot(bharat_stromal_no_doublet, features = "HBM",group.by = "seurat_clusters",flip=TRUE,log=FALSE)


VlnPlot(bharat_stromal_no_doublet, features ="COL1A1" ,group.by = "seurat_clusters",flip=TRUE,log=FALSE)

VlnPlot(bharat_stromal_no_doublet, features = "HBM",group.by = "seurat_clusters",flip=TRUE,log=FALSE,idents = 10)


metat_after_clustering=bharat_stromal_no_doublet@meta.data



table(metat_after_clustering$Phenotype)

covid=subset(metat_after_clustering,metat_after_clustering$Sample_id=="SC108" | metat_after_clustering$Sample_id=="SC147" | metat_after_clustering$Sample_id=="SC151")
no_covid=subset(metat_after_clustering,metat_after_clustering$Sample_id=="SC86" | metat_after_clustering$Sample_id=="SC89" | metat_after_clustering$Sample_id=="SC143")
covid$Status="True"
no_covid$Status="False"
nrow(covid)+nrow(no_covid) == nrow(metat_after_clustering)

metat_after_clustering=rbind(covid,no_covid)
metat_after_clustering$Count_per_cell=1

ggplot(metat_after_clustering, aes(fill=Status, y=Count_per_cell, x=seurat_clusters)) + 
  geom_bar(position="stack",stat="identity") +
  ggtitle("Cluster by COVID19") +
  xlab("Cluster Id")



bharat_stromal_no_doublet@meta.data=metat_after_clustering
saveRDS(bharat_stromal_no_doublet,file="Bharat_Stromal_cells_leiden_clustered_for_further_analysis.rds")





genes_paper=c("ACTA2", "PLIN2", "TWIST1", "LTBP2", "VEGFA", "ITGA8", "ABL2", "MYC", "FAP","PDGFRA", "PDGFRB", "CSF1", "CD34", "FBN1", "FBLN2", "VIT", "NEBL", "HAS2", "CXCL14")


heatmap.object=DoHeatmap(bharat_stromal_no_doublet, features = genes_paper,group.by = "seurat_clusters",slot = 'data')
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_stromal_heatmap_paper_mentioned_genes.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)

dev.off()


