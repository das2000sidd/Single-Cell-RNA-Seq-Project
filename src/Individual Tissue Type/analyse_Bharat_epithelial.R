setwd("/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells")


library(Seurat)




SC85=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC85_seurat_min_3_cells_min_200_genes.rds')
SC88=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC88_seurat_min_3_cells_min_200_genes.rds')
SC107=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC107_seurat_min_3_cells_min_200_genes.rds')
SC141=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC141_seurat_min_3_cells_min_200_genes.rds')
SC145=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC145_seurat_min_3_cells_min_200_genes.rds')
SC152=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC152_seurat_min_3_cells_min_200_genes.rds')


library(scDblFinder)


SC85.sce <- as.SingleCellExperiment(SC85)
SC88.sce <- as.SingleCellExperiment(SC88)
SC107.sce <- as.SingleCellExperiment(SC107)
SC141.sce <- as.SingleCellExperiment(SC141)
SC145.sce <- as.SingleCellExperiment(SC145)
SC152.sce <- as.SingleCellExperiment(SC152)



SC85.sce=scDblFinder(SC85.sce,verbose = TRUE)
SC88.sce=scDblFinder(SC88.sce,verbose = TRUE)
SC107.sce=scDblFinder(SC107.sce,verbose = TRUE)
SC141.sce=scDblFinder(SC141.sce,verbose = TRUE)
SC145.sce=scDblFinder(SC145.sce,verbose = TRUE)
SC152.sce=scDblFinder(SC152.sce,verbose = TRUE)




SC85.doublet.info=data.frame(colnames(SC85.sce),SC85.sce$scDblFinder.class)
SC88.doublet.info=data.frame(colnames(SC88.sce),SC88.sce$scDblFinder.class)
SC107.doublet.info=data.frame(colnames(SC107.sce),SC107.sce$scDblFinder.class)
SC141.doublet.info=data.frame(colnames(SC141.sce),SC141.sce$scDblFinder.class)
SC145.doublet.info=data.frame(colnames(SC145.sce),SC145.sce$scDblFinder.class)
SC152.doublet.info=data.frame(colnames(SC152.sce),SC152.sce$scDblFinder.class)




rownames(SC85.doublet.info)=SC85.doublet.info$colnames.SC85.sce.
rownames(SC88.doublet.info)=SC88.doublet.info$colnames.SC88.sce.
rownames(SC107.doublet.info)=SC107.doublet.info$colnames.SC107.sce.
rownames(SC141.doublet.info)=SC141.doublet.info$colnames.SC141.sce.
rownames(SC145.doublet.info)=SC145.doublet.info$colnames.SC145.sce.
rownames(SC152.doublet.info)=SC152.doublet.info$colnames.SC152.sce.




SC85=AddMetaData(SC85,SC85.doublet.info,col.name = c("colnames.SC85.sce.","Doublet_Identification_Call"))
SC88=AddMetaData(SC88,SC88.doublet.info,col.name = c("colnames.SC88.sce.","Doublet_Identification_Call"))
SC107=AddMetaData(SC107,SC107.doublet.info,col.name = c("colnames.SC107.sce.","Doublet_Identification_Call"))
SC141=AddMetaData(SC141,SC141.doublet.info,col.name = c("colnames.SC141.sce.","Doublet_Identification_Call"))
SC145=AddMetaData(SC145,SC145.doublet.info,col.name = c("colnames.SC145.sce.","Doublet_Identification_Call"))
SC152=AddMetaData(SC152,SC152.doublet.info,col.name = c("colnames.SC152.sce.","Doublet_Identification_Call"))


head(SC85@meta.data)
head(SC88@meta.data)


colnames(SC85@meta.data)[4]="Cell_id"
colnames(SC88@meta.data)[4]="Cell_id"
colnames(SC107@meta.data)[4]="Cell_id"
colnames(SC141@meta.data)[4]="Cell_id"
colnames(SC145@meta.data)[4]="Cell_id"
colnames(SC152@meta.data)[4]="Cell_id"



SC85@meta.data$Sample_id="SC85"
SC88@meta.data$Sample_id="SC88"
SC107@meta.data$Sample_id="SC107"
SC141@meta.data$Sample_id="SC141"
SC145@meta.data$Sample_id="SC145"
SC152@meta.data$Sample_id="SC152"


SC85@meta.data$Phenotype="Control(Bharat)/Epithelial"
SC88@meta.data$Phenotype="Control(Bharat)/Epithelial"
SC107@meta.data$Phenotype="COVID19/Epithelial"
SC141@meta.data$Phenotype="Transplant Donor/Epithelial"
SC145@meta.data$Phenotype="Transplant Recipient/Epithelial"
SC152@meta.data$Phenotype="COVID19/Epithelial"



bharat_epithelial=merge(SC85, y = c(SC88,SC107,SC141,SC145,SC152), add.cell.ids = c("SC85", "SC88","SC107","SC141","SC145","SC152"), project = "Epithelial")
rm(SC85,SC88,SC107,SC141,SC145,SC152)
rm(SC85.sce,SC88.sce,SC107.sce,SC141.sce,SC145.sce,SC152.sce)

View(bharat_epithelial@meta.data)
x=bharat_epithelial@meta.data
x_incmpl=x[!complete.cases(x),]

habermann=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/All_Habermann_etal_samples_Seurat_object_after_various_QC_ONLY_USE_THIS.rds')


head(habermann@meta.data)


bharat_epithelial <- NormalizeData(bharat_epithelial, normalization.method = "LogNormalize", scale.factor = 10000,verbose = TRUE)
bharat_epithelial <- FindVariableFeatures(bharat_epithelial, selection.method = "vst", nfeatures = 2000,verbose = TRUE)
variable.genes=VariableFeatures(bharat_epithelial,selection.method = "vst",verbose=TRUE)
bharat_epithelial <- ScaleData(bharat_epithelial,verbose = TRUE,features = variable.genes)
bharat_epithelial <- RunPCA(bharat_epithelial, features = variable.genes,verbose = TRUE,ndims.print = 1:50)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_Cells_Elbowplot_no_doublet_Removal.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

ElbowPlot(bharat_epithelial,ndims=40,reduction = "pca")

dev.off()
vars.per.PC=bharat_epithelial@reductions$pca@stdev*100/sum(bharat_epithelial@reductions$pca@stdev)
vars.per.PC*100/sum(vars.per.PC)
sum(vars.per.PC[1:40])

bharat_epithelial <- RunUMAP(bharat_epithelial, dims=1:40,verbose = TRUE,reduction = "pca")


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_UMAP_epithelial_to_see_doublet_effect.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
DimPlot(bharat_epithelial, group.by = "Doublet_Identification_Call",raster = FALSE,reduction = "umap")
dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_UMAP_phenotype_epithelial_Cells.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
DimPlot(bharat_epithelial, group.by = "Phenotype",raster = FALSE,reduction = "umap")
dev.off()

bharat_epithelial[["percent.mt"]] <- PercentageFeatureSet(bharat_epithelial, pattern = "^MT-")
VlnPlot(bharat_epithelial, features = c("percent.mt"),group.by = "Sample_id")


library(cowplot)
bharat_epithelial_no_doublet=subset(bharat_epithelial, subset = Doublet_Identification_Call!="doublet")
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_Cells_no_doublet_mito_percent.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

VlnPlot(bharat_epithelial_no_doublet, features = c("percent.mt"),group.by = "Sample_id")
dev.off()


bharat_epithelial_no_SC88_no_doublet=subset(bharat_epithelial_no_doublet,subset = Sample_id!="SC88")


bharat_epithelial_no_SC88_no_doublet <- NormalizeData(bharat_epithelial_no_SC88_no_doublet, normalization.method = "LogNormalize", scale.factor = 10000,verbose = TRUE)
bharat_epithelial_no_SC88_no_doublet <- FindVariableFeatures(bharat_epithelial_no_SC88_no_doublet, selection.method = "vst", nfeatures = 2000,verbose = TRUE)
variable.genes=VariableFeatures(bharat_epithelial_no_SC88_no_doublet,selection.method = "vst",verbose=TRUE)
write.table(variable.genes,file="Bharat_epithelial_variable_genes_after_doublet__SC88_removal.txt",col.names = F,row.names = F,sep="\t",quote = F)
bharat_epithelial_no_SC88_no_doublet <- ScaleData(bharat_epithelial_no_SC88_no_doublet,verbose = TRUE,features = variable.genes)
bharat_epithelial_no_SC88_no_doublet <- RunPCA(bharat_epithelial_no_SC88_no_doublet, features = variable.genes,verbose = TRUE,ndims.print = 1:50)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_Cells_Elbowplot_doublet_removed.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

ElbowPlot(bharat_epithelial_no_SC88_no_doublet)

dev.off()
vars.per.PC=bharat_epithelial_no_SC88_no_doublet@reductions$pca@stdev*100/sum(bharat_epithelial_no_SC88_no_doublet@reductions$pca@stdev)
vars.per.PC*100/sum(vars.per.PC)
sum(vars.per.PC[1:40])

bharat_epithelial_no_SC88_no_doublet <- RunUMAP(bharat_epithelial_no_SC88_no_doublet, dims=1:40,verbose = TRUE,reduction = "pca")


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_SC88_removed.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
x=DimPlot(bharat_epithelial_no_SC88_no_doublet, group.by = "Doublet_Identification_Call",raster = FALSE,reduction = "umap")
dev.off()

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_SC88_removed_phenotype.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
DimPlot(bharat_epithelial_no_SC88_no_doublet, group.by = "Phenotype",raster = FALSE,reduction = "umap")
dev.off()

saveRDS(bharat_epithelial_no_SC88_no_doublet,file='bharat_epithelial_No_doublet_No_SC88_for_clustering.rds')

### RUN ON CLUSTER
bharat_epithelial_no_SC88_no_doublet=FindNeighbors(bharat_epithelial_no_SC88_no_doublet,verbose = TRUE,dims = 1:40)

bharat_epithelial_no_SC88_no_doublet <- FindClusters(bharat_epithelial_no_SC88_no_doublet, resolution = 0.5,algorithm = 4,method = "igraph",verbose = TRUE) ## Leiden algorithm
### RUN ON CLUSTER

table(Idents(bharat_epithelial_no_SC88_no_doublet))



saveRDS(bharat_epithelial_no_SC88_no_doublet,file="bharat_epithelial_cells_leiden_clustered_for_further_analysis.rds")


bharat_epithelial_no_SC88_no_doublet=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_No_doublet_with_cluster_ids.rds')

bharat_epithelial_no_SC88_no_doublet=BuildClusterTree(bharat_epithelial_no_SC88_no_doublet,verbose = TRUE)
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_No_doublet_cluster_dendogram.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

PlotClusterTree(bharat_epithelial_no_SC88_no_doublet)
dev.off()

library(cowplot)
p2 <- DimPlot(bharat_epithelial_no_SC88_no_doublet, reduction = "umap", label = TRUE, repel = TRUE,label.size = 9)
p3 = DimPlot(bharat_epithelial_no_SC88_no_doublet, group.by = "Phenotype",raster = FALSE,reduction = "umap")



pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_removed_leiden_clustered.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(p2)

dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_removed_leiden_clustered_possibly_pheno_on_cluster.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(p3)

dev.off()


library(ggplot2)

goi=c("SCGB3A1", "MUC5B", "SFTPC", "LAMP3", "HOPX", "AGER", "KRT5","FOXJ1", "S100A2", "KRT17", "CLDN18", "PEG10", "CLIC5", "MMP10")
heatmap.object=DoHeatmap(bharat_epithelial_no_SC88_no_doublet, features = goi,group.by = "seurat_clusters",slot = 'data')
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_removed_leiden_clustered_marker_gene_by_cluster_heatmap_blue_white_red.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)

dev.off()


goi=c("TP63", "KRT5", "KRT17", "LAMB3", "LAMC2", "VIM", "CHD2", "FN1", "COL1A1", "TNC", "HMGA2", "CDKN1A", "CDKN2A", "CCND1", "CCND2", "MDM2", "SERPINE1", "TGFB1", "ITGAV", "ITGB6")


heatmap.object=DoHeatmap(bharat_epithelial_no_SC88_no_doublet, features = goi,group.by = "seurat_clusters",slot = 'data')
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_removed_leiden_clustered_marker_gene_from_Habermann_to_predict_cell_type.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)

dev.off()



View(bharat_epithelial_no_SC88_no_doublet@meta.data)
meta_no_doublet=bharat_epithelial_no_SC88_no_doublet@meta.data
meta_no_doublet$Count_per_cell=1
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_removed_leiden_clustered_per_cluster_COVID_status_count.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

ggplot(meta_no_doublet, aes(fill=Phenotype, y=Count_per_cell, x=seurat_clusters)) + 
  geom_bar(position="stack",stat="identity") +
  ggtitle("Cluster by COVID19") +
  xlab("Cluster Id")
dev.off()


goi_habermann=c("TP63", "KRT5", "KRT17", "LAMB3", "LAMC2", 
                "VIM", "CHD2", "FN1", "COL1A1", "TNC", 
                "HMGA2", "CDKN1A", "CDKN2A", "CCND1", "CCND2", 
                "MDM2", "SERPINE1", "TGFB1", "ITGAV", "ITGB6")


heatmap.object=DoHeatmap(bharat_epithelial_no_SC88_no_doublet, features = goi_habermann,group.by = "seurat_clusters",slot = 'data')
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_removed_leiden_clustered_genes_heatmap_meitioned_in_Habermann.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)
dev.off()


dotplot.object=DotPlot(bharat_epithelial_no_SC88_no_doublet, features = goi_habermann,group.by = "seurat_clusters",cluster.idents = TRUE,scale.by = "size",idents = 1:10)
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_violin_plot_doublet_removed_leiden_clustered_genes_dotplot_meitioned_in_Habermann_C1_to_C10.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 20)


plot_grid(dotplot.object)
dev.off()


p3 = DimPlot(bharat_epithelial_no_SC88_no_doublet, group.by = "Phenotype",raster = FALSE,reduction = "umap")

metadata_no_dub=bharat_epithelial_no_SC88_no_doublet@meta.data
covid=subset(metadata_no_dub,metadata_no_dub$Phenotype=="COVID19/Epithelial" | metadata_no_dub$Phenotype=="Transplant Recipient/Epithelial")
not_covid=subset(metadata_no_dub,metadata_no_dub$Phenotype!="Control(Bharat)/Epithelial" & metadata_no_dub$Phenotype!="Transplant Donor/Epithelial")
covid$Covid_Status="Yes"
not_covid$Covid_Status="No"
metadata_no_dub=rbind(covid,not_covid)
metadata_no_dub$Count_per_cell=1

ggplot(metadata_no_dub, aes(fill=Covid_Status, y=Count_per_cell, x=seurat_clusters)) + 
  geom_bar(position="stack",stat="identity") +
  ggtitle("Cluster by COVID19") +
  xlab("Cluster Id")


bharat_epithelial_no_SC88@meta.data=metadata_no_dub

p3 = DimPlot(bharat_epithelial_no_SC88_no_doublet, group.by = "Covid_Status",raster = FALSE,reduction = "umap")

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_removed_leiden_clustered_covid_status.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(p3)
dev.off()



pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_removed_leiden_clustered_marker_gene_by_cluster_violin_plot.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

VlnPlot(bharat_epithelial_no_SC88_no_doublet, features = goi,group.by = "seurat_clusters",flip=TRUE,log=TRUE)

dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_removed_leiden_clustered_dendogram.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

PlotClusterTree(bharat_epithelial_no_SC88_no_doublet)

dev.off()





Idents(bharat_epithelial_no_doublet)

heatmap.object=DoHeatmap(bharat_epithelial_no_SC88_no_doublet, features = goi,group.by = "seurat_clusters",slot = 'data',group.colors = c("blue","red"))  #+ scale_fill_gradient(colors = c("blue", "white", "red"))
Idents(bharat_epithelial_no_SC88_no_doublet)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_UMAP_doublet_removed_leiden_clustered_marker_gene_by_cluster_heatmap.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)


dev.off()



#### Use the expression pattern of NM to combine clusters




metadata_no_dub=bharat_epithelial_no_SC88_no_doublet@meta.data
table(metadata_no_dub$Phenotype)

covid=subset(metadata_no_dub,metadata_no_dub$Phenotype=="COVID19/Epithelial" | metadata_no_dub$Phenotype=="Transplant Recipient/Epithelial")
no_covid=subset(metadata_no_dub,metadata_no_dub$Phenotype=="Control(Bharat)/Epithelial" | metadata_no_dub$Phenotype=="Transplant Donor/Epithelial")
covid$Status="True"
no_covid$Status="False"
nrow(covid)+nrow(no_covid) == nrow(metadata_no_dub)

metadata_no_dub=rbind(covid,no_covid)
metadata_no_dub$Count_per_cell=1

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_post_clustering_phenotype_cell_count_per_cluster.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

ggplot(metadata_no_dub, aes(fill=Status, y=Count_per_cell, x=seurat_clusters)) + 
  geom_bar(position="stack",stat="identity") +
  ggtitle("Cluster by COVID19") +
  xlab("Cluster Id")
dev.off()




bharat_epithelial_no_SC88_no_doublet@meta.data=metat_after_clustering
saveRDS(bharat_epithelial_no_SC88_no_doublet,file="bharat_epithelial_cells_leiden_clustered_for_further_analysis.rds")





genes_paper=c()


heatmap.object=DoHeatmap(bharat_epithelial_no_SC88, features = genes_paper,group.by = "seurat_clusters",slot = 'data')
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_heatmap_paper_mentioned_genes.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)

dev.off()


habermann_epithelial_genes=c("CCL2", "VCAN", "VEGFA", "FCN1", "SELENOP", "SLC40A1","SPP1", "IL1RN", "MMP9", "CHI3L1", "MARCKS", "PLA2G7","PPARG", "MARCO", "MRC1", "INHBA","SLC9B","CKB", "SIGLEC15", "AK5")

heatmap.object=DoHeatmap(bharat_epithelial_no_SC88, features = habermann_epithelial_genes,group.by = "seurat_clusters",slot = 'data')
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_heatmap_habermann_paper_mentioned_genes.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)

dev.off()


saveRDS(bharat_epithelial_no_SC88_no_doublet,file="bharat_epithelial_cells_leiden_clustered_no_doublet_no_SC88_for_Habermann_integration.rds")





integrated.clustered=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_cluster_id_habermann_cell_type_epithelial_integrated_with_cluster_ids.rds')


head(integrated.clustered@meta.data)



table(integrated.clustered@meta.data$integrated_snn_res.0.5) ## 22
table(integrated.clustered@meta.data$integrated_snn_res.0.3) ## 17
table(integrated.clustered@meta.data$integrated_snn_res.0.8) ## 31
table(integrated.clustered@meta.data$integrated_snn_res.1.2) ## 36


DimPlot(integrated.clustered,raster = FALSE,reduction = "umap")
Idents(integrated.clustered)=integrated.clustered@meta.data$integrated_snn_res.0.3

library(Seurat)
integrated.clustered=BuildClusterTree(integrated.clustered,dims = 1:40)
PlotClusterTree(integrated.clustered,direction = "rightwards")

original_clustering_Table=integrated.clustered@meta.data
write.csv(original_clustering_Table,file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Epithelial_integrated_post_clustering_table.csv",col.names = T,row.names = T,quote = F)

integrated.clustered@active.assay


DefaultAssay(integrated.clustered) <- "RNA"

comp3_10 <- FindMarkers(integrated.clustered, ident.1 =3, ident.2 = 10, verbose = TRUE,test.use = "wilcox")

comp3_10_different=subset(comp3_10,comp3_10$avg_log2FC>1 & comp3_10$p_val_adj < 0.01)
comp3_10_different$Gene=rownames(comp3_10_different)

comp3_10_conserved <- FindConservedMarkers(integrated.clustered, ident.1 =3, ident.2 = 10, verbose = TRUE,grouping.var = "")


marker_genes=c("SCGB3A1", "MUC5B", "SFTPC", "LAMP3", "HOPX", "AGER", "KRT5", 
        "FOXJ1", "S100A2", "KRT17", "CLDN18", "PEG10", "CLIC5", "MMP10")
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_habermann_epithelial_integrated_with_epithelial_marker_Genes.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

DoHeatmap(integrated.clustered,features = marker_genes,group.by = "ident",angle = 0.5,raster = FALSE,slot = 'data')
dev.off()

bharat_epithelial_no_SC88_no_doublet=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_No_doublet_with_cluster_ids.rds')

Idents(bharat_epithelial_no_SC88_no_doublet)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_with_epithelial_marker_Genes.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)


DoHeatmap(bharat_epithelial_no_SC88_no_doublet,features = marker_genes,group.by = "ident",angle = 0.5,raster = FALSE,slot = 'data')
dev.off()

genes_habermann_Celltype=c("TP63", "KRT5", "KRT17", "LAMB3", "LAMC2", 
                           "VIM", "CHD2", "FN1", "COL1A1", "TNC", 
                           "HMGA2", "CDKN1A", "CDKN2A", "CCND1", "CCND2", 
                           "MDM2", "SERPINE1", "TGFB1", "ITGAV", "ITGB6")

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_clusters_with_Habermann_Celltype_marker_Genes_scaled.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)


DoHeatmap(bharat_epithelial_no_SC88_no_doublet,features = genes_habermann_Celltype,group.by = "ident",angle = 0.5,raster = FALSE,slot = 'scale.data')
dev.off()

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_clusters_violin_plot_with_Habermann_Celltype_marker_Genes.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

VlnPlot(bharat_epithelial_no_SC88_no_doublet,features = genes_habermann_Celltype,log = TRUE)
dev.off()

cumsum(prop.table(table(Idents(bharat_epithelial_no_SC88_no_doublet)))*100)

original_clustering=bharat_epithelial_no_SC88_no_doublet@meta.data
  
## keep till 13 clusters which is 73% of data

subsetted_data=subset(bharat_epithelial_no_SC88_no_doublet, idents= c(1:13))



prop.table(table(Idents(bharat_epithelial_no_SC88_no_doublet)))*100




pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_epithelial_larest_13_clusters_heatmap_with_Habermann_Celltype_marker_Genes.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

DoHeatmap(subsetted_data,features = genes_habermann_Celltype,group.by = "ident",angle = 0.5,raster = FALSE,slot = 'data')
dev.off()




Average_by_cluster=AverageExpression(bharat_epithelial_no_SC88_no_doublet,features = genes_habermann_Celltype,verbose = TRUE,group.by = 'ident')


avg_per_cluster=Average_by_cluster$RNA
avg_per_cluster_t=t(avg_per_cluster)
avg_per_cluster_t=as.data.frame(avg_per_cluster_t)
avg_per_cluster_t$cluster_id=rownames(avg_per_cluster_t)
boxplot(avg_per_cluster_t$TP63~avg_per_cluster_t$cluster_id)





  
write.csv(original_clustering,file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Epithelial_Bharat_post_doublet_removal_clustering_metadat.csv",col.names = T,row.names = T,quote = F)

original_clustering$seurat_clusters[original_clustering$seurat_clusters==10]=3
original_clustering$seurat_clusters[original_clustering$seurat_clusters==14]=2
original_clustering$seurat_clusters[original_clustering$seurat_clusters==11]=5
original_clustering$seurat_clusters[original_clustering$seurat_clusters==15]=8
original_clustering$seurat_clusters[original_clustering$seurat_clusters==17]=6
original_clustering$seurat_clusters[original_clustering$seurat_clusters==13]=7

table(original_clustering$RNA_snn_res.1.2)

bharat_epithelial_no_SC88_no_doublet=BuildClusterTree(bharat_epithelial_no_SC88_no_doublet,dims = 1:40)
PlotClusterTree(bharat_epithelial_no_SC88_no_doublet,direction = "rightwards")


FeaturePlot(bharat_epithelial_no_SC88_no_doublet, features = genes_habermann_Celltype)
FeaturePlot(bharat_epithelial_no_SC88_no_doublet, features = "KRT5")
FeaturePlot(bharat_epithelial_no_SC88_no_doublet, features = "KRT17")
FeaturePlot(bharat_epithelial_no_SC88_no_doublet, features = "LAMB3")
FeaturePlot(bharat_epithelial_no_SC88_no_doublet, features = "LAMC2")
FeaturePlot(bharat_epithelial_no_SC88_no_doublet, features = "VIM")
FeaturePlot(bharat_epithelial_no_SC88_no_doublet, features = "CHD2")
FeaturePlot(bharat_epithelial_no_SC88_no_doublet, features = c("LAMB3","CCND1","ITGB6"))
FeaturePlot(bharat_epithelial_no_SC88_no_doublet, features = c("LAMB3","VIM","CHD2","CDKN1A","CCND1","ITGAV","ITGB6"))
FeaturePlot(bharat_epithelial_no_SC88_no_doublet, features = c("KRT5","KRT17"))
FeaturePlot(bharat_epithelial_no_SC88_no_doublet, features = c("LAMB3","CCND1","ITGB6")) ## AT1
VlnPlot(bharat_epithelial_no_SC88_no_doublet, features = c("LAMB3","LAMC2","ITGB6")) ## AT1


library(cowplot)
p3 = DimPlot(bharat_epithelial_no_SC88_no_doublet, group.by = "ident",raster = FALSE,reduction = "umap",label = TRUE)
plot_grid(p3)

table(Idents(bharat_epithelial_no_SC88_no_doublet))





DoHeatmap(subsetted_data,features = genes_habermann_Celltype,group.by = "ident",angle = 0.5,raster = FALSE,slot = 'scale.data')





immune=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_No_doublet_with_cluster_ids.rds')
View(immune@meta.data)
immune.markers=c("MRC1", "MARCO", "FABP4", "SPP1", "CCL2", "CD3E", "CD4", "CD8A", 
                 "JCHAIN", "CD79A", "KIT", "CCR7", "CD1C", "FCER1A", "CLEC9A", "CLEC4C",
                 "VWF")
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Testing_immune_marker.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

DoHeatmap(immune,features = immune.markers,group.by = "ident",angle = 0.5,raster = FALSE,slot = 'data')

dev.off()




habermann_genes=c("CCL2", "VCAN", "VEGFA", "FCN1", "SELENOP", "SLC40A1",
"SPP1", "IL1RN", "MMP9", "CHI3L1", "MARCKS", "PLA2G7",
"PPARG", "MARCO", "MRC1", "INHBA",
"SLC9B2", "CKB", "SIGLEC15", "AK5")


DoHeatmap(immune,features = habermann_genes,group.by = "ident",angle = 0.5,raster = FALSE,slot = 'data')
DoHeatmap(immune,features = "CD23",group.by = "ident",angle = 0.5,raster = FALSE,slot = 'data')


rownames(immune)[grep("CD203",rownames(immune))]


immune_meta=immune@meta.data
immune_meta$seurat_clusters=as.character(immune_meta$seurat_clusters)
immune_meta$Predicted_Cell_Type=immune_meta$seurat_clusters

immune_meta$Predicted_Cell_Type[immune_meta$seurat_clusters=="9"]="MoAM-1"
immune_meta$Predicted_Cell_Type[immune_meta$seurat_clusters=="5"]="TRAM-1"
immune_meta$Predicted_Cell_Type[immune_meta$seurat_clusters=="3" | immune_meta$seurat_clusters=="2"]="MoAM-2"
immune_meta$Predicted_Cell_Type[immune_meta$seurat_clusters=="6"]="Plasma cells"
immune_meta$Predicted_Cell_Type[immune_meta$seurat_clusters=="1" | immune_meta$seurat_clusters=="4"]="MoAM-3"
immune_meta$Predicted_Cell_Type[immune_meta$seurat_clusters=="12"]="T cells"
immune_meta$Predicted_Cell_Type[immune_meta$seurat_clusters=="22"]="B cells"

immune@meta.data=immune_meta

saveRDS(immune,file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_No_doublet_with_cluster_ids_and_Bharat_Cell_Types.rds")

immune_meta$Count=1

library(ggplot2)
ggplot(immune_meta, aes(fill=Predicted_Cell_Type, y=Count, x=seurat_clusters)) + 
  geom_bar(position="stack",stat="identity") +
  ggtitle("Per Cluster Cell COunt") +
  xlab("Cluster Id")


prop.table(table(immune_meta$Predicted_Cell_Type))*100


