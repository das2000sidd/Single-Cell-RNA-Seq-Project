setwd("/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells")


library(Seurat)




SC84=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC86_seurat_min_3_cells_min_200_genes.rds')
SC87=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC87_seurat_min_3_cells_min_200_genes.rds')
SC109=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC109_seurat_min_3_cells_min_200_genes.rds')
SC142=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC142_seurat_min_3_cells_min_200_genes.rds')
SC144=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC144_seurat_min_3_cells_min_200_genes.rds')
SC146=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC146_seurat_min_3_cells_min_200_genes.rds')
SC148=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC148_seurat_min_3_cells_min_200_genes.rds')
SC149=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC149_seurat_min_3_cells_min_200_genes.rds')
SC150=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC150_seurat_min_3_cells_min_200_genes.rds')


## SC86 is Donor1_3
## SC89 is Donor1_6
## SC108 is PB1_2
## SC143 is TX_Donor_3
## SC147 is TX_Recipient_3
## SC151 is PB2_3


library(scDblFinder)


SC84.sce <- as.SingleCellExperiment(SC84)
SC87.sce <- as.SingleCellExperiment(SC87)
SC109.sce <- as.SingleCellExperiment(SC109)
SC142.sce <- as.SingleCellExperiment(SC142)
SC144.sce <- as.SingleCellExperiment(SC144)
SC146.sce <- as.SingleCellExperiment(SC146)
SC148.sce <- as.SingleCellExperiment(SC148)
SC149.sce <- as.SingleCellExperiment(SC149)
SC150.sce <- as.SingleCellExperiment(SC150)



SC84.sce=scDblFinder(SC84.sce,verbose = TRUE)
SC87.sce=scDblFinder(SC87.sce,verbose = TRUE)
SC109.sce=scDblFinder(SC109.sce,verbose = TRUE)
SC142.sce=scDblFinder(SC142.sce,verbose = TRUE)
SC144.sce=scDblFinder(SC144.sce,verbose = TRUE)
SC146.sce=scDblFinder(SC146.sce,verbose = TRUE)
SC148.sce=scDblFinder(SC148.sce,verbose = TRUE)
SC149.sce=scDblFinder(SC149.sce,verbose = TRUE)
SC150.sce=scDblFinder(SC150.sce,verbose = TRUE)


#SC86.seurat.with.doublet=as.Seurat(SC86.sce)

SC84.doublet.info=data.frame(colnames(SC84.sce),SC84.sce$scDblFinder.class)
SC87.doublet.info=data.frame(colnames(SC87.sce),SC87.sce$scDblFinder.class)
SC109.doublet.info=data.frame(colnames(SC109.sce),SC109.sce$scDblFinder.class)
SC142.doublet.info=data.frame(colnames(SC142.sce),SC142.sce$scDblFinder.class)
SC144.doublet.info=data.frame(colnames(SC144.sce),SC144.sce$scDblFinder.class)
SC146.doublet.info=data.frame(colnames(SC146.sce),SC146.sce$scDblFinder.class)
SC148.doublet.info=data.frame(colnames(SC148.sce),SC148.sce$scDblFinder.class)
SC149.doublet.info=data.frame(colnames(SC149.sce),SC149.sce$scDblFinder.class)
SC150.doublet.info=data.frame(colnames(SC150.sce),SC150.sce$scDblFinder.class)


rownames(SC84.doublet.info)=SC84.doublet.info$colnames.SC84.sce.
rownames(SC87.doublet.info)=SC87.doublet.info$colnames.SC87.sce.
rownames(SC109.doublet.info)=SC109.doublet.info$colnames.SC109.sce.
rownames(SC142.doublet.info)=SC142.doublet.info$colnames.SC142.sce.
rownames(SC144.doublet.info)=SC144.doublet.info$colnames.SC144.sce.
rownames(SC146.doublet.info)=SC146.doublet.info$colnames.SC146.sce.
rownames(SC148.doublet.info)=SC148.doublet.info$colnames.SC148.sce.
rownames(SC149.doublet.info)=SC149.doublet.info$colnames.SC149.sce.
rownames(SC150.doublet.info)=SC150.doublet.info$colnames.SC150.sce.



View(SC84@meta.data)

SC84=AddMetaData(SC84,SC84.doublet.info,col.name = c("colnames.SC84.sce.","Doublet_Identification_Call"))
SC87=AddMetaData(SC87,SC87.doublet.info,col.name = c("colnames.SC87.sce.","Doublet_Identification_Call"))
SC109=AddMetaData(SC109,SC109.doublet.info,col.name = c("colnames.SC109.sce.","Doublet_Identification_Call"))
SC142=AddMetaData(SC142,SC142.doublet.info,col.name = c("colnames.SC142.sce.","Doublet_Identification_Call"))
SC144=AddMetaData(SC144,SC144.doublet.info,col.name = c("colnames.SC144.sce.","Doublet_Identification_Call"))
SC146=AddMetaData(SC146,SC146.doublet.info,col.name = c("colnames.SC146.sce.","Doublet_Identification_Call"))
SC148=AddMetaData(SC148,SC148.doublet.info,col.name = c("colnames.SC148.sce.","Doublet_Identification_Call"))
SC149=AddMetaData(SC149,SC149.doublet.info,col.name = c("colnames.SC149.sce.","Doublet_Identification_Call"))
SC150=AddMetaData(SC150,SC150.doublet.info,col.name = c("colnames.SC150.sce.","Doublet_Identification_Call"))


head(SC86@meta.data)
head(SC89@meta.data)


colnames(SC84@meta.data)[4]="Cell_id"
colnames(SC87@meta.data)[4]="Cell_id"
colnames(SC109@meta.data)[4]="Cell_id"
colnames(SC142@meta.data)[4]="Cell_id"
colnames(SC144@meta.data)[4]="Cell_id"
colnames(SC146@meta.data)[4]="Cell_id"
colnames(SC148@meta.data)[4]="Cell_id"
colnames(SC149@meta.data)[4]="Cell_id"
colnames(SC150@meta.data)[4]="Cell_id"




SC84@meta.data$Sample_id="SC84"
SC87@meta.data$Sample_id="SC87"
SC109@meta.data$Sample_id="SC109"
SC142@meta.data$Sample_id="SC142"
SC144@meta.data$Sample_id="SC144"
SC146@meta.data$Sample_id="SC146"
SC148@meta.data$Sample_id="SC148"
SC149@meta.data$Sample_id="SC149"
SC150@meta.data$Sample_id="SC150"


SC84@meta.data$Phenotype="Control(Bharat)/Myeloid"
SC87@meta.data$Phenotype="Control(Bharat)/Myeloid"
SC109@meta.data$Phenotype="COVID19/Myeloid"
SC142@meta.data$Phenotype="Transplant Donor/Myeloid"
SC144@meta.data$Phenotype="Transplant Donor/CD31"
SC146@meta.data$Phenotype="Transplant Recipient/Myeloid"
SC148@meta.data$Phenotype="Transplant Recipient/CD31"
SC149@meta.data$Phenotype="COVID19/Myeloid"
SC150@meta.data$Phenotype="COVID19/CD31"



bharat_immune=merge(SC84, y = c(SC87,SC109,SC142,SC144,SC146,SC148,SC149,SC150), add.cell.ids = c("SC84", "SC87","SC109","SC142","SC144","SC146","SC148","SC149","SC150"), project = "Immune")
rm(SC84,SC87,SC109,SC142,SC144,SC146,SC148,SC149,SC150)
rm(SC84.sce,SC87.sce,SC109.sce,SC142.sce,SC144.sce,SC146.sce,SC148.sce,SC149.sce,SC150.sce)

View(bharat_immune@meta.data)
x=bharat_immune@meta.data
x_incmpl=x[!complete.cases(x),]

habermann=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/All_Habermann_etal_samples_Seurat_object_after_various_QC_ONLY_USE_THIS.rds')


head(habermann@meta.data)


bharat_immune <- NormalizeData(bharat_immune, normalization.method = "LogNormalize", scale.factor = 10000,verbose = TRUE)
bharat_immune <- FindVariableFeatures(bharat_immune, selection.method = "vst", nfeatures = 2000,verbose = TRUE)
variable.genes=VariableFeatures(bharat_immune,selection.method = "vst",verbose=TRUE)
bharat_immune <- ScaleData(bharat_immune,verbose = TRUE,features = variable.genes)
bharat_immune <- RunPCA(bharat_immune, features = variable.genes,verbose = TRUE,ndims.print = 1:50)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_Cells_Elbowplot_no_doublet_Removal.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

ElbowPlot(bharat_immune)

dev.off()
vars.per.PC=bharat_immune@reductions$pca@stdev*100/sum(bharat_immune@reductions$pca@stdev)
vars.per.PC*100/sum(vars.per.PC)
sum(vars.per.PC[1:40])

bharat_immune <- RunUMAP(bharat_immune, dims=1:40,verbose = TRUE,reduction = "pca")


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_UMAP_immune_to_see_doublet_effect.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
DimPlot(bharat_immune, group.by = "Doublet_Identification_Call",raster = FALSE,reduction = "umap")
dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_UMAP_phenotype_immune_Cells.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
DimPlot(bharat_immune, group.by = "Phenotype",raster = FALSE,reduction = "umap")
dev.off()

bharat_immune[["percent.mt"]] <- PercentageFeatureSet(bharat_immune, pattern = "^MT-")
VlnPlot(bharat_immune, features = c("percent.mt"),group.by = "Sample_id")


library(cowplot)
bharat_immune_no_doublet=subset(bharat_immune, subset = Doublet_Identification_Call!="doublet")
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_Cells_no_doublet_mito_percent.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

VlnPlot(bharat_immune_no_doublet, features = c("percent.mt"),group.by = "Sample_id")
dev.off()



bharat_immune_no_doublet <- NormalizeData(bharat_immune_no_doublet, normalization.method = "LogNormalize", scale.factor = 10000,verbose = TRUE)
bharat_immune_no_doublet <- FindVariableFeatures(bharat_immune_no_doublet, selection.method = "vst", nfeatures = 2000,verbose = TRUE)
variable.genes=VariableFeatures(bharat_immune_no_doublet,selection.method = "vst",verbose=TRUE)
bharat_immune_no_doublet <- ScaleData(bharat_immune_no_doublet,verbose = TRUE,features = variable.genes)
bharat_immune_no_doublet <- RunPCA(bharat_immune_no_doublet, features = variable.genes,verbose = TRUE,ndims.print = 1:50)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_Cells_Elbowplot_doublet_removed.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

ElbowPlot(bharat_immune_no_doublet)

dev.off()
vars.per.PC=bharat_immune_no_doublet@reductions$pca@stdev*100/sum(bharat_immune_no_doublet@reductions$pca@stdev)
vars.per.PC*100/sum(vars.per.PC)
sum(vars.per.PC[1:40])

bharat_immune_no_doublet <- RunUMAP(bharat_immune_no_doublet, dims=1:40,verbose = TRUE,reduction = "pca")


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_UMAP_doublet_removed.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
x=DimPlot(bharat_immune_no_doublet, group.by = "Doublet_Identification_Call",raster = FALSE,reduction = "umap")
dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_UMAP_doublet_removed_phenotype.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
DimPlot(bharat_immune_no_doublet, group.by = "Phenotype",raster = FALSE,reduction = "umap")
dev.off()

saveRDS(bharat_immune_no_doublet,file='bharat_immune_No_doublet_for_clustering.rds')

bharat_immune_no_doublet=FindNeighbors(bharat_immune_no_doublet,verbose = TRUE,dims = 1:40)

bharat_immune_no_doublet <- FindClusters(bharat_immune_no_doublet, resolution = 0.5,algorithm = 4,method = "igraph",verbose = TRUE) ## Leiden algorithm

table(Idents(bharat_immune_no_doublet))



saveRDS(bharat_immune_no_doublet,file="bharat_immune_cells_leiden_clustered_for_further_analysis.rds")


bharat_immune_no_doublet=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_No_doublet_with_cluster_ids.rds')

table(Idents(bharat_immune_no_doublet))
bharat_immune_no_doublet=BuildClusterTree(bharat_immune_no_doublet,verbose = TRUE)
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_No_doublet_cluster_dendogram.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

PlotClusterTree(bharat_immune_no_doublet)
dev.off()

library(cowplot)
p2 <- DimPlot(bharat_immune_no_doublet, reduction = "umap", label = TRUE, repel = TRUE,label.size = 9)
p3 = DimPlot(bharat_immune_no_doublet, group.by = "Phenotype",raster = FALSE,reduction = "umap")



pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_UMAP_doublet_removed_leiden_clustered.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(p2)

dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_UMAP_doublet_removed_leiden_clustered_possibly_pheno_on_cluster.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(p3)

dev.off()


library(ggplot2)

goi=c("MRC1", "MARCO", "FABP4", "SPP1","CCL2", "CD3E", "CD4", "CD8A", "JCHAIN", "CD79A", "KIT", "CCR7", "CD1C", "FCER1A", "CLEC9A", "CLEC4C","VWF")
heatmap.object=DoHeatmap(bharat_immune_no_doublet, features = goi,group.by = "seurat_clusters",slot = 'data')
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_UMAP_doublet_removed_leiden_clustered_marker_gene_by_cluster_heatmap_blue_white_red.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)

dev.off()

View(bharat_immune_no_doublet@meta.data)
meta_plot_covid=bharat_immune_no_doublet@meta.data
meta_plot_covid$Count_per_cell=1
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_UMAP_doublet_removed_leiden_clustered_per_cluster_COVID_status_count.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

ggplot(meta_plot_covid, aes(fill=Covid_Status, y=Count_per_cell, x=seurat_clusters)) + 
  geom_bar(position="stack",stat="identity") +
  ggtitle("Cluster by COVID19") +
  xlab("Cluster Id")
dev.off()


heatmap_cluste_id <- DimPlot(bharat_immune_no_doublet, reduction = "umap", label = TRUE, repel = TRUE,label.size = 9,raster = FALSE)
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_UMAP_doublet_removed_leiden_clustered_UMAP.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap_cluste_id)
dev.off()

p3 = DimPlot(bharat_immune_no_doublet, group.by = "Phenotype",raster = FALSE,reduction = "umap")

metadata_no_dub=bharat_immune_no_doublet@meta.data
covid=subset(metadata_no_dub,metadata_no_dub$Phenotype=="COVID19/CD31" | metadata_no_dub$Phenotype=="COVID19/Myeloid")
not_covid=subset(metadata_no_dub,metadata_no_dub$Phenotype!="COVID19/CD31" & metadata_no_dub$Phenotype!="COVID19/Myeloid")
covid$Covid_Status="Yes"
not_covid$Covid_Status="No"
metadata_no_dub=rbind(covid,not_covid)
bharat_immune_no_doublet@meta.data=metadata_no_dub

p3 = DimPlot(bharat_immune_no_doublet, group.by = "Covid_Status",raster = FALSE,reduction = "umap")

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_UMAP_doublet_removed_leiden_clustered_covid_status.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(p3)
dev.off()



pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_UMAP_doublet_removed_leiden_clustered_marker_gene_by_cluster_violin_plot.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

VlnPlot(bharat_immune_no_doublet, features = goi,group.by = "seurat_clusters",flip=TRUE,log=TRUE)

dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_UMAP_doublet_removed_leiden_clustered_dendogram.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

PlotClusterTree(bharat_immune_no_doublet)

dev.off()





Idents(bharat_immune_no_doublet)

heatmap.object=DoHeatmap(bharat_immune_no_doublet, features = goi,group.by = "seurat_clusters",slot = 'data',group.colors = c("blue","red"))  #+ scale_fill_gradient(colors = c("blue", "white", "red"))


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_UMAP_doublet_removed_leiden_clustered_marker_gene_by_cluster_heatmap.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)


dev.off()



#### Use the expression pattern of NM to combine clusters


heatmap.object=DoHeatmap(bharat_immune_no_doublet, features = "HBM",group.by = "seurat_clusters",slot = 'data',group.colors = c("blue","red"))  #+ scale_fill_gradient(colors = c("blue", "white", "red"))

VlnPlot(bharat_immune_no_doublet, features = "HBM",group.by = "seurat_clusters",flip=TRUE,log=FALSE)
VlnPlot(bharat_immune_no_doublet, features = "COL1A1",group.by = "seurat_clusters",flip=TRUE,log=FALSE)
VlnPlot(bharat_immune_no_doublet, features = "FBLN1",group.by = "seurat_clusters",flip=TRUE,log=FALSE)
VlnPlot(bharat_immune_no_doublet, features = "MMP2",group.by = "seurat_clusters",flip=TRUE,log=FALSE)
VlnPlot(bharat_immune_no_doublet, features = "HBM",group.by = "seurat_clusters",flip=TRUE,log=FALSE)


VlnPlot(bharat_immune_no_doublet, features ="COL1A1" ,group.by = "seurat_clusters",flip=TRUE,log=FALSE)

VlnPlot(bharat_immune_no_doublet, features = "HBM",group.by = "seurat_clusters",flip=TRUE,log=FALSE,idents = 10)


metat_after_clustering=bharat_immune_no_doublet@meta.data



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



bharat_immune_no_doublet@meta.data=metat_after_clustering
saveRDS(bharat_immune_no_doublet,file="bharat_immune_cells_leiden_clustered_for_further_analysis.rds")





genes_paper=c("ACTA2", "PLIN2", "TWIST1", "LTBP2", "VEGFA", "ITGA8", "ABL2", "MYC", "FAP","PDGFRA", "PDGFRB", "CSF1", "CD34", "FBN1", "FBLN2", "VIT", "NEBL", "HAS2", "CXCL14")


heatmap.object=DoHeatmap(bharat_immune_no_doublet, features = genes_paper,group.by = "seurat_clusters",slot = 'data')
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_heatmap_paper_mentioned_genes.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)

dev.off()


habermann_immune_genes=c("CCL2", "VCAN", "VEGFA", "FCN1", "SELENOP", "SLC40A1","SPP1", "IL1RN", "MMP9", "CHI3L1", "MARCKS", "PLA2G7","PPARG", "MARCO", "MRC1", "INHBA","SLC9B","CKB", "SIGLEC15", "AK5")

heatmap.object=DoHeatmap(bharat_immune_no_doublet, features = habermann_immune_genes,group.by = "seurat_clusters",slot = 'data')
pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/bharat_immune_heatmap_habermann_paper_mentioned_genes.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

plot_grid(heatmap.object)

dev.off()


saveRDS(bharat_immune_no_doublet,file="bharat_immune_cells_leiden_clustered_for_Habermann_integration.rds")

