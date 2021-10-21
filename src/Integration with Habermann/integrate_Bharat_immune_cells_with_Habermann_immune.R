setwd("/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering")


bharat=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_immune_clustered_0.2_resolution_with_cell_type_ids.rds')
immunne_doublet_NM=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Immune_doublets_detected_by_Bharat.csv",header = T,stringsAsFactors = F)
habermann=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/All_Habermann_etal_samples_Seurat_object_after_various_QC_ONLY_USE_THIS.rds')


bharat_clustered_meta=bharat@meta.data
VlnPlot(bharat,features = "percent.mt",group.by = "Sample_id",pt.size = 0)
bharat@meta.data$sample_cell_id=rownames(bharat@meta.data)
length(intersect(bharat_clustered_meta$sample_cell_id,immunne_doublet_NM$Doublet_id)) ## 64953 of the 104719 cells I am using are doublet as per NM analysis

`%nin%` = Negate(`%in%`)
bharat_no_NM_doublets=subset(bharat,subset = sample_cell_id %nin% immunne_doublet_NM$Doublet_id)

View(bharat@meta.data)
#bharat=subset(bharat,subset = Cell_type %in% c("MoAM1/MoAM2","MoAM3","TRAM1/TRAM2/TRAM3"))
habermann=subset(habermann,subset= population %in% c("Immune"))

original.metadata.bharat=bharat@meta.data
original.metadata.habermann=habermann@meta.data
#metadata.habermann.imune=subset(original.metadata.habermann,original.metadata.habermann$population=="Immune")

prop.table(table(original.metadata.bharat$Cell_type))*100


bharat.metadata.use=original.metadata.bharat[,c(2,3,4,7,8,17)]
habermann.metadata.use=original.metadata.habermann[,c(1,2,3,4,5,10)]
colnames(bharat.metadata.use)
colnames(habermann.metadata.use)
bharat.metadata.use$orig.ident=bharat.metadata.use$Sample_id
bharat.metadata.use$Sample_id=NULL
colnames(habermann.metadata.use)[1]="orig.ident"
write.csv(bharat.metadata.use,file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_immune_cells_for_integration_metadata.csv",col.names =T,row.names = T,quote = F)
write.csv(habermann.metadata.use,file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Habermann_immune_cells_for_integration_metadata.csv",col.names =T,row.names = T,quote = F)

#bharat.metadata.use=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_immune_cells_for_integration_metadata.csv",header = T,stringsAsFactors = F)
#habermann.metadata.use=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Habermann_immune_cells_for_integration_metadata.csv",header = T,stringsAsFactors = F)


bharat.metadata.use$Study="Bharat"
habermann.metadata.use$Study="Habermann"

head(bharat.metadata.use)
head(habermann.metadata.use)
bharat.metadata.use=bharat.metadata.use[,c(6,1,2,3,4,5,7)]
colnames(bharat.metadata.use)=colnames(habermann.metadata.use)
#head(bharat@meta.data)
#head(habermann@meta.data)
bharat@meta.data=bharat.metadata.use
habermann@meta.data=habermann.metadata.use


library(Seurat)
immune_markers=c("MRC1", "MARCO", "FABP4", "SPP1", "CCL2", "CD3E", "CD4", "CD8A", 
                 "JCHAIN", "CD79A", "KIT", "CCR7", "CD1C", "FCER1A", "CLEC9A", "CLEC4C",
                 "VWF")


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharar_immune_feature_plot_for_all_immune_marker_genes.pdf",   # The directory you want to save the file in
    width = 30, # The width of the plot in inches
    height = 30)

FeaturePlot(bharat,features = immune_markers,reduction = "umap",keep.scale = "all",label = TRUE,raster = FALSE)
dev.off()
### Bharat cluster annotation
## 6 is B cells
## 9 is T cells
## 13 is DC cells
## 3 is MoAM1 cells
## 11 is plasma cells
## not sure of TRAM2
## 

FeaturePlot(bharat,features = c("KIT","FCER1A"),reduction = "umap",keep.scale = "all",label = TRUE,raster = FALSE)
FeaturePlot(bharat,features = c("CD4"),reduction = "umap",keep.scale = "all",label = TRUE,raster = FALSE)
FeaturePlot(bharat,features = c("CD1C","FCER1A"),reduction = "umap",keep.scale = "all",label = TRUE,raster = FALSE)

bharat@meta.data=original.metadata.bharat
VlnPlot(bharat,features = c("KIT","FCER1A"),group.by = "RNA_snn_res.0.2",pt.size = 0.5)
covid=subset(original.metadata.bharat,Phenotype=="COVID19/CD31" | Phenotype=="COVID19/Myeloid" | Phenotype=="Transplant Recipient/CD31" | Phenotype=="Transplant Recipient/Myeloid")
not_covid=subset(original.metadata.bharat,Phenotype!="COVID19/CD31" & Phenotype!="COVID19/Myeloid" & Phenotype!="Transplant Recipient/CD31" & Phenotype!="Transplant Recipient/Myeloid")
covid$Covid_Status="COVID19"
not_covid$Covid_Status="Control"
original.metadata.bharat=rbind(covid,not_covid)


original.metadata.bharat$Cell_Count=1
ggplot(original.metadata.bharat, aes(fill=Covid_Status, y=Cell_Count, x=`RNA_snn_res.0.2`)) + 
  geom_bar(position="stack",stat="identity") +
  ggtitle("Cluster by Phenotype") +
  xlab("Cluster Id")

## 6= Bcells, 2,4,8 are TRAMs probably, 9= T cells, 3=MoAM3, 
original.metadata.bharat$Cell_type_new[original.metadata.bharat$RNA_snn_res.0.2 %in% c(3)]="MoAM3/MoAM2" ## most covid+ve
original.metadata.bharat$Cell_type_new[original.metadata.bharat$RNA_snn_res.0.2 %in% c(6)]="B cells" ## most covid+ve
original.metadata.bharat$Cell_type_new[original.metadata.bharat$RNA_snn_res.0.2 %in% c(1)]="MoAM1" ## most covid+ve
original.metadata.bharat$Cell_type_new[original.metadata.bharat$RNA_snn_res.0.2 %in% c(4)]="TRAM2" ## most covid+ve
original.metadata.bharat$Cell_type_new[original.metadata.bharat$RNA_snn_res.0.2 %in% c(2)]="TRAM1" ## most covid+ve
original.metadata.bharat$Cell_type_new[original.metadata.bharat$RNA_snn_res.0.2 %in% c(8)]="TRAM3" ## most covid+ve
original.metadata.bharat$Cell_type_new[original.metadata.bharat$RNA_snn_res.0.2 %in% c(11)]="Plasma Cells" ## most covid+ve
original.metadata.bharat$Cell_type_new[original.metadata.bharat$RNA_snn_res.0.2 %in% c(9)]="T Cells" ## most covid+ve
original.metadata.bharat$Cell_type_new[original.metadata.bharat$RNA_snn_res.0.2 %in% c(13)]="Mast Cells" ## most covid+ve
colnames(original.metadata.bharat)[17]="Cell_type_not_to_use"
table(Not_to_use=original.metadata.bharat$Cell_type_not_to_use,To_use=original.metadata.bharat$Cell_type_new)
bharat@meta.data=original.metadata.bharat

moam_tram_bharat = subset(bharat,subset = Cell_type_new %in% c("MoAM1","MoAM3/MoAM2","TRAM1","TRAM2","TRAM3"))
habermannn_genes=c("CCL2", "VCAN", "VEGFA", "FCN1", "SELENOP", "SLC40A1","SPP1", "IL1RN", "MMP9", "CHI3L1", "MARCKS", "PLA2G7","PPARG", "MARCO", "MRC1", "INHBA","SLC9B2", "CKB", "SIGLEC15", "AK5")

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Testing_MoAM_TRAm_Habermann_genes.pdf",   # The directory you want to save the file in
    width = 30, # The width of the plot in inches
    height = 30)

DoHeatmap(moam_tram_bharat, features =habermannn_genes,group.by = "Cell_type_new",group.bar = TRUE,draw.lines = FALSE,slot='data')
dev.off()
### MoAM looks good, TRAM3 to be excluded looks weird

#View(bharat_use@meta.data)
#View(habermann@meta.data)


#habermann_immune = subset(habermann,subset = population=="Immune")
#head(bharat_use@meta.data)
#head(habermann_immune@meta.data)
#table(bharat_use@meta.data$Doublet_Identification_Call)


#bharat_use@meta.data$RNA_snn_res.0.4=NULL
#bharat_use@meta.data$RNA_snn_res.0.8=NULL
#bharat_use@meta.data$RNA_snn_res.1.2=NULL
#bharat_use@meta.data$Doublet_Identification_Call=NULL
#bharat_use@meta.data$population="Immune"
#bharat_use@meta.data$orig.ident=bharat_use@meta.data$Sample_id
#colnames(bharat_use@meta.data)[6]="Diagnosis"
#colnames(bharat_use@meta.data)[8]="celltype"
#colnames(bharat_use@meta.data)[5]="Sample_Name"




#rm(bharat,habermann)
#### UNINTEGRATED ANALYSIS



all_immune_merged=merge(bharat, y = c(habermann), add.cell.ids = c("Bharat", "Habermann"), project = "immune")

View(all_immune_merged@meta.data)



#all_immune_merged@meta.data$nCount_SCT=NULL
#all_immune_merged@meta.data$nFeature_SCT=NULL



all_immune_merged <- NormalizeData(all_immune_merged, normalization.method = "LogNormalize", scale.factor = 10000,verbose = TRUE)
all_immune_merged <- FindVariableFeatures(all_immune_merged, selection.method = "vst", nfeatures = 3000,verbose = TRUE)
top10 <- head(VariableFeatures(all_immune_merged), 10)
variable.genes=VariableFeatures(all_immune_merged)
plot1 <- VariableFeaturePlot(all_immune_merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(all_immune_merged)
all_immune_merged <- ScaleData(all_immune_merged,verbose = TRUE)
all_immune_merged <- RunPCA(all_immune_merged, features = variable.genes,verbose = TRUE)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Habermann_Bharat_immune_not_integrated_Elbowplot.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)

ElbowPlot(all_immune_merged)

dev.off()

DimPlot(all_immune_merged, reduction = "pca",raster = FALSE,group.by = "Study")

x=all_immune_merged@reductions$pca@stdev*100/sum(all_immune_merged@reductions$pca@stdev)
cumsum(x)

all_immune_merged <- RunUMAP(all_immune_merged, dims=1:40,verbose = TRUE,reduction = "pca")


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Habermann_Bharat_immune_not_integrated_UMAP.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

DimPlot(all_immune_merged, group.by = "Study",label = TRUE,raster = FALSE,reduction = "umap")

dev.off()




#### INTEGRATED ANALYSIS

table(habermann@meta.data$Cell_type)
habermann_for_integ=subset(habermann,subset = Cell_type %in% c("Macrophages","Monocytes","Proliferating Macrophages"))


habermann_reduced_meta=habermann@meta.data[,c(1,2,3,7)]
bharat_reduced_meta=bharat_no_NM_doublets@meta.data[,c("nCount_RNA","nFeature_RNA","Sample_id","sample_cell_id")]


samples.seurat.list=list(bharat,habermann)

samples.seurat.list <- lapply(X = samples.seurat.list, FUN = function(x) {
  x <- SCTransform(x,verbose = TRUE)
})



# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = samples.seurat.list,verbose = TRUE,nfeatures =3000) ## trying with jkust 2000 features to see how well it works
write.table(features,file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Genes_being_used_for_immune_integrative_analysis.txt",col.names = F,row.names = F,sep="\t",quote = F)
## Perform PCA before integration using the 2000 variable features
samples.seurat.list[[1]]=RunPCA(samples.seurat.list[[1]],verbose = TRUE,features=features)
samples.seurat.list[[2]]=RunPCA(samples.seurat.list[[2]],verbose = TRUE,features=features)


### Running PrepSCTIntegration for check before integration
samples.seurat.list=PrepSCTIntegration(samples.seurat.list,verbose = TRUE,anchor.features=features)
### IF FEATURES DOES NOT LOOK GOOD THEN USE PCs


immune.anchors <- FindIntegrationAnchors(object.list = samples.seurat.list,normalization.method = "SCT",reduction = "rpca",verbose = TRUE,anchor.features = features)

immune.integrated <- IntegrateData(anchorset = immune.anchors,verbose = TRUE,normalization.method = "SCT")



saveRDS(immune.integrated,file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_macrophages_monocytes_integrated.rds")





saveRDS(immune.integrated,file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_macrophages_monocytes_integrated_for_clustering.rds")

integrated_not_clustered=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_macrophages_monocytes_integrated_for_clustering.rds')
View(integrated_not_clustered@meta.data)
integrated_not_clustered@meta.data$Cell_id_new = rownames(integrated_not_clustered@meta.data)



DefaultAssay(object = immune.integrated) <- "integrated"
immune.integrated <- ScaleData(object = immune.integrated, verbose = TRUE)
immune.integrated <- RunPCA(object = immune.integrated, npcs = 50, verbose = TRUE)

deviation=immune.integrated@reductions$pca@stdev*100/(sum(immune.integrated@reductions$pca@stdev))
ElbowPlot(immune.integrated,ndims=50,reduction = "pca")
immune.integrated <- RunUMAP(immune.integrated, dims=1:45,verbose = TRUE,reduction = "pca")




meta.integ=immune.integrated@meta.data



write.csv(meta.integ,file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Immune_integrated_metadata_to_Edit.csv",col.names = T,row.names = T,quote = F)

#meta_bharat=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Immune_metadata_Bharat_batch_ccorrected.csv",row.names = 1,header = TRUE,stringsAsFactors = F)
#meta_habermann=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Immune_metadata_Habermann_batch_ccorrected.csv",row.names = 1,header = TRUE,stringsAsFactors = F)

#meta_bharat$Study="Bharat"
#meta_habermann$Study="Habermann"

#combined_corrected_meta=rbind(meta_bharat,meta_habermann)

#length(intersect(rownames(combined_corrected_meta),colnames(immune.integrated)))
#immune.integrated@meta.data=combined_corrected_meta

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Habermann_Bharat_immune_integrated_UMAP.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

DimPlot(immune.integrated, group.by = "Study",label = TRUE,raster = FALSE,reduction = "umap")

dev.off()



immune.integrated=FindNeighbors(immune.integrated,verbose = TRUE,dims = 1:45)

immune.integrated <- FindClusters(immune.integrated, resolution = 0.2,algorithm = 4,method = "igraph",verbose = TRUE) ## Leiden algorithm


saveRDS(immune.integrated,file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_macrophages_monocytes_integrated_clustered.rds")


immune.integrated=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_macrophages_monocytes_integrated_clustered.rds')

Idents(immune.integrated)
View(immune.integrated@meta.data)

immune.integrated=BuildClusterTree(immune.integrated,verbose = TRUE)
PlotClusterTree(immune.integrated) ## seems like 9,11,5,13,3,1,2,4 are distinct enough
table(Idents(immune.integrated))
table(immune.integrated@meta.data$celltype)

immune.integrated@meta.data$Cell_id_new = rownames(immune.integrated@meta.data)
immune.integrated[["percent.mt"]] <- PercentageFeatureSet(immune.integrated, pattern = "^MT-", assay = 'RNA')


integrated_clustered_meta=immune.integrated@meta.data
VlnPlot(immune.integrated,features = "nFeature_RNA",group.by = "seurat_clusters",pt.size = 0)
VlnPlot(immune.integrated,features = "percent.mt",group.by = "seurat_clusters",pt.size = 0) ## 10 seems to have high mitochondria
### dropping 10 as too high mitochondria and 3,7 as feature count too low

### Finally keeping 9,11,5,13,1,2,4,3
subsetted_data=subset(immune.integrated, idents = c(1,2,3,4,5,9,11,13))

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_immune_integrated_select_clusters_paper_genes_heatmap.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

DoHeatmap(subsetted_data, features =goi,group.by = "ident",group.bar = TRUE,draw.lines = FALSE,slot='data')
dev.off()

##Just based on counts,7 is  MoAM2, 3 is MoAM3, 4 is Monocytes, 1 is MoAM1, 5,6,8 and 12 is AM1 and --- is AM2


VlnPlot(immune.integrated,features = c("FABP4","INHBA","PPARG"),group.by = "seurat_clusters",pt.size = 0,sort="increasing") ## mostly 5,6,8,12,14 AM2/TRAM
VlnPlot(immune.integrated,features = c("FCN1","CD300E","APOBEC3A"),group.by = "seurat_clusters",pt.size = 0) ## mostly 1,2,4,14 Monocytes
VlnPlot(immune.integrated,features = c("MRC1","SPP1","PLA2G7","CCL2","AK5","SIGLEC15","CKB","SLC9B2","CD163","SLC40A1","MERTK","PLTP","ABCA1","MRC1","CLEC10A","F13A1"),group.by = "seurat_clusters",pt.size = 0) ## mostly 1,2,4,14 Monocytes
VlnPlot(immune.integrated,features = c("AK5","SIGLEC15","CKB","SLC9B2"),group.by = "seurat_clusters",pt.size = 0) ## mostly 1,2,4,14 MoAM3
VlnPlot(immune.integrated,features = c("AK5"),group.by = "seurat_clusters",pt.size = 1,sort="increasing") ## mostly 1,2,4,14 MoAM3
VlnPlot(immune.integrated,features = c("SIGLEC15"),group.by = "seurat_clusters",pt.size = 1,sort="increasing") ## mostly 1,2,14 MoAM3
VlnPlot(immune.integrated,features = c("PLA2G7","SPP1","CCL2"),group.by = "seurat_clusters",pt.size = 1,sort="decreasing") ## mostly 

VlnPlot(immune.integrated,features = c("FABP4","INHBA","PPARG"),group.by = "seurat_clusters",pt.size = 0,sort="increasing") ## mostly 5,6,8,12,14 AM2/TRAM based on Bharat cluster 1,5,7,8
VlnPlot(immune.integrated,features = c("SPP1","PLA2G7","CCL2","AK5","SIGLEC15","CKB","SLC9B2"),group.by = "seurat_clusters",pt.size = 0,sort="increasing") ##All MoAM genes based on Bharat CLuter 0
VlnPlot(immune.integrated,features = c("CD163","SLC40A1","MERTK","PLTP","ABCA1"),group.by = "seurat_clusters",pt.size = 0,sort="increasing") ## mostly 15,1,7 All case 1 MoAM based onn Bhafrat cluster 3
VlnPlot(immune.integrated,features = c("CLEC10A","F13A1"),group.by = "seurat_clusters",pt.size = 0.1,sort="increasing") ## mostly 15,1,4,2 All moAM based on Bharat cluster 6 
VlnPlot(immune.integrated,features = c("FCN1","CD300E","APOBEC3A","FABP4","INHBA","PPARG"),group.by = "seurat_clusters",pt.size = 0.1,sort="decreasing") ## mostly 1,2,4,14 Monocytes

#### "FCN1","CD300E","APOBEC3A" dominate in 1,2,4 and are monocytes
### "FABP4","INHBA","PPARG" dominate in 5,6,8 and are AM2
### 


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_immune_integrated_select_clusters_paper_genes_sorted_violin_plot.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

VlnPlot(immune.integrated,features = goi,group.by = "seurat_clusters",pt.size = 0,sort="increasing") ## mostly 
dev.off()


## MoM4 would be only SPP1, PLA2G7 and CCL2 strongest. mostly 2
## 


DimPlot(subsetted_data,label = TRUE,raster = FALSE,reduction = "umap") 




immune.integrated.no.c10=subset(immune.integrated, subset = seurat_clusters!=10)
integrated_clustered_meta=immune.integrated.no.c10@meta.data


combined_pheno=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_Combined_phenotype.csv",header = T,stringsAsFactors = F)

combined_pheno=combined_pheno[,c(1,3)]
colnames(combined_pheno)[2]="Phenotype"
library(dplyr)

integrated_clustered_meta_pheno=left_join(integrated_clustered_meta,combined_pheno,by=c("orig.ident"="Sample"))
rownames(integrated_clustered_meta_pheno)=integrated_clustered_meta_pheno$Cell_id_new
integrated_clustered_meta_pheno$Cell_Count=1

library(ggplot2)
ggplot(integrated_clustered_meta_pheno, aes(fill=celltype, y=Cell_Count, x=seurat_clusters)) + 
  geom_bar(position="stack",stat="identity") +
  ggtitle("Cluster by CellType") +
  xlab("Cluster Id")

ggplot(integrated_clustered_meta_pheno, aes(fill=Phenotype, y=Cell_Count, x=seurat_clusters)) + 
  geom_bar(position="stack",stat="identity") +
  ggtitle("Cluster by Phenotype") +
  xlab("Cluster Id")

table(integrated_clustered_meta_pheno$Phenotype,integrated_clustered_meta_pheno$Study)


goi=c("MRC1", "FABP4","INHBA","PPARG","SPP1","PLA2G7","CCL2","AK5","SIGLEC15","CKB","SLC9B2","CD163","SLC40A1","MERTK","PLTP","ABCA1","MRC1","CLEC10A","F13A1","FCN1","APOBEC3A","CD300E")



pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_immune_integrated_paper_genes_heatmap.pdf",   # The directory you want to save the file in
    width = 30, # The width of the plot in inches
    height = 30)

DoHeatmap(immune.integrated, features =goi,group.by = "seurat_clusters",group.bar = TRUE,slot = 'data',draw.lines = FALSE)
dev.off()


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_immune_integrated_paper_genes_violin_plot.pdf",   # The directory you want to save the file in
    width = 30, # The width of the plot in inches
    height = 30)

VlnPlot(immune.integrated,features = goi,group.by = "seurat_clusters",pt.size = 0)
dev.off()

immune.integrated@meta.data=integrated_clustered_meta_pheno

DimPlot(immune.integrated,raster = FALSE,reduction = "umap",group.by = 'Phenotype',label = FALSE,label.size = NA)


DimPlot(immune.integrated,raster = FALSE,reduction = "umap",group.by = 'seurat_clusters',label = TRUE,label.size = NA)

sample_id_info=read.csv(file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_Combined_phenotype.csv",header = T,stringsAsFactors = F)
sample_id_info=sample_id_info[,c("Sample","Id_by_pheno")]

integrated_clustered_meta_pheno_with_sample_id=left_join(integrated_clustered_meta_pheno,sample_id_info,by=c("orig.ident"="Sample"))

rownames(integrated_clustered_meta_pheno_with_sample_id)=integrated_clustered_meta_pheno_with_sample_id$Cell_id_new

immune.integrated@meta.data=integrated_clustered_meta_pheno_with_sample_id


DefaultAssay(immune.integrated)="integrated"

immune.integrated=NormalizeData(immune.integrated,verbose = TRUE,normalization.method = "LogNormalize")
avg_expn=AverageExpression(immune.integrated,features =goi,group.by = c("seurat_clusters","Id_by_pheno"),slot="data")


avg_expn_sct=avg_expn$SCT

library(gplots)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)


heatmap.2(avg_expn_sct,scale = 'row',Rowv = "NA",Colv = "Rowv",dendrogram = "none",col = my_palette,lhei=c(2,5),lwid=c(2,5))
### still need to plot the mean average expn by cluster for each phenotype


library(Seurat)



### CLuster 3 mostly has control by phenotype but by celltype seems to be MoAM2 which is weird
### checking MoAM2 specific genes in all the clusters MRC1, MARCO,CD4


DefaultAssay(immune.integrated)="integrated"
immune.integrated=ScaleData(immune.integrated,verbose = TRUE,normalization.method = "LogNormalize")
DoHeatmap(immune.integrated, features =c("MRC1","MARCO","CD4"),group.by = "seurat_clusters",group.bar = TRUE,draw.lines = FALSE)


cluster1=subset(immune.integrated, subset = seurat_clusters==1)



metadata=immune.integrated@meta.data


write.csv(metadata,file="/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Immune_integrated_macrophahes_monocytes_only_metadata_info.csv",col.names = T,row.names = T,quote = F)


### taking 1,2,4,5,6,8

table(Idents(immune.integrated))














































### NOT USEFUL AS CELL COUNT TOO BECOMES TOO LOW COMPARED TO NM ANALYSIS

immune.integrated_subset = subset(immune.integrated , subset = seurat_clusters %in% c(1,2,4,5,6,8))

DefaultAssay(immune.integrated_subset)="integrated"
immune.integrated_subset <- ScaleData(object = immune.integrated_subset, verbose = TRUE)
immune.integrated_subset <- RunPCA(object = immune.integrated_subset, npcs = 50, verbose = TRUE)
x=immune.integrated_subset@reductions$pca@stdev*100/sum(immune.integrated_subset@reductions$pca@stdev)
cumsum(x)
PCAPlot(object = immune.integrated_subset, dims=c(1,2))

ElbowPlot(immune.integrated_subset,ndims = 50)
immune.integrated_subset=FindNeighbors(immune.integrated_subset,verbose = TRUE,dims = 1:40)

immune.integrated_subset <- FindClusters(immune.integrated_subset, resolution =0.5,algorithm = 4,method = "igraph",verbose = TRUE) ## Leiden algorithm
table(Idents(immune.integrated_subset))

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_immune_testing_second_clustering_scaled_data.pdf",   # The directory you want to save the file in
    width = 30, # The width of the plot in inches
    height = 30)

DoHeatmap(immune.integrated_subset, features =goi,group.by = "ident",group.bar = TRUE,draw.lines = FALSE,slot='scale.data')
dev.off()


immune.integrated_subset_meta=immune.integrated_subset@meta.data

integrated_clustered_meta_pheno=left_join(immune.integrated_subset_meta,combined_pheno,by=c("orig.ident"="Sample"))
rownames(integrated_clustered_meta_pheno)=integrated_clustered_meta_pheno$Cell_id_new
integrated_clustered_meta_pheno$Cell_Count=1

library(ggplot2)
ggplot(integrated_clustered_meta_pheno, aes(fill=celltype, y=Cell_Count, x=seurat_clusters)) + 
  geom_bar(position="stack",stat="identity") +
  ggtitle("Cluster by CellType") +
  xlab("Cluster Id")

ggplot(integrated_clustered_meta_pheno, aes(fill=Phenotype, y=Cell_Count, x=seurat_clusters)) + 
  geom_bar(position="stack",stat="identity") +
  ggtitle("Cluster by Phenotype") +
  xlab("Cluster Id")

bharat_meta=subset(integrated_clustered_meta_pheno,integrated_clustered_meta_pheno$Study=="Bharat")
habermann_meta=subset(integrated_clustered_meta_pheno,integrated_clustered_meta_pheno$Study=="Habermann")
rownames(bharat_meta)=paste(bharat_meta$orig.ident,bharat_meta$Cell_id,sep="_")
rownames(habermann_meta)=habermann_meta$Cell_id

integrated_clustered_meta_pheno_new=rbind(bharat_meta,habermann_meta)
immune.integrated_subset@meta.data=integrated_clustered_meta_pheno_new

rownames(integrated_clustered_meta_pheno)=paste(integrated_clustered_meta_pheno$orig.ident,integrated_clustered_meta_pheno$Cell_id)

DimPlot(immune.integrated_subset,raster = FALSE,reduction = "umap",group.by = 'Phenotype',label = FALSE,label.size = NA)

DimPlot(immune.integrated_subset,raster = FALSE,reduction = "umap",group.by = 'seurat_clusters',label = FALSE,label.size = NA)


immune.integrated_subset=BuildClusterTree(immune.integrated_subset,verbose = TRUE)
PlotClusterTree(immune.integrated_subset)

first8_clusters=subset(immune.integrated_subset,subset = seurat_clusters %in% c(1:8))

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_Habermann_immune_testing_second_clustering_violin_plot.pdf",   # The directory you want to save the file in
    width = 30, # The width of the plot in inches
    height = 30)

VlnPlot(immune.integrated_subset,features = goi,group.by = "seurat_clusters",pt.size = 0)
dev.off()

## 9,11,16 is monocyte
### 4,6,7 is AM2
### 2, 5 is MoM3
### 



SC84=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC86_seurat_min_3_cells_min_200_genes.rds')
SC87=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC87_seurat_min_3_cells_min_200_genes.rds')
SC109=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC109_seurat_min_3_cells_min_200_genes.rds')
SC142=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC142_seurat_min_3_cells_min_200_genes.rds')
SC144=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC144_seurat_min_3_cells_min_200_genes.rds')
SC146=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC146_seurat_min_3_cells_min_200_genes.rds')
SC148=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC148_seurat_min_3_cells_min_200_genes.rds')
SC149=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC149_seurat_min_3_cells_min_200_genes.rds')
SC150=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/SC150_seurat_min_3_cells_min_200_genes.rds')

SC84@meta.data$Sample="SC84"
SC87@meta.data$Sample="SC87"
SC109@meta.data$Sample="SC109"
SC142@meta.data$Sample="SC142"
SC144@meta.data$Sample="SC144"
SC146@meta.data$Sample="SC146"
SC148@meta.data$Sample="SC148"
SC149@meta.data$Sample="SC149"
SC150@meta.data$Sample="SC150"

SC84@meta.data$composite_id=paste(SC84@meta.data$Sample,rownames(SC84@meta.data),sep="_")
SC87@meta.data$composite_id=paste(SC87@meta.data$Sample,rownames(SC87@meta.data),sep="_")
SC109@meta.data$composite_id=paste(SC109@meta.data$Sample,rownames(SC109@meta.data),sep="_")
SC142@meta.data$composite_id=paste(SC142@meta.data$Sample,rownames(SC142@meta.data),sep="_")
SC144@meta.data$composite_id=paste(SC144@meta.data$Sample,rownames(SC144@meta.data),sep="_")
SC146@meta.data$composite_id=paste(SC146@meta.data$Sample,rownames(SC146@meta.data),sep="_")
SC148@meta.data$composite_id=paste(SC148@meta.data$Sample,rownames(SC148@meta.data),sep="_")
SC149@meta.data$composite_id=paste(SC149@meta.data$Sample,rownames(SC149@meta.data),sep="_")
SC150@meta.data$composite_id=paste(SC150@meta.data$Sample,rownames(SC150@meta.data),sep="_")


View(SC84@meta.data)

all_merged=merge(SC84,y=c(SC87,SC109,SC142,SC144,SC146,SC148,SC149,SC150))

View(all_merged@meta.data)


all_merged[["percent.mt"]] <- PercentageFeatureSet(all_merged, pattern = "^MT-")


VlnPlot(all_merged,features = "percent.mt",group.by = "Sample",pt.size = 0)+ggtitle("Mitochondria percentage by sample before filtering")
VlnPlot(bharat,features = "percent.mt",group.by = "Sample_id",pt.size = 0)+ggtitle("Mitochondria percentage by sample after filtering")





