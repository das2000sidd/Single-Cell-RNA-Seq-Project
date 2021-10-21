setwd("/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering")


bharat=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/Bharat_all_cells_with_singlet_doublet_status_with_leiden_cluster_ids.rds')
habermann=readRDS('/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/All_Habermann_etal_samples_Seurat_object_after_various_QC_ONLY_USE_THIS.rds')


original.metadata.bharat=bharat@meta.data
original.metadata.habermann=habermann@meta.data
metadata.habermann.imune=subset(original.metadata.habermann,original.metadata.habermann$population=="Epithelial")


bharat.metadata.use=original.metadata.bharat[,c(1,2,3,4,5)]
habermann.metadata.use=original.metadata.habermann[,c(1,2,3,4,9)]
colnames(habermann.metadata.use)[1]="orig.ident"


bharat.metadata.use$Study="Bharat"
habermann.metadata.use$Study="Habermann"

bharat.metadata.use$orig.ident=bharat.metadata.use$Sample_Id
bharat.metadata.use$Sample_Id=NULL
bharat.metadata.use$cell_oid=rownames(bharat.metadata.use)
colnames(bharat.metadata.use)[4]="population"

#head(bharat.metadata.use)
#head(habermann.metadata.use)


bharat@meta.data=bharat.metadata.use
habermann@meta.data=habermann.metadata.use


#head(bharat@meta.data)
#head(habermann@meta.data)


table(bharat@meta.data$population)
table(habermann@meta.data$population)

library(Seurat)

View(bharat@meta.data)
View(habermann@meta.data)


bharat_epithelial = subset(bharat,subset = population=="Epithelial")
habermann_epithelial = subset(habermann,subset = population=="Epithelial")

table(habermann_epithelial@meta.data$celltype)

#rm(bharat,habermann)
#### UNINTEGRATED ANALYSIS

all_epithelial_merged=merge(bharat_epithelial, y = c(habermann_epithelial), add.cell.ids = c("Bharat", "Habermann"), project = "Epithelial")

View(all_epithelial_merged@meta.data)



all_epithelial_merged@meta.data$nCount_SCT=NULL
all_epithelial_merged@meta.data$nFeature_SCT=NULL



all_epithelial_merged <- NormalizeData(all_epithelial_merged, normalization.method = "LogNormalize", scale.factor = 10000,verbose = TRUE)
all_epithelial_merged <- FindVariableFeatures(all_epithelial_merged, selection.method = "vst", nfeatures = 2000,verbose = TRUE)
top10 <- head(VariableFeatures(all_epithelial_merged), 10)
variable.genes=VariableFeatures(all_epithelial_merged)
plot1 <- VariableFeaturePlot(all_epithelial_merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(all_epithelial_merged)
all_epithelial_merged <- ScaleData(all_epithelial_merged,verbose = TRUE)
all_epithelial_merged <- RunPCA(all_epithelial_merged, features = variable.genes,verbose = TRUE)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Habermann_Bharat_epithelial_not_integrated_Elbowplot.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)

ElbowPlot(all_epithelial_merged)

dev.off()

DimPlot(all_epithelial_merged, reduction = "pca",raster = FALSE,group.by = "Study")
all_epithelial_merged <- RunUMAP(all_epithelial_merged, dims=1:20,verbose = TRUE,reduction = "pca")


pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Habermann_Bharat_epithelial_not_integrated_UMAP.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)

DimPlot(all_epithelial_merged, group.by = "Study",label = TRUE,raster = FALSE,reduction = "umap")

dev.off()




#### INTEGRATED ANALYSIS
samples.seurat.list=list(bharat_epithelial,habermann_epithelial)

samples.seurat.list <- lapply(X = samples.seurat.list, FUN = function(x) {
  x <- SCTransform(x,verbose = TRUE)
})



# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = samples.seurat.list,verbose = TRUE,nfeatures =1000) ## trying with jkust 2000 features to see how well it works

## Perform PCA before integration using the 2000 variable features
samples.seurat.list[[1]]=RunPCA(samples.seurat.list[[1]],verbose = TRUE,features=features)
samples.seurat.list[[2]]=RunPCA(samples.seurat.list[[2]],verbose = TRUE,features=features)



### Running PrepSCTIntegration for check before integration
samples.seurat.list=PrepSCTIntegration(samples.seurat.list,verbose = TRUE,anchor.features=features)



epithelial.anchors <- FindIntegrationAnchors(object.list = samples.seurat.list,normalization.method = "SCT",reduction = "rpca",verbose = TRUE,anchor.features = features)

epithelial.integrated <- IntegrateData(anchorset = epithelial.anchors,verbose = TRUE,normalization.method = "SCT")







# switch to integrated assay. The variable features of this assay are automatically set during
DefaultAssay(object = epithelial.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
epithelial.integrated <- ScaleData(object = epithelial.integrated, verbose = TRUE)
epithelial.integrated <- RunPCA(object = epithelial.integrated, npcs = 30, verbose = TRUE)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Habermann_Bharat_epithelial_integrated_Elbowplot.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)

ElbowPlot(epithelial.integrated)

dev.off()

epithelial.integrated <- RunUMAP(object = epithelial.integrated, reduction = "pca", dims = 1:20,verbose = TRUE)
View(epithelial.integrated@meta.data)

pdf(file = "/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Habermann_Bharat_epithelial_integrated_UMAP_By_Study_type.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
DimPlot(epithelial.integrated, group.by = "Study",raster = FALSE)
dev.off()

samples.seurat.list[[2]]@meta.data=metadata.habermann.imune
colnames(metadata.habermann.imune)[1]="orig.ident"
samples.seurat.list[[2]]@meta.data=metadata.habermann.imune

#colnames(samples.seurat.list[[2]]@meta.data)[1]="orig.ident"
transfer.achors <- FindTransferAnchors(reference =  epithelial.integrated,query = samples.seurat.list[[1]], dims = 1:20,verbose = TRUE,reduction = "rpca",normalization.method = 'SCT')
predictions <- TransferData(anchorset = transfer.achors, reference = samples.seurat.list[[2]],query=samples.seurat.list[[1]],dims = 1:20,refdata = "celltype",verbose = TRUE)


table(predictions@meta.data$predicted.id)


saveRDS(predictions,'/projects/b1025/sdi0596/Covid19_Single_Cell_RNASeq_project/all_h5ad_after_filtering/New_Analysis_Using_Bharat_all_cells/Bharat_epithelial_cells_predicted_cell_type_USING_ONLY_Habermann_epithelial_cell_type.rds')



