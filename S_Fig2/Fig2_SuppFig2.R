library(Seurat) 
library(ggplot2) 
library(dplyr)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(grid)
library(fishplot)
library("randomcoloR")
library("xlsx")

#### 1. Load the reads and QC ####
data <- Read10X(data.dir = "~/302_GEX/filtered_feature_bc_matrix/")
Seurat.obj <- CreateSeuratObject(counts = data, project = "")
Seurat.obj
#Add metadata
mito.genes <- grep(pattern = "^MT-", x = rownames(x = Seurat.obj), value = TRUE)
percent.counts <- Matrix::colSums(Seurat.obj@assays$RNA@data[mito.genes,])
percent.mito <- Matrix::colSums(Seurat.obj@assays$RNA@data[mito.genes,])/Matrix::colSums(Seurat.obj@assays$RNA@data)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = percent.mito, col.name = "percent.mito")
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = percent.counts, col.name = "percent.counts")
Seurat.obj$log10_nCount_RNA <- log10(Seurat.obj$nCount_RNA)
Seurat.obj$Sample_id <- "Seurat.obj"

upper_bound <- mean(Seurat.obj$log10_nCount_RNA) + 2*sd(Seurat.obj$log10_nCount_RNA)
lower_bound <- mean(Seurat.obj$log10_nCount_RNA) - 2*sd(Seurat.obj$log10_nCount_RNA)
p1 <- FeatureScatter(object = Seurat.obj, feature1 = "log10_nCount_RNA", feature2 = "percent.mito")
p2 <- FeatureScatter(object = Seurat.obj, feature1 = "log10_nCount_RNA", feature2 = "nFeature_RNA")
p1 + geom_vline(xintercept=lower_bound) + geom_vline(xintercept=3) + geom_vline(xintercept=4.6) + geom_hline(yintercept=0.2)  + geom_hline(yintercept=0.1)
p2 + geom_vline(xintercept=lower_bound) + geom_hline(yintercept=200) + geom_hline(yintercept=5000)

# Fiilter out cells with very few reads or too much (possible doublets)
Seurat.obj <- subset(x = Seurat.obj, subset = nCount_RNA > 1000 & nCount_RNA < 40000 & nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 0.1)
Seurat.obj

Seurat.obj$nGenes <- Seurat.obj$nFeature_RNA
Seurat.obj$nUMI <- Seurat.obj$nCount_RNA
pdf(file="~/SuppFig2_a.pdf",  width = 6, height = 5)
SuppFig2_a <- VlnPlot(object = Seurat.obj, features = c("nGenes", "nUMI", "percent.mito"), ncol = 3)
SuppFig2_a
dev.off()

#### 2. Normalize the data ####
Seurat.obj <- NormalizeData(object = Seurat.obj, normalization.method = "LogNormalize", scale.factor = 1e4)
Seurat.obj <- FindVariableFeatures(Seurat.obj, selection.method = "vst", nfeatures = 2000)
Seurat.obj <- ScaleData(Seurat.obj,vars.to.regress = c("nCount_RNA","percent.mito"), features = VariableFeatures(Seurat.obj))

#### 3. PCA and markers ####
Seurat.obj <- RunPCA(Seurat.obj, features = VariableFeatures(Seurat.obj), npcs = 40)
ElbowPlot(Seurat.obj,ndims = 40)
#Use first 9 PC
n_PC = 9
Seurat.obj <- FindNeighbors(Seurat.obj, dims = 1:n_PC)
Seurat.obj <- FindClusters(Seurat.obj, resolution = 0.5)

n_cells = ncol(Seurat.obj)
grob <- grobTree(textGrob(paste0("n = ",n_cells), x=0.05,  y=0.95, hjust=0,gp=gpar(fontsize=22)))

#PCA
DimPlot(Seurat.obj, reduction = "pca")

#UMAP
Seurat.obj <- RunUMAP(Seurat.obj, dims = 1:n_PC)

markers_all_cells <- FindAllMarkers(Seurat.obj, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(markers_all_cells,file="~/markers.clusters.csv")

#### 4. Rename the clusters ####
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colors = getPalette(12)
new.cluster.ids <- c("1 Cyto Pre-exh", "2 Cyto Exh", "3 CD4 Naive", "4 CD8 Naive", "5 CD8 Cyto", "6", "7 CD4 Eff/Mem", "8 CD8 Pre-exh/Exh", "9 Degraded")
names(new.cluster.ids) <- levels(Seurat.obj)
Seurat.obj <- RenameIdents(Seurat.obj, new.cluster.ids)
DimPlot(Seurat.obj, reduction = "umap",  pt.size = 1.2) + scale_color_manual(values = colors[c(1:6,10:12)]) + annotation_custom(grob)

pdf("~/Fig2c_heatmap_plot.pdf",width = 12,height = 8)
top8 <- markers_all_cells %>% group_by(cluster) %>% top_n(n = 8, wt = avg_logFC)
Fig2_c <- DoHeatmap(Seurat.obj, features = top8$gene, group.colors = colors[c(1:6,10:12)]) + NoLegend()
plot(Fig2_c)
dev.off()

#### 5. Take cluster 6 and recluster ####
Seurat.obj.6 <- subset(Seurat.obj, idents="6")

Seurat.obj.6 <- NormalizeData(object = Seurat.obj.6, normalization.method = "LogNormalize", scale.factor = 1e4)
Seurat.obj.6 <- FindVariableFeatures(Seurat.obj.6, selection.method = "vst", nfeatures = 2000)
Seurat.obj.6 <- ScaleData(Seurat.obj.6,vars.to.regress = c("nCount_RNA","percent.mito"))

Seurat.obj.6 <- RunPCA(Seurat.obj.6, features = VariableFeatures(Seurat.obj.6), npcs = 40)
ElbowPlot(Seurat.obj.6,ndims = 40)

#Use first 9 PC
n_PC = 9
Seurat.obj.6 <- FindNeighbors(Seurat.obj.6, dims = 1:n_PC)
Seurat.obj.6 <- FindClusters(Seurat.obj.6, resolution = 0.5)
#PCA
DimPlot(Seurat.obj.6, reduction = "pca")
n_cells = ncol(Seurat.obj.6)
grob <- grobTree(textGrob(paste0("n = ",n_cells), x=0.05,  y=0.95, hjust=0,gp=gpar(fontsize=22)))
#UMAP
Seurat.obj.6 <- RunUMAP(Seurat.obj.6, dims = 1:n_PC)
markers_6 <- FindAllMarkers(Seurat.obj.6, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(markers_6,file="~/markers.Subcluster6.csv")

new.cluster.ids <- c("6.0 CD4 Interm N-Eff", "6.1 Treg", "6.2 CD8 Cyto", "6.3 CD8 Prolif")
names(new.cluster.ids) <- levels(Seurat.obj.6)
Seurat.obj.6 <- RenameIdents(Seurat.obj.6, new.cluster.ids)

pdf("~/Fig2c_heatmap_subplot.pdf",width = 12,height = 8)
top10 <- markers_6 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
Fig2_c_subplot <- DoHeatmap(Seurat.obj.6, features = top10$gene, group.colors = colors[c(6:9)]) + NoLegend()
plot(Fig2_c_subplot)
dev.off()

#### 6. Reassign this clustering to the original dataset ####
old_clustering <- as.character(Idents(Seurat.obj))
old_clustering[which(names(Idents(Seurat.obj))%in%names(Idents(Seurat.obj.6)))] <- as.character(Idents(Seurat.obj.6))
names(old_clustering) <- names(Idents(Seurat.obj))
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = old_clustering, col.name = "new_clustering")

pdf("~/Fig2b.pdf",width = 12,height = 8)
Fig2_b <- DimPlot(Seurat.obj, reduction = "umap",  pt.size = 1.2, group.by = "new_clustering") + scale_color_manual(values = colors) + annotation_custom(grob)
plot(Fig2_b)
dev.off()

#### Plot expression of marker genes
genes_T_cell <- c("CD4","CD8A","CD8B","CD3E")
genes_naive <- c("SELL","CCR7","TCF7","IL7R")
genes_exh <- c("PDCD1", "LAG3", "HAVCR2", "KLRG1")
genes_exh2 <- c( "TIGIT", "CD244", "CD160", "BTLA")
genes_exh3 <- c( "CTLA4", "ENTPD1", "CD160", "ID2")
genes_exh4 <- c("CD274")
cyt_genes <- c("GZMB", "GZMK", "GZMH", "GZMA")
cyt_genes2 <- c("GZMM", "PRF1", "GNLY", "NKG7")
prolif_genes <- c("TOP2A","MKI67","TNF", "IFNG")
reg_genes <- c("MKI67","FOXP3","CD4","CD8A")

pdf("~/SuppFig2_c.pdf")
FeaturePlot(Seurat.obj,features = genes_T_cell, pt.size = 0.2, cols=c("lightgrey", "red"))
FeaturePlot(Seurat.obj,features = genes_naive, pt.size = 0.2, cols=c("lightgrey", "red"))
FeaturePlot(Seurat.obj,features = genes_exh, pt.size = 0.2, cols=c("lightgrey", "red"))
FeaturePlot(Seurat.obj,features = genes_exh2, pt.size = 0.2, cols=c("lightgrey", "red"))
FeaturePlot(Seurat.obj,features = genes_exh3, pt.size = 0.2, cols=c("lightgrey", "red"))
FeaturePlot(Seurat.obj,features = cyt_genes, pt.size = 0.2, cols=c("lightgrey", "red"))
FeaturePlot(Seurat.obj,features = cyt_genes2, pt.size = 0.2, cols=c("lightgrey", "red"))
FeaturePlot(Seurat.obj,features = reg_genes, pt.size = 0.2, cols=c("lightgrey", "red"))
FeaturePlot(Seurat.obj,features = genes_exh4, pt.size = 0.2, cols=c("lightgrey", "red"))
dev.off()

pdf("~/SuppFig2_b.pdf")
VlnPlot(Seurat.obj,features = c("CD27"), ncol=1, group.by = "new_clustering")
dev.off()

#### 7. Test the exhaustion signature ####
#Remove the degraded cluster
Seurat.obj.f <- subset(Seurat.obj, idents=levels(Idents(Seurat.obj))[c(1:8)])
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colors = getPalette(12)
exh_genes <- c("PDCD1", "LAG3", "HAVCR2", "KLRG1", "TIGIT", "CD244", "CD160", "BTLA", "CTLA4", "ENTPD1", "CD160", "ID2")
Seurat.obj.f <- AddModuleScore(object = Seurat.obj.f, features = list(exh_genes), name="Exh signature")

#Plot the score signatures
Seurat.obj.f$umap1 <- Seurat.obj.f@reductions$umap@cell.embeddings[,1]
Seurat.obj.f$umap2 <- Seurat.obj.f@reductions$umap@cell.embeddings[,2]
aux_df1 <- data.frame(umap1=Seurat.obj.f$umap1,umap2=Seurat.obj.f$umap2,clusters=Seurat.obj.f$new_clustering,Exh.signature=Seurat.obj.f$Exh.signature1)
ggplot(aux_df1, aes(x=umap1, y=umap2)) + geom_point(aes(color=Exh.signature),size = 1.2) + ggtitle("Exh.signature") + scale_colour_gradient(low = "white", high = "darkred") + NoAxes() + theme_bw()
pdf("~/Fig2c_exhausation_boxplot.pdf",width = 12,height = 8)
Fig2_c_exh <- ggplot(aux_df1, aes(x=fct_reorder(clusters, Exh.signature, .fun = median, .desc = TRUE), y=Exh.signature,fill=clusters)) + geom_boxplot(outlier.shape = NA) + theme_classic() + theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16)) + scale_fill_manual(values=colors[1:11]) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_cartesian(ylim=c(-0.4,0.5)) + xlab("")
plot(Fig2_c_exh)
dev.off()


#### 8. Incoroprate the TCR info ####
TCR_dict <- Dict$new()
Clonotype_size_dict <- Dict$new()
Clonotype_id_dict <- Dict$new()
TCR_UMIs_dict <- Dict$new()

# Control
Control_TCR <- read.csv(file =  "~/302_TCR/outs/filtered_contig_annotations.csv")
#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]

#Nº cells with TRA or TRB
Control_TCR_valid_n <- length(unique(Control_TCR_valid$barcode))
Control_TCR_clonotype <- read.csv(file="~/302_TCR/outs/clonotypes.csv")
Control_TCR_valid_clonotype <- merge(Control_TCR_valid,Control_TCR_clonotype,by.x="raw_clonotype_id",by.y="clonotype_id")
Control_TCR_valid_unique <- as.data.frame(table(Control_TCR_valid$barcode,Control_TCR_valid$chain))
Control_TCR_valid_unique2 <- Control_TCR_valid_unique[which(Control_TCR_valid_unique$Freq!=0),]

full_cont = 0
i <- 1
for(i in 1:nrow(Control_TCR_valid_unique2)){
  barcode_id <- paste0(substr(as.character(Control_TCR_valid_unique2[i,1]),1,nchar(as.character(Control_TCR_valid_unique2[i,1]))-2),"-1")
  tcr_info <- as.character(Control_TCR_valid_unique2[i,2])
  clone_size_info <- unique(Control_TCR_valid_clonotype$frequency[which(as.character(Control_TCR_valid_clonotype$barcode)%in%as.character(Control_TCR_valid_unique2[i,1]))])
  clonotype_id <- unique(Control_TCR_valid_clonotype$cdr3s_aa[which(as.character(Control_TCR_valid_clonotype$barcode)%in%as.character(Control_TCR_valid_unique2[i,1]))])
  UMIs <- unique(Control_TCR_valid_clonotype$umis[which(as.character(Control_TCR_valid_clonotype$barcode)%in%as.character(Control_TCR_valid_unique2[i,1]))])
  #Update the info related with this barcode
  #Check if the barcode is already in the dict
  if(TCR_dict$has(barcode_id)){
    TCR_dict$set(barcode_id,paste0(TCR_dict$peek(barcode_id),"/",tcr_info))
    full_cont = full_cont + 1
  } else{
    TCR_dict$add(barcode_id,tcr_info)
  }
  if(Clonotype_size_dict$has(barcode_id)){
    Clonotype_size_dict$set(barcode_id,clone_size_info)
  } else{
    Clonotype_size_dict$add(barcode_id,clone_size_info)
  }
  if(Clonotype_id_dict$has(barcode_id)){
    Clonotype_id_dict$set(barcode_id,clonotype_id)
  } else{
    Clonotype_id_dict$add(barcode_id,clonotype_id)
  }
  if(TCR_UMIs_dict$has(barcode_id)){
    TCR_UMIs_dict$set(barcode_id,UMIs)
  } else{
    TCR_UMIs_dict$add(barcode_id,UMIs)
  }
}
paste0("Number of cells with TCR (alpha or beta) reconstruction: ",i)
paste0("Number of cells with full TCR (alpha and beta) reconstruction: ",full_cont)

#Add the info to the Seurat object
n_cells <- length(colnames(Seurat.obj))

#For each barcode, check its status in TCR_dict
TCR_list <- NULL
Clono_size_list <- NULL
Clono_id_list <- NULL
TCR_UMIs_list <- NULL
i <- 1
cont <- 0
for(i in 1:ncol(Seurat.obj)){
  id <- colnames(Seurat.obj)[i]
  if(is.null(TCR_dict$peek(id))){
    TCR_list <- c(TCR_list,"No TCR") 
    Clono_size_list <- c(Clono_size_list,0)  
    Clono_id_list <- c(Clono_id_list,"No clonotype")  
    TCR_UMIs_list <- c(TCR_UMIs_list,0) 
  } else{
    # print(i)
    Clono_size <- Clonotype_size_dict$get(id)
    Clono_size_list <- c(Clono_size_list,Clono_size)
    Clono_id <- Clonotype_id_dict$get(id)
    Clono_id_list <- c(Clono_id_list,as.character(Clono_id))
    #Update TCR status according to Clono_id
    TCR_status <- ifelse(grepl("TRA",Clono_id),
                         ifelse(grepl("TRB",Clono_id),
                                "TRA/TRB",
                                "TRA"),
                         ifelse(grepl("TRB",Clono_id),
                                "TRB",
                                "No TCR"))
    # TCR_status <- TCR_dict$get(id)
    TCR_list <- c(TCR_list,TCR_status)
    TCR_UMIs_list <- c(TCR_UMIs_list,round(mean(TCR_UMIs_dict$get(id)),0)) 
    cont = cont + 1
  }
}

#Add this info to the seurat object
names(TCR_list) <- names(Seurat.obj$orig.ident)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = TCR_list, col.name = "TCR_info")
names(Clono_size_list) <- names(Seurat.obj$orig.ident)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = Clono_size_list, col.name = "Clono_size")
names(Clono_id_list) <- names(Seurat.obj$orig.ident)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = Clono_id_list, col.name = "Clono_id")
names(TCR_UMIs_list) <- names(Seurat.obj$orig.ident)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = TCR_UMIs_list, col.name = "TCR_UMIs")

actual_sizes <- table(Clono_id_list)
#Create another object with the actual sizes of the clonotypes present
Clono_size_list_updated <- Clono_size_list
# i <- 1
for(i in 1:length(Clono_size_list)){
  # print(i)
  id_cell <- names(Clono_size_list_updated[i])
  #Take the clonotype_id 
  clonotype_id <- Clono_id_list[which(names(Clono_id_list)%in%id_cell)]
  #Take the real size of each cell
  if(clonotype_id=="No clonotype"){
    new_size <- 0
  } else{
    new_size <- actual_sizes[which(names(actual_sizes)%in%clonotype_id)]
  }
  #Update the size
  Clono_size_list_updated[i] <- new_size
}

names(Clono_size_list_updated) <- names(Seurat.obj$orig.ident)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = Clono_size_list_updated, col.name = "Clono_size_updated")

#### 9. Plot which cells comes from M1 or Mets ####

### M1
#### Load the alpha and beta sequences
alpha_seq <- read.csv(file="/media/trincadojl/data/Projects/DEMATTOS/76_alpha.csv")
beta_seq <- read.csv(file="/media/trincadojl/data/Projects/DEMATTOS/76_beta.csv")

#### Load the clonotypes info from the sc
Control_TCR <- read.csv(file =  "/media/trincadojl/data/Projects/DEMATTOS/302_TCR_v2/outs/filtered_contig_annotations.csv")
#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]
#Separate between alpha and beta regions
Alpha_TCR_valid <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRA"),]
Beta_TCR_valid <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRB"),]

#### Compare clonotypes
alpha_merge <- merge(Alpha_TCR_valid,alpha_seq,by.x="cdr3",by.y="CDR3_AA")
beta_merge <- merge(Beta_TCR_valid,beta_seq,by.x="cdr3",by.y="CDR3_AA")
#How many cells in our Seurat.obj has an alpha/beta region coincidence in the bulk data?
match <- rep("No match",length(colnames(Seurat.obj)))
names(match) <- colnames(Seurat.obj)
match[which(names(match)%in%unique(alpha_merge$barcode))] <- "Alpha match"
match[which(names(match)%in%unique(beta_merge$barcode))] <- "Beta match"
match[which(names(match)%in%both_merge_unique_bc)] <- "Alpha/Beta match"
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = match, col.name = "M1_match")

### Mets
#### Load the sequences
metastasis_sequences <- read.table(file="/media/trincadojl/data/Projects/DEMATTOS/302_SingleFile_MiXCR_Para_JuanLu.tsv",sep="\t",header = TRUE)
TR_sequences <- metastasis_sequences[which(grepl("TR",metastasis_sequences$allVHitsWithScore)),]
TR_sequences$aaSeqCDR3_formatted <- gsub("_|\\*","",TR_sequences$aaSeqCDR3)
alpha_sequences <- metastasis_sequences[which(grepl("TRA",metastasis_sequences$allVHitsWithScore)),]
alpha_sequences$aaSeqCDR3_formatted <- gsub("_|\\*","",alpha_sequences$aaSeqCDR3)
beta_sequences <- metastasis_sequences[which(grepl("TRB",metastasis_sequences$allVHitsWithScore)),]
beta_sequences$aaSeqCDR3_formatted <- gsub("_|\\*","",beta_sequences$aaSeqCDR3)

#### Load the clonotypes info from the sc
Control_TCR <- read.csv(file =  "/media/trincadojl/data/Projects/DEMATTOS/302_TCR_v2/outs/filtered_contig_annotations.csv")
#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]
#Separate between alpha and beta regions
Alpha_TCR_valid <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRA"),]
Beta_TCR_valid <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRB"),]

#### Compare clonotypes
alpha_merge <- merge(Alpha_TCR_valid,alpha_sequences,by.x="cdr3",by.y="aaSeqCDR3_formatted")
beta_merge <- merge(Beta_TCR_valid,beta_sequences,by.x="cdr3",by.y="aaSeqCDR3_formatted")
#How many cells in our Seurat.obj has an alpha/beta region coincidence in the metastasis data?
match <- rep("No match",length(colnames(Seurat.obj)))
names(match) <- colnames(Seurat.obj)
match[which(names(match)%in%unique(alpha_merge$barcode))] <- "Alpha match"
match[which(names(match)%in%unique(beta_merge$barcode))] <- "Beta match"
match[which(names(match)%in%both_merge_unique_bc)] <- "Alpha/Beta match"
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = match, col.name = "Met_match")

### Plot on top of the UMAP
aux_df <- data.frame(M1=Seurat.obj$M1_match,Mets=Seurat.obj$Met_match)
aux_df$Origin <- "No match"
aux_df$Origin[which(aux_df$M1!="No match")] <- "M1"
aux_df$Origin[which(aux_df$Mets!="No match")] <- "Mets"
aux_df$Origin[which(aux_df$M1!="No match" & aux_df$Mets!="No match")] <- "M1/Mets"
Origin <- aux_df$Origin
names(Origin) <- rownames(aux_df)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = Origin, col.name = "Origin")
Seurat.obj$Origin <- factor(Seurat.obj$Origin,levels=c("No match","M1","Mets","M1/Mets"))

#Plot where these dots are in the UMAP
colors2 <- c('#bdbdbd', "#328bd9", '#cfa900', '#0ec244')
n_cells = ncol(Seurat.obj)
pt_sizes <- rep(0.6,n_cells)
pt_sizes[which(Seurat.obj$Origin!="No match")] <- 2
n_cells2 = length(which(Seurat.obj$Origin!="No match"))
grob2 <- grobTree(textGrob(paste0("n = ",n_cells2), x=0.05,  y=0.95, hjust=0,gp=gpar(fontsize=22)))
Fig2f <- DimPlot(Seurat.obj, reduction = "umap",  pt.size = pt_sizes, group.by = "Origin") + annotation_custom(grob2) + scale_color_manual(values = colors2)
ggsave(Fig2f, file='~/Fig2f.svg', device = "svg", width = 6, height = 6, scale = 0.6, dpi=300)

#### 10. Plot the expanded, singleton and non-tcr cells ####
n_cells = ncol(Seurat.obj)
Expanded <- rep("Non-TCR",n_cells)
Expanded[which(Seurat.obj$Clono_size>1)] <- "Expanded"
Expanded[which(Seurat.obj$Clono_size==1)] <- "Singleton"
table(Expanded)
names(Expanded) <- colnames(Seurat.obj)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = Expanded, col.name = "Expanded")

pt_sizes <- rep(0.6,n_cells)
pt_sizes[which(Seurat.obj$TCR_info!="No TCR")] <- 2
pt_sizes <- rep(0.4,n_cells)
pt_sizes[which(Seurat.obj$Clono_size==1)] <- 0.6
pt_sizes[which(Seurat.obj$Clono_size>=2)] <- 0.8
pt_sizes[which(Seurat.obj$Clono_size>=10)] <- 1
pt_sizes[which(Seurat.obj$Clono_size>=50)] <- 1.2
pt_sizes[which(Seurat.obj$Clono_size>=100)] <- 1.4
pt_sizes[which(Seurat.obj$Clono_size>=200)] <- 1.6

Seurat.obj$Expanded <- factor(Seurat.obj$Expanded, levels=c("Non-TCR","Singleton","Expanded"))
tplot = DimPlot(Seurat.obj, reduction = "umap",  group.by = "Expanded",  pt.size = pt_sizes) + 
  scale_color_manual(values = c("#A0A0A0","#F3B041","#900C3F"))
tplot[[1]]$layers[[1]]$aes_params$alpha = .8
pdf(file="~/Fig2_e_UMAP.svg",width = 10, height = 8)
plot(tplot)
dev.off()

#Create a barplot with the proportion of expanded clones from each tiempoint

#### Load the alpha and beta sequences (M1)
alpha_seq_M1 <- read.csv(file="~/Alpha_TCR_M1.csv")
beta_seq_M1 <- read.csv(file="~/Beta_TCR_M1.csv")
#Label them as expanded and non expanded
alpha_seq_M1$Expanded <- "Expanded"
alpha_seq_M1$Expanded[which(alpha_seq_M1$count==1)] <- "Singleton"
beta_seq_M1$Expanded <- "Expanded"
beta_seq_M1$Expanded[which(beta_seq_M1$count==1)] <- "Singleton"

#### Load the clonotypes info from the sc
Control_TCR <- read.csv(file =  "~/302_TCR/outs/filtered_contig_annotations.csv")
#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]
#We need to know if they are associated with an expanded clone
Control_TCR_clonotype <- read.csv(file="~/302_TCR_v2/outs/clonotypes.csv")
Control_TCR_valid_clonotype <- merge(Control_TCR_valid,Control_TCR_clonotype,by.x="raw_clonotype_id",by.y="clonotype_id")
Control_TCR_valid_clonotype$Expanded <- "Expanded"
Control_TCR_valid_clonotype$Expanded[which(Control_TCR_valid_clonotype$frequency==1)] <- "Singleton"
#Separate between alpha and beta regions
Alpha_TCR_valid_sc <- Control_TCR_valid_clonotype[which(Control_TCR_valid_clonotype$chain=="TRA"),]
Beta_TCR_valid_sc <- Control_TCR_valid_clonotype[which(Control_TCR_valid_clonotype$chain=="TRB"),]

#### Load the sequences (Mets)
metastasis_sequences <- read.table(file="~/TCR_Mets.tsv",sep="\t",header = TRUE)
metastasis_sequences$aaSeqCDR3_formatted <- gsub("_|\\*","",metastasis_sequences$aaSeqCDR3)
metastasis_sequences$Expanded <- "Expanded"
metastasis_sequences$Expanded[which(metastasis_sequences$cloneCount==1)] <- "Singleton"
TR_sequences <- metastasis_sequences[which(grepl("TR",metastasis_sequences$allVHitsWithScore)),]
alpha_sequences_Mets <- metastasis_sequences[which(grepl("TRA",metastasis_sequences$allVHitsWithScore)),]
alpha_sequences_Mets$aaSeqCDR3_formatted <- gsub("_|\\*","",alpha_sequences_Mets$aaSeqCDR3)
beta_sequences_Mets <- metastasis_sequences[which(grepl("TRB",metastasis_sequences$allVHitsWithScore)),]
beta_sequences_Mets$aaSeqCDR3_formatted <- gsub("_|\\*","",beta_sequences_Mets$aaSeqCDR3)

#### Do the barplot
alpha_M1 <- data.frame(Expanded=alpha_seq_M1$Expanded,Timepoint="M1")
beta_M1 <- data.frame(Expanded=beta_seq_M1$Expanded,Timepoint="M1")
alpha_sc <- data.frame(Expanded=Alpha_TCR_valid_sc$Expanded,Timepoint="sc")
beta_sc <- data.frame(Expanded=Beta_TCR_valid_sc$Expanded,Timepoint="sc")
alpha_Mets <- data.frame(Expanded=alpha_sequences_Mets$Expanded,Timepoint="Mets")
beta_Mets <- data.frame(Expanded=beta_sequences_Mets$Expanded,Timepoint="Mets")

all_alpha <- rbind(alpha_M1,alpha_sc,alpha_Mets)
all_alpha$Timepoint <- factor(all_alpha$Timepoint,levels = c("M1","sc","Mets"))
all_beta <- rbind(beta_M1,beta_sc,beta_Mets)
all_beta$Timepoint <- factor(all_beta$Timepoint,levels = c("M1","sc","Mets"))

barplot_alpha <- ggplot(all_alpha,aes(Timepoint,fill=Expanded)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#900C3F","#F3B041")) +
  theme_classic() +
  ylab("") + xlab("") +
  ggtitle("Alpha")

ggsave(barplot_alpha, file='~/Fig2e_barplot_alpha.svg', device = "svg", width = 6, height = 6, scale = 0.6, dpi=300)

barplot_beta <- ggplot(all_beta,aes(Timepoint,fill=Expanded)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#900C3F","#F3B041")) +
  theme_classic() +
  ylab("") + xlab("") +
  ggtitle("Beta")

ggsave(barplot_beta, file='~/Fig2e_barplot_beta.svg', device = "svg", width = 6, height = 6, scale = 0.6, dpi=300)

#### 11. Do a fishplot ####

#Take only the expanded clones

#Count the number of cells sharing each TCR region on sc
Control_TCR <- read.csv(file =  "~/302_TCR/outs/filtered_contig_annotations.csv")
#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]
#Separate between alpha and beta regions
Alpha_TCR_valid_sc <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRA"),]
Beta_TCR_valid_sc <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRB"),]
Alpha_TCR_valid_sc_prop <- Alpha_TCR_valid_sc %>% 
  group_by(cdr3) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
Beta_TCR_valid_sc_prop <- Beta_TCR_valid_sc %>% 
  group_by(cdr3) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
#There are too many non expanded clones. Remove the non expanded 
Alpha_TCR_valid_sc_prop_ex <- Alpha_TCR_valid_sc_prop[which(Alpha_TCR_valid_sc_prop$n>1),]
Beta_TCR_valid_sc_prop_ex <- Beta_TCR_valid_sc_prop[which(Beta_TCR_valid_sc_prop$n>1),]
#Get the proportions again
total = sum(Alpha_TCR_valid_sc_prop_ex$n)
Alpha_TCR_valid_sc_prop_ex$prop <- Alpha_TCR_valid_sc_prop_ex$n*100/total
total = sum(Beta_TCR_valid_sc_prop_ex$n)
Beta_TCR_valid_sc_prop_ex$prop <- Beta_TCR_valid_sc_prop_ex$n*100/total

#Get the proportion for each clone on M1
alpha_seq_M1_f <- unique(alpha_seq_M1[,c("CDR3_AA","count")])
total = sum(alpha_seq_M1_f$count)
alpha_seq_M1_f$prop <- alpha_seq_M1_f$count*100/total
beta_seq_M1_f <- unique(beta_seq_M1[,c("CDR3_AA","count")])
total = sum(beta_seq_M1_f$count)
beta_seq_M1_f$prop <- beta_seq_M1_f$count*100/total
#There are too many non expanded clones. Remove the non expanded 
alpha_seq_M1_f_ex <- alpha_seq_M1_f[which(alpha_seq_M1_f$count>1),]
total = sum(alpha_seq_M1_f_ex$count)
alpha_seq_M1_f_ex$prop <- alpha_seq_M1_f_ex$count*100/total
beta_seq_M1_f_ex <- beta_seq_M1_f[which(beta_seq_M1_f$count>1),]
total = sum(beta_seq_M1_f_ex$count)
beta_seq_M1_f_ex$prop <- beta_seq_M1_f_ex$count*100/total

#Get the proportion for each clone on Mets
alpha_sequences_Mets_f <- unique(alpha_sequences_Mets[,c("aaSeqCDR3_formatted","cloneCount")])
total = sum(alpha_sequences_Mets_f$cloneCount)
alpha_sequences_Mets_f$prop <- alpha_sequences_Mets_f$cloneCount*100/total
beta_sequences_Mets_f <- unique(beta_sequences_Mets[,c("aaSeqCDR3_formatted","cloneCount")])
total = sum(beta_sequences_Mets_f$cloneCount)
beta_sequences_Mets_f$prop <- beta_sequences_Mets_f$cloneCount*100/total
#There are too many non expanded clones. Remove the non expanded 
alpha_sequences_Mets_f_ex <- alpha_sequences_Mets_f[which(alpha_sequences_Mets_f$cloneCount>1),]
total = sum(alpha_sequences_Mets_f_ex$cloneCount)
alpha_sequences_Mets_f_ex$prop <- alpha_sequences_Mets_f_ex$cloneCount*100/total
beta_seq_M1_f_ex <- beta_seq_M1_f[which(beta_seq_M1_f$count>1),]
total = sum(beta_seq_M1_f_ex$count)
beta_seq_M1_f_ex$prop <- beta_seq_M1_f_ex$count*100/total

#There are still too many clones in M1 and Mets and very few on sc. Fishplot are goona look horrible. It's better just restrict on those clones that are shared between the sc clones and M1 and Mets
aux_df <- data.frame(M1=Seurat.obj$M1_match,Mets=Seurat.obj$Met_match,Clone=Seurat.obj$Clono_id)
aux_df$Origin <- "No match"
aux_df$Origin[which(aux_df$M1!="No match")] <- "M1"
aux_df$Origin[which(aux_df$Mets!="No match")] <- "Mets"
aux_df$Origin[which(aux_df$M1!="No match" & aux_df$Mets!="No match")] <- "M1/Mets"
#Take only the ones matching with either M1 or Mets
aux_df_f <- aux_df[which(aux_df$Origin!="No match"),]
#From the whole list of clones select only the ones on aux_df_f

#sc
Alpha_TCR_valid_sc_prop_f <- Alpha_TCR_valid_sc_prop[which(unlist(lapply(Alpha_TCR_valid_sc_prop$cdr3,function(x)any((grepl(x,aux_df_f$Clone))==TRUE)))),]
total = sum(Alpha_TCR_valid_sc_prop_f$n)
Alpha_TCR_valid_sc_prop_f$prop <- Alpha_TCR_valid_sc_prop_f$n*100/total
Beta_TCR_valid_sc_prop_f <- Beta_TCR_valid_sc_prop[which(unlist(lapply(Beta_TCR_valid_sc_prop$cdr3,function(x)any((grepl(x,aux_df_f$Clone))==TRUE)))),]
total = sum(Beta_TCR_valid_sc_prop_f$n)
Beta_TCR_valid_sc_prop_f$prop <- Beta_TCR_valid_sc_prop_f$n*100/total

#M1
alpha_seq_M1_f <- alpha_seq_M1[which(unlist(lapply(alpha_seq_M1$CDR3_AA,function(x)any((grepl(x,aux_df_f$Clone))==TRUE)))),]
total = sum(alpha_seq_M1_f$count)
alpha_seq_M1_f$prop <- alpha_seq_M1_f$count*100/total
beta_seq_M1_f <- beta_seq_M1[which(unlist(lapply(beta_seq_M1$CDR3_AA,function(x)any((grepl(x,aux_df_f$Clone))==TRUE)))),]
total = sum(beta_seq_M1_f$count)
beta_seq_M1_f$prop <- beta_seq_M1_f$count*100/total

#Mets
alpha_sequences_Mets_f <- alpha_sequences_Mets[which(unlist(lapply(alpha_sequences_Mets$aaSeqCDR3_formatted,function(x)any((grepl(x,aux_df_f$Clone))==TRUE)))),]
total = sum(alpha_sequences_Mets_f$cloneCount)
alpha_sequences_Mets_f$prop <- alpha_sequences_Mets_f$cloneCount*100/total
beta_sequences_Mets_f <- beta_sequences_Mets[which(unlist(lapply(beta_sequences_Mets$aaSeqCDR3_formatted,function(x)any((grepl(x,aux_df_f$Clone))==TRUE)))),]
total = sum(beta_sequences_Mets_f$cloneCount)
beta_sequences_Mets_f$prop <- beta_sequences_Mets_f$cloneCount*100/total

#Merge the clones from all timepoints
#Alpha
merge1 <- merge(Alpha_TCR_valid_sc_prop_f,alpha_seq_M1_f,by.x="cdr3",by.y="CDR3_AA",all=TRUE)
merge1 <- merge1[,c("cdr3", "n", "count")]
colnames(merge1) <- c("cdr3", "cloneCount.sc", "cloneCount.M1")
merge2 <- merge(merge1,alpha_sequences_Mets_f,by.x="cdr3",by.y="aaSeqCDR3_formatted",all=TRUE)
merge2 <- merge2[,c("cdr3", "cloneCount.sc", "cloneCount.M1", "cloneCount")]
colnames(merge2) <- c("cdr3", "cloneCount.sc", "cloneCount.M1", "cloneCount.Mets")
merge2[is.na(merge2)] <- 0
#There are some clones only appearing on a timepoint. Remove them
merge2_f <- merge2[which(!((merge2$cloneCount.sc!=0 & merge2$cloneCount.M1==0 & merge2$cloneCount.Mets==0) |
                             (merge2$cloneCount.sc==0 & merge2$cloneCount.M1!=0 & merge2$cloneCount.Mets==0) |
                             (merge2$cloneCount.sc==0 & merge2$cloneCount.M1==0 & merge2$cloneCount.Mets!=0))),]
#Label clones if they are showing on all timepoints or only in one
merge2_f$Origin <- "M1/Mets"
merge2_f$Origin[which(merge2_f$cloneCount.M1==0 & merge2_f$cloneCount.Mets!=0)] <- "Mets"
merge2_f$Origin[which(merge2_f$cloneCount.M1!=0 & merge2_f$cloneCount.Mets==0)] <- "M1"

#Get the proportions again
total.sc = sum(merge2_f$cloneCount.sc)
merge2_f$prop.sc <- merge2_f$cloneCount.sc*100/total.sc
total.M1 = sum(merge2_f$cloneCount.M1)
merge2_f$prop.M1 <- merge2_f$cloneCount.M1*100/total.M1
total.Mets = sum(merge2_f$cloneCount.Mets)
merge2_f$prop.Mets <- merge2_f$cloneCount.Mets*100/total.Mets

timepoints=c(0,45,90)      
frac.table = matrix(c(merge2_f$prop.M1,merge2_f$prop.sc,merge2_f$prop.Mets),
                    ncol=length(timepoints))
parents =  replicate(nrow(frac.table), 0) #all clonotypes are independent, parent 0
fish = createFishObject(frac.table,parents,timepoints=timepoints)
fish = layoutClones(fish, separate.independent.clones=T)
colors_cs <- randomColor(count = nrow(frac.table))
palette_color <-  c("#328bd9", '#cfa900', '#0ec244')
names(palette_color) <- c("M1","Mets","M1/Mets")
colors_cs <- palette_color[merge2_f$Origin]
#Color according to the status of the clone

fish = setCol(fish,colors_cs)
vlines=c(0,45,90)
pdf("~/Fig2g_Alpha_FishPlot.pdf",width = 8, height = 12)
fishPlot(fish, shape="spline", vlines=vlines, vlab=vlines,cex.vlab=0.9)
dev.off()

#Beta
merge1 <- merge(Beta_TCR_valid_sc_prop_f,beta_seq_M1_f,by.x="cdr3",by.y="CDR3_AA",all=TRUE)
merge1 <- merge1[,c("cdr3", "n", "count")]
colnames(merge1) <- c("cdr3", "cloneCount.sc", "cloneCount.M1")
merge2 <- merge(merge1,beta_sequences_Mets_f,by.x="cdr3",by.y="aaSeqCDR3_formatted",all=TRUE)
merge2 <- merge2[,c("cdr3", "cloneCount.sc", "cloneCount.M1", "cloneCount")]
colnames(merge2) <- c("cdr3", "cloneCount.sc", "cloneCount.M1", "cloneCount.Mets")
merge2[is.na(merge2)] <- 0
#There are some clones only appearing on a timepoint. Remove them
merge2_f <- merge2[which(!((merge2$cloneCount.sc!=0 & merge2$cloneCount.M1==0 & merge2$cloneCount.Mets==0) |
                             (merge2$cloneCount.sc==0 & merge2$cloneCount.M1!=0 & merge2$cloneCount.Mets==0) |
                             (merge2$cloneCount.sc==0 & merge2$cloneCount.M1==0 & merge2$cloneCount.Mets!=0))),]
#Label clones if they are showing on all timepoints or only in one
merge2_f$Origin <- "M1/Mets"
merge2_f$Origin[which(merge2_f$cloneCount.M1==0 & merge2_f$cloneCount.Mets!=0)] <- "Mets"
merge2_f$Origin[which(merge2_f$cloneCount.M1!=0 & merge2_f$cloneCount.Mets==0)] <- "M1"

#Get the proportions again
total.sc = sum(merge2_f$cloneCount.sc)
merge2_f$prop.sc <- merge2_f$cloneCount.sc*100/total.sc
total.M1 = sum(merge2_f$cloneCount.M1)
merge2_f$prop.M1 <- merge2_f$cloneCount.M1*100/total.M1
total.Mets = sum(merge2_f$cloneCount.Mets)
merge2_f$prop.Mets <- merge2_f$cloneCount.Mets*100/total.Mets

#Do the fishplot
timepoints=c(0,45,90)      
frac.table = matrix(c(merge2_f$prop.M1,merge2_f$prop.sc,merge2_f$prop.Mets),
                    ncol=length(timepoints))
parents =  replicate(nrow(frac.table), 0) #all clonotypes are independent, parent 0
fish = createFishObject(frac.table,parents,timepoints=timepoints)
fish = layoutClones(fish, separate.independent.clones=T)
colors_cs <- randomColor(count = nrow(frac.table))
palette_color <-  c("#328bd9", '#cfa900', '#0ec244')
names(palette_color) <- c("M1","Mets","M1/Mets")
colors_cs <- palette_color[merge2_f$Origin]
#Color according to the status of the clone

fish = setCol(fish,colors_cs)
vlines=c(0,45,90)
pdf("~/Fig2g_Beta_FishPlot.pdf",width = 8, height = 12)
fishPlot(fish, shape="spline", vlines=vlines, vlab=vlines,cex.vlab=0.9)
dev.off()










# #### 7. Create a table with the name of the cells, their % content, ident ... 
# aux_df <- data.frame(cellBC=colnames(Seurat.obj),cluster=Seurat.obj$new_clustering,UMIcounts=Seurat.obj$nCount_RNA,Ngenes=Seurat.obj$nFeature_RNA,mito_content=Seurat.obj$percent.mito,cell_cycle=Seurat.obj$Phase,TCR_seq=Seurat.obj$TCR_info,CDR3=Seurat.obj$Clono_id,TCR_UMIs=Seurat.obj$TCR_UMIs,Clone_size=Seurat.obj$Clono_size)
# write.csv(aux_df,file="~/sc_info_cells.csv",row.names = FALSE, quote = FALSE)


#### 12. Diversity ####
#### get a measure of the clonal diversity
entropy <- function(target) {
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}

#### Load the alpha and beta sequences from M1
alpha_TCR_M1 <- read.csv(file="~/Alpha_TCR_M1.csv")
beta_TCR_M1 <- read.csv(file="~/Beta_TCR_M1.csv")
#Get only CDR3 that starts with C and ends with F (like the sc TCR regions)
alpha_TCR_M1_f <- alpha_TCR_M1
beta_TCR_M1_f <- beta_TCR_M1

#### Load the alpha and beta sequences from the metastasis
metastasis_sequences <- read.table(file="~/TCR_Mets.tsv",sep="\t",header = TRUE)
TR_sequences <- metastasis_sequences[which(grepl("TR",metastasis_sequences$allVHitsWithScore)),]
TR_sequences$aaSeqCDR3_formatted <- gsub("_|\\*","",TR_sequences$aaSeqCDR3)
alpha_TCR_Mets <- metastasis_sequences[which(grepl("TRA",metastasis_sequences$allVHitsWithScore)),]
alpha_TCR_Mets$aaSeqCDR3_formatted <- gsub("_|\\*","",alpha_TCR_Mets$aaSeqCDR3)
# alpha_TCR_Mets_f <- alpha_TCR_Mets[which((substr(alpha_TCR_Mets$aaSeqCDR3,1,1)=="C") & (substr(alpha_TCR_Mets$aaSeqCDR3,nchar(alpha_TCR_Mets$aaSeqCDR3),nchar(alpha_TCR_Mets$aaSeqCDR3))=="F")),]
alpha_TCR_Mets_f <- alpha_TCR_Mets
beta_TCR_Mets <- metastasis_sequences[which(grepl("TRB",metastasis_sequences$allVHitsWithScore)),]
beta_TCR_Mets$aaSeqCDR3_formatted <- gsub("_|\\*","",beta_TCR_Mets$aaSeqCDR3)
# beta_TCR_Mets_f <- beta_TCR_Mets[which((substr(beta_TCR_Mets$aaSeqCDR3,1,1)=="C") & (substr(beta_TCR_Mets$aaSeqCDR3,nchar(beta_TCR_Mets$aaSeqCDR3),nchar(beta_TCR_Mets$aaSeqCDR3))=="F")),]
beta_TCR_Mets_f <- beta_TCR_Mets

#### Load the clonotypes info from the sc
Control_TCR <- read.csv(file =  "~/302_TCR/outs/filtered_contig_annotations.csv")
#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]
#Separate between alpha and beta regions
Alpha_TCR_sc <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRA"),]
Beta_TCR_sc <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRB"),]

#### Get the entropies
alpha_TCR_Mets_f2 <- alpha_TCR_Mets_f %>%
  group_by(Sample) %>%
  mutate(entropy_per_sample = entropy(aaSeqCDR3_formatted))
alpha_TCR_Mets_f2 <- unique(data.frame(sample=alpha_TCR_Mets_f2$Sample,entropy=alpha_TCR_Mets_f2$entropy_per_sample))
beta_TCR_Mets_f2 <- beta_TCR_Mets_f %>%
  group_by(Sample) %>%
  mutate(entropy_per_sample = entropy(aaSeqCDR3_formatted))
beta_TCR_Mets_f2 <- unique(data.frame(sample=beta_TCR_Mets_f2$Sample,entropy=beta_TCR_Mets_f2$entropy_per_sample))
all_TCR_Mets <- rbind(alpha_TCR_Mets_f2,beta_TCR_Mets_f2)
all_TCR_Mets$TCR <- c(rep(c("alpha"),18),rep(c("beta"),17))

#By tissue and sample of origin (M1/sc/Mets)
all_TCR_Mets
df_tissue2 <- data.frame(alpha_sc=entropy(Alpha_TCR_sc$cdr3),beta_sc=entropy(Beta_TCR_sc$cdr3),alpha_M1=entropy(alpha_TCR_M1_f$CDR3_AA),beta_M1=entropy(beta_TCR_M1_f$CDR3_AA))
melt.df_tissue2 <- melt(df_tissue2)
melt.df_tissue2$TCR <- rep(c("alpha","beta"),2)
melt.df_tissue2$variable <- factor(c("sc","sc","M1","M1"),levels=c("M1","sc"))
colnames(melt.df_tissue2)[1:2] <- c("sample","entropy")
melt.all <- rbind(all_TCR_Mets,melt.df_tissue2)
# melt.all$sample <- factor(melt.all$sample,levels=c("M1","sc","302_005_A","302_005_B","302_005_C","302_006_A","302_006_B","302_007","302_008_A","302_008_B","302_009","302_010_C","302_010_D","302_010_E","302_011_A","302_011_B","302_012_A","302_015","302_016","302_018"))
melt.all$sample <- factor(melt.all$sample,levels=c("M1","sc","M5A","M5B","M5C","M6A","M6B","M7","M8A","M8B","M9","M10C","M10D","M10E","M11A","M11B","M12A","M15","M16","M18"))
SuppfFig2_e <- ggplot(melt.all,aes(x=sample,y=entropy,color=TCR,group=TCR)) +
  geom_point() +
  geom_line(linetype="dashed") +
  theme_classic() +
  ylab("Shannon entropy") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
pdf("~/SuppFig2_e.pdf")
plot(SuppfFig2_e)
dev.off()

#### 13. Expansion ####

#### Load the alpha and beta sequences from M1
alpha_TCR_M1 <- read.csv(file="~/Alpha_TCR_M1.csv")
beta_TCR_M1 <- read.csv(file="~/Beta_TCR_M1.csv")
#Get only CDR3 that starts with C and ends with F (like the sc TCR regions)
alpha_TCR_M1_f <- alpha_TCR_M1
beta_TCR_M1_f <- beta_TCR_M1

#### Load the alpha and beta sequences from the metastasis
metastasis_sequences <- read.table(file="~/TCR_Mets.tsv",sep="\t",header = TRUE)
TR_sequences <- metastasis_sequences[which(grepl("TR",metastasis_sequences$allVHitsWithScore)),]
TR_sequences$aaSeqCDR3_formatted <- gsub("_|\\*","",TR_sequences$aaSeqCDR3)
alpha_TCR_Mets <- metastasis_sequences[which(grepl("TRA",metastasis_sequences$allVHitsWithScore)),]
alpha_TCR_Mets$aaSeqCDR3_formatted <- gsub("_|\\*","",alpha_TCR_Mets$aaSeqCDR3)
# alpha_TCR_Mets_f <- alpha_TCR_Mets[which((substr(alpha_TCR_Mets$aaSeqCDR3,1,1)=="C") & (substr(alpha_TCR_Mets$aaSeqCDR3,nchar(alpha_TCR_Mets$aaSeqCDR3),nchar(alpha_TCR_Mets$aaSeqCDR3))=="F")),]
alpha_TCR_Mets_f <- alpha_TCR_Mets
beta_TCR_Mets <- metastasis_sequences[which(grepl("TRB",metastasis_sequences$allVHitsWithScore)),]
beta_TCR_Mets$aaSeqCDR3_formatted <- gsub("_|\\*","",beta_TCR_Mets$aaSeqCDR3)
# beta_TCR_Mets_f <- beta_TCR_Mets[which((substr(beta_TCR_Mets$aaSeqCDR3,1,1)=="C") & (substr(beta_TCR_Mets$aaSeqCDR3,nchar(beta_TCR_Mets$aaSeqCDR3),nchar(beta_TCR_Mets$aaSeqCDR3))=="F")),]
beta_TCR_Mets_f <- beta_TCR_Mets

#### Load the clonotypes info from the sc
Control_TCR <- read.csv(file =  "~/302_TCR/outs/filtered_contig_annotations.csv")
#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]
#Separate between alpha and beta regions
Alpha_TCR_sc <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRA"),]
Beta_TCR_sc <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRB"),]

#### Get the clonal size per tiempoint
n_clonotypes.alpha_M1<- length(as.character(unique(alpha_TCR_M1_f$CDR3_AA)))
n_clonotypes.beta_M1<- length(as.character(unique(beta_TCR_M1_f$CDR3_AA)))
n_clonotypes.alpha_sc<- length(as.character(unique(Alpha_TCR_sc$cdr3)))
n_clonotypes.beta_sc<- length(as.character(unique(Beta_TCR_sc$cdr3)))
n_clonotypes.alpha_mets<- length(as.character(unique(alpha_TCR_Mets_f$aaSeqCDR3_formatted)))
n_clonotypes.beta_mets<- length(as.character(unique(beta_TCR_Mets_f$aaSeqCDR3_formatted)))
n_clonotypes <- data.frame(Timepoint=c("M1","M1","sc","sc","Mets","Mets"),Type=rep(c("alpha","beta"),3),n_clones=c(n_clonotypes.alpha_M1,n_clonotypes.beta_M1,n_clonotypes.alpha_sc,n_clonotypes.beta_sc,n_clonotypes.alpha_mets,n_clonotypes.beta_mets))
n_clonotypes$Timepoint <- factor(n_clonotypes$Timepoint,levels=c("M1","sc","Mets"))
SuppFig2d <- ggplot(n_clonotypes,aes(x=Timepoint,y=n_clones,fill=Type)) +
  geom_bar(stat="identity",position=position_dodge()) +
  theme(axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank()) +
  xlab("")
pdf("~/SuppFig2d.pdf",width = 6, height = 4)
plot(SuppFig2d)
dev.off()

#### Get a degree of expansion
#For sc: get the number of cells are associated on each cdr3
Alpha_TCR_sc_counts <- Alpha_TCR_sc %>% 
  group_by(cdr3) %>%
  tally()
Beta_TCR_sc_counts <- Beta_TCR_sc %>% 
  group_by(cdr3) %>%
  tally()

#Count number of expanded clones on each condition
threshold = 2
n_expanded.alpha_M1<- length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.beta_M1<- length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.alpha_sc<- nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),])
n_expanded.beta_sc<- nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),])
n_expanded.alpha_mets<- length(as.character(unique(alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))
n_expanded.beta_mets<- length(as.character(unique(beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))

norm.alpha_M1 <- n_expanded.alpha_M1/n_TCR_cells.alpha_M1
norm.beta_M1 <- n_expanded.beta_M1/n_TCR_cells.beta_M1
norm.alpha_sc <- n_expanded.alpha_sc/n_TCR_cells.alpha_sc
norm.beta_sc <- n_expanded.beta_sc/n_TCR_cells.beta_sc
norm.alpha_mets <- n_expanded.alpha_mets/n_TCR_cells.alpha_mets
norm.beta_mets <- n_expanded.beta_mets/n_TCR_cells.beta_mets
n_expanded_2 <- data.frame(Timepoint=c("M1","M1","sc","sc","Mets","Mets"),Type=rep(c("alpha","beta"),3),n_clones=c(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets))
n_expanded_2$Timepoint <- factor(n_expanded_2$Timepoint,levels=c("M1","sc","Mets"))

threshold = 3
n_expanded.alpha_M1<- length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.beta_M1<- length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.alpha_sc<- nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),])
n_expanded.beta_sc<- nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),])
n_expanded.alpha_mets<- length(as.character(unique(alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))
n_expanded.beta_mets<- length(as.character(unique(beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))

norm.alpha_M1 <- n_expanded.alpha_M1/n_TCR_cells.alpha_M1
norm.beta_M1 <- n_expanded.beta_M1/n_TCR_cells.beta_M1
norm.alpha_sc <- n_expanded.alpha_sc/n_TCR_cells.alpha_sc
norm.beta_sc <- n_expanded.beta_sc/n_TCR_cells.beta_sc
norm.alpha_mets <- n_expanded.alpha_mets/n_TCR_cells.alpha_mets
norm.beta_mets <- n_expanded.beta_mets/n_TCR_cells.beta_mets
n_expanded_3 <- data.frame(Timepoint=c("M1","M1","sc","sc","Mets","Mets"),Type=rep(c("alpha","beta"),3),n_clones=c(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets))
n_expanded_3$Timepoint <- factor(n_expanded_3$Timepoint,levels=c("M1","sc","Mets"))

threshold = 4
n_expanded.alpha_M1<- length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.beta_M1<- length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.alpha_sc<- nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),])
n_expanded.beta_sc<- nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),])
n_expanded.alpha_mets<- length(as.character(unique(alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))
n_expanded.beta_mets<- length(as.character(unique(beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))

norm.alpha_M1 <- n_expanded.alpha_M1/n_TCR_cells.alpha_M1
norm.beta_M1 <- n_expanded.beta_M1/n_TCR_cells.beta_M1
norm.alpha_sc <- n_expanded.alpha_sc/n_TCR_cells.alpha_sc
norm.beta_sc <- n_expanded.beta_sc/n_TCR_cells.beta_sc
norm.alpha_mets <- n_expanded.alpha_mets/n_TCR_cells.alpha_mets
norm.beta_mets <- n_expanded.beta_mets/n_TCR_cells.beta_mets
n_expanded_4 <- data.frame(Timepoint=c("M1","M1","sc","sc","Mets","Mets"),Type=rep(c("alpha","beta"),3),n_clones=c(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets))
n_expanded_4$Timepoint <- factor(n_expanded_4$Timepoint,levels=c("M1","sc","Mets"))

threshold = 5
n_expanded.alpha_M1<- length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.beta_M1<- length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.alpha_sc<- nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),])
n_expanded.beta_sc<- nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),])
n_expanded.alpha_mets<- length(as.character(unique(alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))
n_expanded.beta_mets<- length(as.character(unique(beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))

norm.alpha_M1 <- n_expanded.alpha_M1/n_TCR_cells.alpha_M1
norm.beta_M1 <- n_expanded.beta_M1/n_TCR_cells.beta_M1
norm.alpha_sc <- n_expanded.alpha_sc/n_TCR_cells.alpha_sc
norm.beta_sc <- n_expanded.beta_sc/n_TCR_cells.beta_sc
norm.alpha_mets <- n_expanded.alpha_mets/n_TCR_cells.alpha_mets
norm.beta_mets <- n_expanded.beta_mets/n_TCR_cells.beta_mets
n_expanded_5 <- data.frame(Timepoint=c("M1","M1","sc","sc","Mets","Mets"),Type=rep(c("alpha","beta"),3),n_clones=c(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets))
n_expanded_5$Timepoint <- factor(n_expanded_5$Timepoint,levels=c("M1","sc","Mets"))

threshold = 10
n_expanded.alpha_M1<- length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.beta_M1<- length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.alpha_sc<- nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),])
n_expanded.beta_sc<- nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),])
n_expanded.alpha_mets<- length(as.character(unique(alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))
n_expanded.beta_mets<- length(as.character(unique(beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))

norm.alpha_M1 <- n_expanded.alpha_M1/n_TCR_cells.alpha_M1
norm.beta_M1 <- n_expanded.beta_M1/n_TCR_cells.beta_M1
norm.alpha_sc <- n_expanded.alpha_sc/n_TCR_cells.alpha_sc
norm.beta_sc <- n_expanded.beta_sc/n_TCR_cells.beta_sc
norm.alpha_mets <- n_expanded.alpha_mets/n_TCR_cells.alpha_mets
norm.beta_mets <- n_expanded.beta_mets/n_TCR_cells.beta_mets
n_expanded_10 <- data.frame(Timepoint=c("M1","M1","sc","sc","Mets","Mets"),Type=rep(c("alpha","beta"),3),n_clones=c(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets))
n_expanded_10$Timepoint <- factor(n_expanded_10$Timepoint,levels=c("M1","sc","Mets"))

threshold = 15
n_expanded.alpha_M1<- length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.beta_M1<- length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.alpha_sc<- nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),])
n_expanded.beta_sc<- nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),])
n_expanded.alpha_mets<- length(as.character(unique(alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))
n_expanded.beta_mets<- length(as.character(unique(beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))

norm.alpha_M1 <- n_expanded.alpha_M1/n_TCR_cells.alpha_M1
norm.beta_M1 <- n_expanded.beta_M1/n_TCR_cells.beta_M1
norm.alpha_sc <- n_expanded.alpha_sc/n_TCR_cells.alpha_sc
norm.beta_sc <- n_expanded.beta_sc/n_TCR_cells.beta_sc
norm.alpha_mets <- n_expanded.alpha_mets/n_TCR_cells.alpha_mets
norm.beta_mets <- n_expanded.beta_mets/n_TCR_cells.beta_mets
n_expanded_15 <- data.frame(Timepoint=c("M1","M1","sc","sc","Mets","Mets"),Type=rep(c("alpha","beta"),3),n_clones=c(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets))
n_expanded_15$Timepoint <- factor(n_expanded_15$Timepoint,levels=c("M1","sc","Mets"))

threshold = 20
n_expanded.alpha_M1<- length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.beta_M1<- length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.alpha_sc<- nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),])
n_expanded.beta_sc<- nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),])
n_expanded.alpha_mets<- length(as.character(unique(alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))
n_expanded.beta_mets<- length(as.character(unique(beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))

norm.alpha_M1 <- n_expanded.alpha_M1/n_TCR_cells.alpha_M1
norm.beta_M1 <- n_expanded.beta_M1/n_TCR_cells.beta_M1
norm.alpha_sc <- n_expanded.alpha_sc/n_TCR_cells.alpha_sc
norm.beta_sc <- n_expanded.beta_sc/n_TCR_cells.beta_sc
norm.alpha_mets <- n_expanded.alpha_mets/n_TCR_cells.alpha_mets
norm.beta_mets <- n_expanded.beta_mets/n_TCR_cells.beta_mets
n_expanded_20 <- data.frame(Timepoint=c("M1","M1","sc","sc","Mets","Mets"),Type=rep(c("alpha","beta"),3),n_clones=c(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets))
n_expanded_20$Timepoint <- factor(n_expanded_20$Timepoint,levels=c("M1","sc","Mets"))

threshold = 25
n_expanded.alpha_M1<- length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.beta_M1<- length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.alpha_sc<- nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),])
n_expanded.beta_sc<- nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),])
n_expanded.alpha_mets<- length(as.character(unique(alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))
n_expanded.beta_mets<- length(as.character(unique(beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))

norm.alpha_M1 <- n_expanded.alpha_M1/n_TCR_cells.alpha_M1
norm.beta_M1 <- n_expanded.beta_M1/n_TCR_cells.beta_M1
norm.alpha_sc <- n_expanded.alpha_sc/n_TCR_cells.alpha_sc
norm.beta_sc <- n_expanded.beta_sc/n_TCR_cells.beta_sc
norm.alpha_mets <- n_expanded.alpha_mets/n_TCR_cells.alpha_mets
norm.beta_mets <- n_expanded.beta_mets/n_TCR_cells.beta_mets
n_expanded_25 <- data.frame(Timepoint=c("M1","M1","sc","sc","Mets","Mets"),Type=rep(c("alpha","beta"),3),n_clones=c(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets))
n_expanded_25$Timepoint <- factor(n_expanded_25$Timepoint,levels=c("M1","sc","Mets"))

threshold = 50
n_expanded.alpha_M1<- length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.beta_M1<- length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA)))
n_expanded.alpha_sc<- nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),])
n_expanded.beta_sc<- nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),])
n_expanded.alpha_mets<- length(as.character(unique(alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))
n_expanded.beta_mets<- length(as.character(unique(beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),]$aaSeqCDR3_formatted)))

norm.alpha_M1 <- n_expanded.alpha_M1/n_TCR_cells.alpha_M1
norm.beta_M1 <- n_expanded.beta_M1/n_TCR_cells.beta_M1
norm.alpha_sc <- n_expanded.alpha_sc/n_TCR_cells.alpha_sc
norm.beta_sc <- n_expanded.beta_sc/n_TCR_cells.beta_sc
norm.alpha_mets <- n_expanded.alpha_mets/n_TCR_cells.alpha_mets
norm.beta_mets <- n_expanded.beta_mets/n_TCR_cells.beta_mets
n_expanded_50 <- data.frame(Timepoint=c("M1","M1","sc","sc","Mets","Mets"),Type=rep(c("alpha","beta"),3),n_clones=c(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets))
n_expanded_50$Timepoint <- factor(n_expanded_50$Timepoint,levels=c("M1","sc","Mets"))

#Merge all the output matrixes
n_expanded_all <- cbind(n_expanded_2$n_clones,n_expanded_3$n_clones,n_expanded_4$n_clones,n_expanded_5$n_clones,n_expanded_10$n_clones,n_expanded_15$n_clones,n_expanded_20$n_clones,n_expanded_25$n_clones,n_expanded_50$n_clones)
colnames(n_expanded_all) <- c("2","3","4","5","10","15","20","25","50")
n_expanded_all <- as.data.frame(n_expanded_all)
n_expanded_all$Timepoint <- c("M1","M1","sc","sc","Mets","Mets") 
n_expanded_all$Timepoint <- factor(n_expanded_all$Timepoint,levels=c("M1","sc","Mets"))
n_expanded_all$Type <- rep(c("alpha","beta"),3)
melt.n_expanded_all <- melt(n_expanded_all)
SuppFig2_f_1 <- ggplot(melt.n_expanded_all,aes(x=variable,y=value,colour=Timepoint,group=Timepoint)) +
  geom_point() +
  geom_line() +
  xlab("Threshold") +
  ylab("Degree of expansion") +
  theme(axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_grid(. ~ Type)
pdf("~/SuppFig2_f_1.pdf",width = 7, height = 4)
plot(SuppFig2_f_1)
dev.off()


#### 14. Expansion from each independet metastasis ####

#### Load the alpha and beta sequences from M1
alpha_TCR_M1 <- read.csv(file="/media/trincadojl/data/Projects/DEMATTOS/76_alpha.csv")
beta_TCR_M1 <- read.csv(file="/media/trincadojl/data/Projects/DEMATTOS/76_beta.csv")
#Get only CDR3 that starts with C and ends with F (like the sc TCR regions)
alpha_TCR_M1_f <- alpha_TCR_M1
beta_TCR_M1_f <- beta_TCR_M1

#### Load the alpha and beta sequences from the metastasis
metastasis_sequences <- read.table(file="/media/trincadojl/data/Projects/DEMATTOS/302_SingleFile_MiXCR_Para_JuanLu.tsv",sep="\t",header = TRUE)
TR_sequences <- metastasis_sequences[which(grepl("TR",metastasis_sequences$allVHitsWithScore)),]
TR_sequences$aaSeqCDR3_formatted <- gsub("_|\\*","",TR_sequences$aaSeqCDR3)
alpha_TCR_Mets <- metastasis_sequences[which(grepl("TRA",metastasis_sequences$allVHitsWithScore)),]
alpha_TCR_Mets$aaSeqCDR3_formatted <- gsub("_|\\*","",alpha_TCR_Mets$aaSeqCDR3)
# alpha_TCR_Mets_f <- alpha_TCR_Mets[which((substr(alpha_TCR_Mets$aaSeqCDR3,1,1)=="C") & (substr(alpha_TCR_Mets$aaSeqCDR3,nchar(alpha_TCR_Mets$aaSeqCDR3),nchar(alpha_TCR_Mets$aaSeqCDR3))=="F")),]
alpha_TCR_Mets_f <- alpha_TCR_Mets
beta_TCR_Mets <- metastasis_sequences[which(grepl("TRB",metastasis_sequences$allVHitsWithScore)),]
beta_TCR_Mets$aaSeqCDR3_formatted <- gsub("_|\\*","",beta_TCR_Mets$aaSeqCDR3)
# beta_TCR_Mets_f <- beta_TCR_Mets[which((substr(beta_TCR_Mets$aaSeqCDR3,1,1)=="C") & (substr(beta_TCR_Mets$aaSeqCDR3,nchar(beta_TCR_Mets$aaSeqCDR3),nchar(beta_TCR_Mets$aaSeqCDR3))=="F")),]
beta_TCR_Mets_f <- beta_TCR_Mets

#### Load the clonotypes info from the sc
Control_TCR <- read.csv(file =  "/media/trincadojl/data/Projects/DEMATTOS/302_TCR_v2/outs/filtered_contig_annotations.csv")
#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]
#Separate between alpha and beta regions
Alpha_TCR_sc <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRA"),]
Beta_TCR_sc <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRB"),]

#### Get the clonal size per tiempoint
n_clonotypes.alpha_M1<- data.frame(Sample="M1",length=length(as.character(unique(alpha_TCR_M1_f$CDR3_AA))),Type="alpha")
n_clonotypes.beta_M1<- data.frame(Sample="M1",length=length(as.character(unique(beta_TCR_M1_f$CDR3_AA))),Type="beta")
n_clonotypes.alpha_sc<- data.frame(Sample="sc",length=length(as.character(unique(Alpha_TCR_sc$cdr3))),Type="alpha")
n_clonotypes.beta_sc<- data.frame(Sample="sc",length=length(as.character(unique(Beta_TCR_sc$cdr3))),Type="beta")
n_clonotypes.alpha_mets <- alpha_TCR_Mets_f %>% 
  group_by(Sample) %>%
  summarise(length=length(unique(aaSeqCDR3_formatted)))
n_clonotypes.alpha_mets$Type <- "alpha"
n_clonotypes.beta_mets <- beta_TCR_Mets_f %>% 
  group_by(Sample) %>%
  summarise(length=length(unique(aaSeqCDR3_formatted)))
n_clonotypes.beta_mets$Type <- "beta"
n_clonotypes <- rbind(n_clonotypes.alpha_M1,n_clonotypes.beta_M1,n_clonotypes.alpha_sc,n_clonotypes.beta_sc,n_clonotypes.alpha_mets,n_clonotypes.beta_mets)

#### Same as previous point, but normalizing to the total TCR cells recovered
n_TCR_cells.alpha_M1 <-  data.frame(Sample="M1",sum=sum(alpha_TCR_M1_f$count),Type="alpha")
n_TCR_cells.beta_M1 <-  data.frame(Sample="M1",sum=sum(beta_TCR_M1_f$count),Type="beta")
n_TCR_cells.alpha_sc <- data.frame(Sample="sc",sum=nrow(Alpha_TCR_sc),Type="alpha") 
n_TCR_cells.beta_sc <- data.frame(Sample="sc",sum=nrow(Beta_TCR_sc),Type="beta") 
n_TCR_cells.alpha_mets <- alpha_TCR_Mets_f %>% 
  group_by(Sample) %>%
  summarise(sum=sum(cloneCount))
n_TCR_cells.alpha_mets$Type <- "alpha"
n_TCR_cells.beta_mets <- beta_TCR_Mets_f %>% 
  group_by(Sample) %>%
  summarise(sum=sum(cloneCount))
n_TCR_cells.beta_mets$Type <- "beta"

norm.alpha_M1 <- cbind(n_clonotypes.alpha_M1,n_TCR_cells.alpha_M1)
norm.alpha_M1$norm <- norm.alpha_M1$length/norm.alpha_M1$sum
norm.beta_M1 <- cbind(n_clonotypes.beta_M1,n_TCR_cells.beta_M1)
norm.beta_M1$norm <- norm.beta_M1$length/norm.beta_M1$sum
norm.alpha_sc <- cbind(n_clonotypes.alpha_sc,n_TCR_cells.alpha_sc)
norm.alpha_sc$norm <- norm.alpha_sc$length/norm.alpha_sc$sum
norm.beta_sc <- cbind(n_clonotypes.beta_sc,n_TCR_cells.beta_sc)
norm.beta_sc$norm <- norm.beta_sc$length/norm.beta_sc$sum
norm.alpha_mets <- cbind(n_clonotypes.alpha_mets,n_TCR_cells.alpha_mets)
norm.alpha_mets$norm <- norm.alpha_mets$length/norm.alpha_mets$sum
norm.beta_mets <- cbind(n_clonotypes.beta_mets,n_TCR_cells.beta_mets)
norm.beta_mets$norm <- norm.beta_mets$length/norm.beta_mets$sum

n_clonotypes <- rbind(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets)
n_clonotypes <- n_clonotypes[,-3]

#### Get a degree of expansion
#Save this info in an external xlsx
file = "~/TCR_clones.xlsx"
write.xlsx(alpha_TCR_M1_f, file, sheetName = "Alpha_M1", col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(beta_TCR_M1_f, file, sheetName = "Beta_M1", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(Alpha_TCR_sc_counts, file, sheetName = "Alpha_sc", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(Beta_TCR_sc_counts, file, sheetName = "Beta_sc", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(alpha_TCR_Mets_f, file, sheetName = "Alpha_Mets", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(beta_TCR_Mets_f, file, sheetName = "Beta_Mets", col.names = TRUE, row.names = TRUE, append = TRUE)

#Count number of expanded clones on each condition
threshold = 2
n_expanded.alpha_M1<- data.frame(Sample="M1",length=length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="alpha")
n_expanded.beta_M1<- data.frame(Sample="M1",length=length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="beta")
n_expanded.alpha_sc<- data.frame(Sample="sc",length=nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),]),Type="alpha")
n_expanded.beta_sc<- data.frame(Sample="sc",length=nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),]),Type="beta")
n_expanded.alpha_mets <- alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.alpha_mets$Type <- "alpha"
n_expanded.beta_mets <- beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.beta_mets$Type <- "beta"

norm.alpha_M1 <- merge(n_expanded.alpha_M1,n_TCR_cells.alpha_M1,by="Sample")
norm.alpha_M1$norm <- norm.alpha_M1$length/norm.alpha_M1$sum
norm.alpha_M1 <- norm.alpha_M1[,-3]
colnames(norm.alpha_M1)[4] <- "Type"
norm.beta_M1 <- merge(n_expanded.beta_M1,n_TCR_cells.beta_M1,by="Sample")
norm.beta_M1$norm <- norm.beta_M1$length/norm.beta_M1$sum
norm.beta_M1 <- norm.beta_M1[,-3]
colnames(norm.beta_M1)[4] <- "Type"
norm.alpha_sc <- merge(n_expanded.alpha_sc,n_TCR_cells.alpha_sc,by="Sample")
norm.alpha_sc$norm <- norm.alpha_sc$length/norm.alpha_sc$sum
norm.alpha_sc <- norm.alpha_sc[,-3]
colnames(norm.alpha_sc)[4] <- "Type"
norm.beta_sc <- merge(n_expanded.beta_sc,n_TCR_cells.beta_sc,by="Sample")
norm.beta_sc$norm <- norm.beta_sc$length/norm.beta_sc$sum
norm.beta_sc <- norm.beta_sc[,-3]
colnames(norm.beta_sc)[4] <- "Type"
norm.alpha_mets <- merge(n_expanded.alpha_mets,n_TCR_cells.alpha_mets,by="Sample")
norm.alpha_mets$norm <- norm.alpha_mets$length/norm.alpha_mets$sum
norm.alpha_mets <- norm.alpha_mets[,-3]
colnames(norm.alpha_mets)[4] <- "Type"
norm.beta_mets <- merge(n_expanded.beta_mets,n_TCR_cells.beta_mets,by="Sample")
norm.beta_mets$norm <- norm.beta_mets$length/norm.beta_mets$sum
norm.beta_mets <- norm.beta_mets[,-3]
colnames(norm.beta_mets)[4] <- "Type"

n_expanded_2 <- rbind(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets)

threshold = 3
n_expanded.alpha_M1<- data.frame(Sample="M1",length=length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="alpha")
n_expanded.beta_M1<- data.frame(Sample="M1",length=length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="beta")
n_expanded.alpha_sc<- data.frame(Sample="sc",length=nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),]),Type="alpha")
n_expanded.beta_sc<- data.frame(Sample="sc",length=nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),]),Type="beta")
n_expanded.alpha_mets <- alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.alpha_mets$Type <- "alpha"
n_expanded.beta_mets <- beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.beta_mets$Type <- "beta"

norm.alpha_M1 <- merge(n_expanded.alpha_M1,n_TCR_cells.alpha_M1,by="Sample")
norm.alpha_M1$norm <- norm.alpha_M1$length/norm.alpha_M1$sum
norm.alpha_M1 <- norm.alpha_M1[,-3]
colnames(norm.alpha_M1)[4] <- "Type"
norm.beta_M1 <- merge(n_expanded.beta_M1,n_TCR_cells.beta_M1,by="Sample")
norm.beta_M1$norm <- norm.beta_M1$length/norm.beta_M1$sum
norm.beta_M1 <- norm.beta_M1[,-3]
colnames(norm.beta_M1)[4] <- "Type"
norm.alpha_sc <- merge(n_expanded.alpha_sc,n_TCR_cells.alpha_sc,by="Sample")
norm.alpha_sc$norm <- norm.alpha_sc$length/norm.alpha_sc$sum
norm.alpha_sc <- norm.alpha_sc[,-3]
colnames(norm.alpha_sc)[4] <- "Type"
norm.beta_sc <- merge(n_expanded.beta_sc,n_TCR_cells.beta_sc,by="Sample")
norm.beta_sc$norm <- norm.beta_sc$length/norm.beta_sc$sum
norm.beta_sc <- norm.beta_sc[,-3]
colnames(norm.beta_sc)[4] <- "Type"
norm.alpha_mets <- merge(n_expanded.alpha_mets,n_TCR_cells.alpha_mets,by="Sample")
norm.alpha_mets$norm <- norm.alpha_mets$length/norm.alpha_mets$sum
norm.alpha_mets <- norm.alpha_mets[,-3]
colnames(norm.alpha_mets)[4] <- "Type"
norm.beta_mets <- merge(n_expanded.beta_mets,n_TCR_cells.beta_mets,by="Sample")
norm.beta_mets$norm <- norm.beta_mets$length/norm.beta_mets$sum
norm.beta_mets <- norm.beta_mets[,-3]
colnames(norm.beta_mets)[4] <- "Type"

n_expanded_3 <- rbind(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets)

threshold = 4
n_expanded.alpha_M1<- data.frame(Sample="M1",length=length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="alpha")
n_expanded.beta_M1<- data.frame(Sample="M1",length=length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="beta")
n_expanded.alpha_sc<- data.frame(Sample="sc",length=nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),]),Type="alpha")
n_expanded.beta_sc<- data.frame(Sample="sc",length=nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),]),Type="beta")
n_expanded.alpha_mets <- alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.alpha_mets$Type <- "alpha"
n_expanded.beta_mets <- beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.beta_mets$Type <- "beta"

norm.alpha_M1 <- merge(n_expanded.alpha_M1,n_TCR_cells.alpha_M1,by="Sample")
norm.alpha_M1$norm <- norm.alpha_M1$length/norm.alpha_M1$sum
norm.alpha_M1 <- norm.alpha_M1[,-3]
colnames(norm.alpha_M1)[4] <- "Type"
norm.beta_M1 <- merge(n_expanded.beta_M1,n_TCR_cells.beta_M1,by="Sample")
norm.beta_M1$norm <- norm.beta_M1$length/norm.beta_M1$sum
norm.beta_M1 <- norm.beta_M1[,-3]
colnames(norm.beta_M1)[4] <- "Type"
norm.alpha_sc <- merge(n_expanded.alpha_sc,n_TCR_cells.alpha_sc,by="Sample")
norm.alpha_sc$norm <- norm.alpha_sc$length/norm.alpha_sc$sum
norm.alpha_sc <- norm.alpha_sc[,-3]
colnames(norm.alpha_sc)[4] <- "Type"
norm.beta_sc <- merge(n_expanded.beta_sc,n_TCR_cells.beta_sc,by="Sample")
norm.beta_sc$norm <- norm.beta_sc$length/norm.beta_sc$sum
norm.beta_sc <- norm.beta_sc[,-3]
colnames(norm.beta_sc)[4] <- "Type"
norm.alpha_mets <- merge(n_expanded.alpha_mets,n_TCR_cells.alpha_mets,by="Sample")
norm.alpha_mets$norm <- norm.alpha_mets$length/norm.alpha_mets$sum
norm.alpha_mets <- norm.alpha_mets[,-3]
colnames(norm.alpha_mets)[4] <- "Type"
norm.beta_mets <- merge(n_expanded.beta_mets,n_TCR_cells.beta_mets,by="Sample")
norm.beta_mets$norm <- norm.beta_mets$length/norm.beta_mets$sum
norm.beta_mets <- norm.beta_mets[,-3]
colnames(norm.beta_mets)[4] <- "Type"

n_expanded_4 <- rbind(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets)

threshold = 5
n_expanded.alpha_M1<- data.frame(Sample="M1",length=length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="alpha")
n_expanded.beta_M1<- data.frame(Sample="M1",length=length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="beta")
n_expanded.alpha_sc<- data.frame(Sample="sc",length=nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),]),Type="alpha")
n_expanded.beta_sc<- data.frame(Sample="sc",length=nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),]),Type="beta")
n_expanded.alpha_mets <- alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.alpha_mets$Type <- "alpha"
n_expanded.beta_mets <- beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.beta_mets$Type <- "beta"

norm.alpha_M1 <- merge(n_expanded.alpha_M1,n_TCR_cells.alpha_M1,by="Sample")
norm.alpha_M1$norm <- norm.alpha_M1$length/norm.alpha_M1$sum
norm.alpha_M1 <- norm.alpha_M1[,-3]
colnames(norm.alpha_M1)[4] <- "Type"
norm.beta_M1 <- merge(n_expanded.beta_M1,n_TCR_cells.beta_M1,by="Sample")
norm.beta_M1$norm <- norm.beta_M1$length/norm.beta_M1$sum
norm.beta_M1 <- norm.beta_M1[,-3]
colnames(norm.beta_M1)[4] <- "Type"
norm.alpha_sc <- merge(n_expanded.alpha_sc,n_TCR_cells.alpha_sc,by="Sample")
norm.alpha_sc$norm <- norm.alpha_sc$length/norm.alpha_sc$sum
norm.alpha_sc <- norm.alpha_sc[,-3]
colnames(norm.alpha_sc)[4] <- "Type"
norm.beta_sc <- merge(n_expanded.beta_sc,n_TCR_cells.beta_sc,by="Sample")
norm.beta_sc$norm <- norm.beta_sc$length/norm.beta_sc$sum
norm.beta_sc <- norm.beta_sc[,-3]
colnames(norm.beta_sc)[4] <- "Type"
norm.alpha_mets <- merge(n_expanded.alpha_mets,n_TCR_cells.alpha_mets,by="Sample")
norm.alpha_mets$norm <- norm.alpha_mets$length/norm.alpha_mets$sum
norm.alpha_mets <- norm.alpha_mets[,-3]
colnames(norm.alpha_mets)[4] <- "Type"
norm.beta_mets <- merge(n_expanded.beta_mets,n_TCR_cells.beta_mets,by="Sample")
norm.beta_mets$norm <- norm.beta_mets$length/norm.beta_mets$sum
norm.beta_mets <- norm.beta_mets[,-3]
colnames(norm.beta_mets)[4] <- "Type"

n_expanded_5 <- rbind(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets)

threshold = 10
n_expanded.alpha_M1<- data.frame(Sample="M1",length=length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="alpha")
n_expanded.beta_M1<- data.frame(Sample="M1",length=length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="beta")
n_expanded.alpha_sc<- data.frame(Sample="sc",length=nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),]),Type="alpha")
n_expanded.beta_sc<- data.frame(Sample="sc",length=nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),]),Type="beta")
n_expanded.alpha_mets <- alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.alpha_mets$Type <- "alpha"
n_expanded.beta_mets <- beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.beta_mets$Type <- "beta"

norm.alpha_M1 <- merge(n_expanded.alpha_M1,n_TCR_cells.alpha_M1,by="Sample")
norm.alpha_M1$norm <- norm.alpha_M1$length/norm.alpha_M1$sum
norm.alpha_M1 <- norm.alpha_M1[,-3]
colnames(norm.alpha_M1)[4] <- "Type"
norm.beta_M1 <- merge(n_expanded.beta_M1,n_TCR_cells.beta_M1,by="Sample")
norm.beta_M1$norm <- norm.beta_M1$length/norm.beta_M1$sum
norm.beta_M1 <- norm.beta_M1[,-3]
colnames(norm.beta_M1)[4] <- "Type"
norm.alpha_sc <- merge(n_expanded.alpha_sc,n_TCR_cells.alpha_sc,by="Sample")
norm.alpha_sc$norm <- norm.alpha_sc$length/norm.alpha_sc$sum
norm.alpha_sc <- norm.alpha_sc[,-3]
colnames(norm.alpha_sc)[4] <- "Type"
norm.beta_sc <- merge(n_expanded.beta_sc,n_TCR_cells.beta_sc,by="Sample")
norm.beta_sc$norm <- norm.beta_sc$length/norm.beta_sc$sum
norm.beta_sc <- norm.beta_sc[,-3]
colnames(norm.beta_sc)[4] <- "Type"
norm.alpha_mets <- merge(n_expanded.alpha_mets,n_TCR_cells.alpha_mets,by="Sample")
norm.alpha_mets$norm <- norm.alpha_mets$length/norm.alpha_mets$sum
norm.alpha_mets <- norm.alpha_mets[,-3]
colnames(norm.alpha_mets)[4] <- "Type"
norm.beta_mets <- merge(n_expanded.beta_mets,n_TCR_cells.beta_mets,by="Sample")
norm.beta_mets$norm <- norm.beta_mets$length/norm.beta_mets$sum
norm.beta_mets <- norm.beta_mets[,-3]
colnames(norm.beta_mets)[4] <- "Type"

n_expanded_10 <- rbind(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets)

threshold = 15
n_expanded.alpha_M1<- data.frame(Sample="M1",length=length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="alpha")
n_expanded.beta_M1<- data.frame(Sample="M1",length=length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="beta")
n_expanded.alpha_sc<- data.frame(Sample="sc",length=nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),]),Type="alpha")
n_expanded.beta_sc<- data.frame(Sample="sc",length=nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),]),Type="beta")
n_expanded.alpha_mets <- alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.alpha_mets$Type <- "alpha"
n_expanded.beta_mets <- beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.beta_mets$Type <- "beta"

norm.alpha_M1 <- merge(n_expanded.alpha_M1,n_TCR_cells.alpha_M1,by="Sample")
norm.alpha_M1$norm <- norm.alpha_M1$length/norm.alpha_M1$sum
norm.alpha_M1 <- norm.alpha_M1[,-3]
colnames(norm.alpha_M1)[4] <- "Type"
norm.beta_M1 <- merge(n_expanded.beta_M1,n_TCR_cells.beta_M1,by="Sample")
norm.beta_M1$norm <- norm.beta_M1$length/norm.beta_M1$sum
norm.beta_M1 <- norm.beta_M1[,-3]
colnames(norm.beta_M1)[4] <- "Type"
norm.alpha_sc <- merge(n_expanded.alpha_sc,n_TCR_cells.alpha_sc,by="Sample")
norm.alpha_sc$norm <- norm.alpha_sc$length/norm.alpha_sc$sum
norm.alpha_sc <- norm.alpha_sc[,-3]
colnames(norm.alpha_sc)[4] <- "Type"
norm.beta_sc <- merge(n_expanded.beta_sc,n_TCR_cells.beta_sc,by="Sample")
norm.beta_sc$norm <- norm.beta_sc$length/norm.beta_sc$sum
norm.beta_sc <- norm.beta_sc[,-3]
colnames(norm.beta_sc)[4] <- "Type"
norm.alpha_mets <- merge(n_expanded.alpha_mets,n_TCR_cells.alpha_mets,by="Sample")
norm.alpha_mets$norm <- norm.alpha_mets$length/norm.alpha_mets$sum
norm.alpha_mets <- norm.alpha_mets[,-3]
colnames(norm.alpha_mets)[4] <- "Type"
norm.beta_mets <- merge(n_expanded.beta_mets,n_TCR_cells.beta_mets,by="Sample")
norm.beta_mets$norm <- norm.beta_mets$length/norm.beta_mets$sum
norm.beta_mets <- norm.beta_mets[,-3]
colnames(norm.beta_mets)[4] <- "Type"

n_expanded_15 <- rbind(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets)

threshold = 20
n_expanded.alpha_M1<- data.frame(Sample="M1",length=length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="alpha")
n_expanded.beta_M1<- data.frame(Sample="M1",length=length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="beta")
n_expanded.alpha_sc<- data.frame(Sample="sc",length=nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),]),Type="alpha")
n_expanded.beta_sc<- data.frame(Sample="sc",length=nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),]),Type="beta")
n_expanded.alpha_mets <- alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.alpha_mets$Type <- "alpha"
n_expanded.beta_mets <- beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.beta_mets$Type <- "beta"

norm.alpha_M1 <- merge(n_expanded.alpha_M1,n_TCR_cells.alpha_M1,by="Sample")
norm.alpha_M1$norm <- norm.alpha_M1$length/norm.alpha_M1$sum
norm.alpha_M1 <- norm.alpha_M1[,-3]
colnames(norm.alpha_M1)[4] <- "Type"
norm.beta_M1 <- merge(n_expanded.beta_M1,n_TCR_cells.beta_M1,by="Sample")
norm.beta_M1$norm <- norm.beta_M1$length/norm.beta_M1$sum
norm.beta_M1 <- norm.beta_M1[,-3]
colnames(norm.beta_M1)[4] <- "Type"
norm.alpha_sc <- merge(n_expanded.alpha_sc,n_TCR_cells.alpha_sc,by="Sample")
norm.alpha_sc$norm <- norm.alpha_sc$length/norm.alpha_sc$sum
norm.alpha_sc <- norm.alpha_sc[,-3]
colnames(norm.alpha_sc)[4] <- "Type"
norm.beta_sc <- merge(n_expanded.beta_sc,n_TCR_cells.beta_sc,by="Sample")
norm.beta_sc$norm <- norm.beta_sc$length/norm.beta_sc$sum
norm.beta_sc <- norm.beta_sc[,-3]
colnames(norm.beta_sc)[4] <- "Type"
norm.alpha_mets <- merge(n_expanded.alpha_mets,n_TCR_cells.alpha_mets,by="Sample")
norm.alpha_mets$norm <- norm.alpha_mets$length/norm.alpha_mets$sum
norm.alpha_mets <- norm.alpha_mets[,-3]
colnames(norm.alpha_mets)[4] <- "Type"
norm.beta_mets <- merge(n_expanded.beta_mets,n_TCR_cells.beta_mets,by="Sample")
norm.beta_mets$norm <- norm.beta_mets$length/norm.beta_mets$sum
norm.beta_mets <- norm.beta_mets[,-3]
colnames(norm.beta_mets)[4] <- "Type"

n_expanded_20 <- rbind(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets)

threshold = 25
n_expanded.alpha_M1<- data.frame(Sample="M1",length=length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="alpha")
n_expanded.beta_M1<- data.frame(Sample="M1",length=length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="beta")
n_expanded.alpha_sc<- data.frame(Sample="sc",length=nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),]),Type="alpha")
n_expanded.beta_sc<- data.frame(Sample="sc",length=nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),]),Type="beta")
n_expanded.alpha_mets <- alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.alpha_mets$Type <- "alpha"
n_expanded.beta_mets <- beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.beta_mets$Type <- "beta"

norm.alpha_M1 <- merge(n_expanded.alpha_M1,n_TCR_cells.alpha_M1,by="Sample")
norm.alpha_M1$norm <- norm.alpha_M1$length/norm.alpha_M1$sum
norm.alpha_M1 <- norm.alpha_M1[,-3]
colnames(norm.alpha_M1)[4] <- "Type"
norm.beta_M1 <- merge(n_expanded.beta_M1,n_TCR_cells.beta_M1,by="Sample")
norm.beta_M1$norm <- norm.beta_M1$length/norm.beta_M1$sum
norm.beta_M1 <- norm.beta_M1[,-3]
colnames(norm.beta_M1)[4] <- "Type"
norm.alpha_sc <- merge(n_expanded.alpha_sc,n_TCR_cells.alpha_sc,by="Sample")
norm.alpha_sc$norm <- norm.alpha_sc$length/norm.alpha_sc$sum
norm.alpha_sc <- norm.alpha_sc[,-3]
colnames(norm.alpha_sc)[4] <- "Type"
norm.beta_sc <- merge(n_expanded.beta_sc,n_TCR_cells.beta_sc,by="Sample")
norm.beta_sc$norm <- norm.beta_sc$length/norm.beta_sc$sum
norm.beta_sc <- norm.beta_sc[,-3]
colnames(norm.beta_sc)[4] <- "Type"
norm.alpha_mets <- merge(n_expanded.alpha_mets,n_TCR_cells.alpha_mets,by="Sample")
norm.alpha_mets$norm <- norm.alpha_mets$length/norm.alpha_mets$sum
norm.alpha_mets <- norm.alpha_mets[,-3]
colnames(norm.alpha_mets)[4] <- "Type"
norm.beta_mets <- merge(n_expanded.beta_mets,n_TCR_cells.beta_mets,by="Sample")
norm.beta_mets$norm <- norm.beta_mets$length/norm.beta_mets$sum
norm.beta_mets <- norm.beta_mets[,-3]
colnames(norm.beta_mets)[4] <- "Type"

n_expanded_25 <- rbind(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets)

threshold = 50
n_expanded.alpha_M1<- data.frame(Sample="M1",length=length(as.character(unique(alpha_TCR_M1_f[which(alpha_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="alpha")
n_expanded.beta_M1<- data.frame(Sample="M1",length=length(as.character(unique(beta_TCR_M1_f[which(beta_TCR_M1_f$count>=threshold),]$CDR3_AA))),Type="beta")
n_expanded.alpha_sc<- data.frame(Sample="sc",length=nrow(Alpha_TCR_sc_counts[which(Alpha_TCR_sc_counts$n>=threshold),]),Type="alpha")
n_expanded.beta_sc<- data.frame(Sample="sc",length=nrow(Beta_TCR_sc_counts[which(Beta_TCR_sc_counts$n>=threshold),]),Type="beta")
n_expanded.alpha_mets <- alpha_TCR_Mets_f[which(alpha_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.alpha_mets$Type <- "alpha"
n_expanded.beta_mets <- beta_TCR_Mets_f[which(beta_TCR_Mets_f$cloneCount>=threshold),] %>% 
  group_by(Sample) %>%
  summarise(length=length(as.character(unique(aaSeqCDR3_formatted))))
n_expanded.beta_mets$Type <- "beta"

norm.alpha_M1 <- merge(n_expanded.alpha_M1,n_TCR_cells.alpha_M1,by="Sample")
norm.alpha_M1$norm <- norm.alpha_M1$length/norm.alpha_M1$sum
norm.alpha_M1 <- norm.alpha_M1[,-3]
colnames(norm.alpha_M1)[4] <- "Type"
norm.beta_M1 <- merge(n_expanded.beta_M1,n_TCR_cells.beta_M1,by="Sample")
norm.beta_M1$norm <- norm.beta_M1$length/norm.beta_M1$sum
norm.beta_M1 <- norm.beta_M1[,-3]
colnames(norm.beta_M1)[4] <- "Type"
norm.alpha_sc <- merge(n_expanded.alpha_sc,n_TCR_cells.alpha_sc,by="Sample")
norm.alpha_sc$norm <- norm.alpha_sc$length/norm.alpha_sc$sum
norm.alpha_sc <- norm.alpha_sc[,-3]
colnames(norm.alpha_sc)[4] <- "Type"
norm.beta_sc <- merge(n_expanded.beta_sc,n_TCR_cells.beta_sc,by="Sample")
norm.beta_sc$norm <- norm.beta_sc$length/norm.beta_sc$sum
norm.beta_sc <- norm.beta_sc[,-3]
colnames(norm.beta_sc)[4] <- "Type"
norm.alpha_mets <- merge(n_expanded.alpha_mets,n_TCR_cells.alpha_mets,by="Sample")
norm.alpha_mets$norm <- norm.alpha_mets$length/norm.alpha_mets$sum
norm.alpha_mets <- norm.alpha_mets[,-3]
colnames(norm.alpha_mets)[4] <- "Type"
norm.beta_mets <- merge(n_expanded.beta_mets,n_TCR_cells.beta_mets,by="Sample")
norm.beta_mets$norm <- norm.beta_mets$length/norm.beta_mets$sum
norm.beta_mets <- norm.beta_mets[,-3]
colnames(norm.beta_mets)[4] <- "Type"

n_expanded_50 <- rbind(norm.alpha_M1,norm.beta_M1,norm.alpha_sc,norm.beta_sc,norm.alpha_mets,norm.beta_mets)

#Merge all the output matrixes
aux1 <- merge(n_expanded_2,n_expanded_3,all=TRUE,by=c("Sample","Type"),suffixes = c("_2","_3"))
aux2 <- merge(aux1,n_expanded_4,all=TRUE,by=c("Sample","Type"))
colnames(aux2)[9:11] <- paste0(colnames(aux2)[9:11],"_4")
aux3 <- merge(aux2,n_expanded_5,all=TRUE,by=c("Sample","Type"))
colnames(aux3)[12:14] <- paste0(colnames(aux3)[12:14],"_5")
aux4 <- merge(aux3,n_expanded_10,all=TRUE,by=c("Sample","Type"))
colnames(aux4)[15:17] <- paste0(colnames(aux4)[15:17],"_6")
aux5 <- merge(aux4,n_expanded_15,all=TRUE,by=c("Sample","Type"))
colnames(aux5)[18:20] <- paste0(colnames(aux5)[18:20],"_7")
aux6 <- merge(aux5,n_expanded_20,all=TRUE,by=c("Sample","Type"))
colnames(aux6)[21:23] <- paste0(colnames(aux6)[21:23],"_8")
# aux7 <- merge(aux6,n_expanded_25,all=TRUE,by=c("Sample","Type"))
# colnames(aux7)[24:26] <- paste0(colnames(aux7)[24:26],"_9")
# aux8 <- merge(aux7,n_expanded_50,all=TRUE,by=c("Sample","Type"))
# colnames(aux8)[27:29] <- paste0(colnames(aux8)[27:29],"_10")
n_expanded_all <- aux6
n_expanded_all[is.na(n_expanded_all)] <- 0
melt.n_expanded_all <- melt(n_expanded_all)
melt.n_expanded_all <- melt.n_expanded_all[which(grepl("norm",melt.n_expanded_all$variable)),]

getPalette = colorRampPalette(brewer.pal(12, "Set1"))
colors = getPalette(20)

linetype<-c(c(1:10),c(1:10))
shape<-c(rep(c(1:6),3),c(1:2))
linetype_df <- data.frame(Sample=unique(melt.n_expanded_all$Sample),linetype=linetype,shape=shape)
melt.n_expanded_all2 <- merge(melt.n_expanded_all,linetype_df,by="Sample")
melt.n_expanded_all2$Sample <- factor(melt.n_expanded_all2$Sample,levels=c("M1","sc","302_005_A","302_005_B","302_005_C","302_006_A","302_006_B","302_007","302_008_A","302_008_B","302_009","302_010_C","302_010_D","302_010_E","302_011_A","302_011_B","302_012_A","302_015","302_016","302_018"))
melt.n_expanded_all2$Sample2 <- "M1"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="sc")] <- "sc"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_005_A")] <- "M5A"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_005_B")] <- "M5B"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_005_C")] <- "M5C"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_006_A")] <- "M6A"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_006_B")] <- "M6B"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_007")] <- "M7"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_008_A")] <- "M8A"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_008_B")] <- "M8B"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_009")] <- "M9"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_010_C")] <- "M10C"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_010_D")] <- "M10D"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_010_E")] <- "M10E"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_011_A")] <- "M11A"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_011_B")] <- "M11B"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_012_A")] <- "M12A"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_015")] <- "M15"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_016")] <- "M16"
melt.n_expanded_all2$Sample2[which(melt.n_expanded_all2$Sample=="302_018")] <- "M18"
melt.n_expanded_all2$Sample2 <- factor(melt.n_expanded_all2$Sample2,levels=c("M1","sc","M5A","M5B","M5C","M6A","M6B","M7","M8A","M8B","M9","M10C","M10D","M10E","M11A","M11B","M12A","M15","M16","M18"))

SuppFig2f_2 <- ggplot(melt.n_expanded_all2,aes(x=variable,y=value,group=Sample2)) +
  geom_point(aes(shape=as.factor(shape),color=Sample2), size=2) +
  geom_line(aes(linetype=as.factor(linetype),color=Sample2), size=1) +
  xlab("Threshold") +
  ylab("Degree of expansion") +
  theme(axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_grid(Type ~.) +
  scale_color_manual(values = colors) + 
  scale_x_discrete(labels=c("2","3","4","5","10","15","20"))

pdf("~/SuppFig2f_2.pdf",width = 8, height = 7)
plot(SuppFig2f_2)
dev.off()


#### 15. Convert to h5ad (this is necessary to run jupyter notebook Fig2.ipynb) ####
#Converting from Seurat to AnnData via h5Seurat
SaveH5Seurat(Seurat.obj, filename = "~/Seurat.obj.h5Seurat",overwrite = TRUE)
Convert("~/Seurat.obj.h5Seurat", dest = "~/Seurat.obj.h5ad",overwrite = TRUE)

#### 16. Create the commands to generate the circa plots with pyCircos (this is necessary to run jupyter notebook Fig2.ipynb) ####
#### Load the alpha and beta sequences (M1)
M1_alpha_seq <- read.csv(file="~/Alpha_TCR_M1.csv")
M1_beta_seq <- read.csv(file="~/Beta_TCR_M1.csv")

#### Load the clonotypes info from the sc
Control_TCR <- read.csv(file =  "~/302_TCR/outs/filtered_contig_annotations.csv")
#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]
#Separate between alpha and beta regions
Sc_alpha_seq <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRA"),]
Sc_beta_seq <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRB"),]
#Get the clone size for each alpha and beta clone
Sc_alpha_seq_clone_sizes <- as.data.frame(table(Sc_alpha_seq$cdr3))
# Sc_alpha_seq2 <- merge(Sc_alpha_seq,Sc_alpha_seq_clone_sizes,by.x="cdr3",by.y="Var1")
Sc_beta_seq_clone_sizes <- as.data.frame(table(Sc_beta_seq$cdr3))
# Sc_beta_seq2 <- merge(Sc_beta_seq,Sc_beta_seq_clone_sizes,by.x="cdr3",by.y="Var1")

#### Load the clonotypes info from the Metastasis
metastasis_sequences <- read.table(file="~/TCR_Mets.tsv",sep="\t",header = TRUE)
TR_sequences <- metastasis_sequences[which(grepl("TR",metastasis_sequences$allVHitsWithScore)),]
TR_sequences$aaSeqCDR3_formatted <- gsub("_|\\*","",TR_sequences$aaSeqCDR3)
Mets_alpha_sequences <- metastasis_sequences[which(grepl("TRA",metastasis_sequences$allVHitsWithScore)),]
Mets_alpha_sequences$aaSeqCDR3_formatted <- gsub("_|\\*","",Mets_alpha_sequences$aaSeqCDR3)
Mets_beta_sequences <- metastasis_sequences[which(grepl("TRB",metastasis_sequences$allVHitsWithScore)),]
Mets_beta_sequences$aaSeqCDR3_formatted <- gsub("_|\\*","",Mets_beta_sequences$aaSeqCDR3)

#### Put all clones together, with their clone sizes
#Sum repeated clones on M1
M1_alpha_seq_f <- as.data.frame(M1_alpha_seq %>% 
                                  group_by(CDR3_AA) %>%
                                  summarise(sum=sum(count)))
colnames(M1_alpha_seq_f) <- c("CDR3_AA","count")
M1_alpha_seq_f$id <- "M1"
M1_beta_seq_f <- as.data.frame(M1_beta_seq %>% 
                                 group_by(CDR3_AA) %>%
                                 summarise(sum=sum(count)))
colnames(M1_beta_seq_f) <- c("CDR3_AA","count")
M1_beta_seq_f$id <- "M1"
Sc_alpha_seq_f <- Sc_alpha_seq_clone_sizes
colnames(Sc_alpha_seq_f) <- c("CDR3_AA","count")
Sc_alpha_seq_f$id <- "Sc"
Sc_beta_seq_f <- Sc_beta_seq_clone_sizes
colnames(Sc_beta_seq_f) <- c("CDR3_AA","count")
Sc_beta_seq_f$id <- "Sc"
Mets_alpha_sequences_f <- unique(Mets_alpha_sequences[,c("aaSeqCDR3_formatted","cloneCount","Sample")])
colnames(Mets_alpha_sequences_f) <- c("CDR3_AA","count","id")
Mets_beta_sequences_f <- unique(Mets_beta_sequences[,c("aaSeqCDR3_formatted","cloneCount","Sample")])
colnames(Mets_beta_sequences_f) <- c("CDR3_AA","count","id")

alpha_clones <- rbind(M1_alpha_seq_f,Sc_alpha_seq_f,Mets_alpha_sequences_f)
beta_clones <- rbind(M1_beta_seq_f,Sc_beta_seq_f,Mets_beta_sequences_f)

#Sum by clonal size
alpha_clones_total_size <- as.data.frame(alpha_clones %>% 
                                           group_by(id) %>%
                                           summarise(sum=sum(count)))
beta_clones_total_size <- beta_clones %>% 
  group_by(id) %>%
  summarise(sum=sum(count))

#### Iterate the file, generating the necessary commands for the circus plot

#Get the prop order of the clones scaling to 10000
#M1
alpha_M1_sum <- sum(M1_alpha_seq_f$count)
absolute_position <- 0
# relative_position <- 0
total_count_accumulated <- NULL
# i <- 3
for(i in 1:nrow(M1_alpha_seq_f)){
  # pos_new <- ceiling(absolute_position*10000/alpha_M1_sum)
  size <- M1_alpha_seq_f$count[i]
  total_count_accumulated <- c(total_count_accumulated,absolute_position)
  absolute_position <- absolute_position + size
  # relative_position <- relative_position + M1_alpha_seq_f$count[i]
}
M1_alpha_seq_f$absolute_position <- total_count_accumulated
#Scale the abs position to 10000
M1_alpha_seq_f$scaled_pos <- unlist(lapply(M1_alpha_seq_f$absolute_position,function(x)ceiling(x*10000/alpha_M1_sum)))

beta_M1_sum <- sum(M1_beta_seq_f$count)
absolute_position <- 0
# relative_position <- 0
total_count_accumulated <- NULL
# i <- 3
for(i in 1:nrow(M1_beta_seq_f)){
  # pos_new <- ceiling(absolute_position*10000/beta_M1_sum)
  size <- M1_beta_seq_f$count[i]
  total_count_accumulated <- c(total_count_accumulated,absolute_position)
  absolute_position <- absolute_position + size
  # relative_position <- relative_position + M1_beta_seq_f$count[i]
}
M1_beta_seq_f$absolute_position <- total_count_accumulated
#Scale the abs position to 10000
M1_beta_seq_f$scaled_pos <- unlist(lapply(M1_beta_seq_f$absolute_position,function(x)ceiling(x*10000/beta_M1_sum)))

#Sc
alpha_Sc_sum <- sum(Sc_alpha_seq_f$count)
absolute_position <- 0
# relative_position <- 0
total_count_accumulated <- NULL
# i <- 3
for(i in 1:nrow(Sc_alpha_seq_f)){
  # pos_new <- ceiling(absolute_position*10000/alpha_Sc_sum)
  size <- Sc_alpha_seq_f$count[i]
  total_count_accumulated <- c(total_count_accumulated,absolute_position)
  absolute_position <- absolute_position + size
  # relative_position <- relative_position + Sc_alpha_seq_f$count[i]
}
Sc_alpha_seq_f$absolute_position <- total_count_accumulated
#Scale the abs position to 10000
Sc_alpha_seq_f$scaled_pos <- unlist(lapply(Sc_alpha_seq_f$absolute_position,function(x)ceiling(x*10000/alpha_Sc_sum)))

beta_Sc_sum <- sum(Sc_beta_seq_f$count)
absolute_position <- 0
# relative_position <- 0
total_count_accumulated <- NULL
# i <- 3
for(i in 1:nrow(Sc_beta_seq_f)){
  # pos_new <- ceiling(absolute_position*10000/beta_Sc_sum)
  size <- Sc_beta_seq_f$count[i]
  total_count_accumulated <- c(total_count_accumulated,absolute_position)
  absolute_position <- absolute_position + size
  # relative_position <- relative_position + Sc_beta_seq_f$count[i]
}
Sc_beta_seq_f$absolute_position <- total_count_accumulated
#Scale the abs position to 10000
Sc_beta_seq_f$scaled_pos <- unlist(lapply(Sc_beta_seq_f$absolute_position,function(x)ceiling(x*10000/beta_Sc_sum)))

#Mets
#Scale the total size from Mets to 10000
alpha_clones_total_size_Mets <- alpha_clones_total_size[which(grepl("302",alpha_clones_total_size$id)),]
alpha_Mets_sum <- sum(alpha_clones_total_size_Mets$sum)
alpha_clones_total_size_Mets$propsum <- unlist(lapply(alpha_clones_total_size_Mets$sum,function(x)round(x*10000/alpha_Mets_sum,0)))
beta_clones_total_size_Mets <- beta_clones_total_size[which(grepl("302",beta_clones_total_size$id)),]
beta_Mets_sum <- sum(beta_clones_total_size_Mets$sum)
beta_clones_total_size_Mets$propsum <- unlist(lapply(beta_clones_total_size_Mets$sum,function(x)round(x*10000/beta_Mets_sum,0)))

output_df <- data.frame(CDR3_AA=NA,count=NA,id=NA,absolute_position=NA,scaled_pos=NA)
#For each sample in Mets
i <- 1
for(i in 1:nrow(alpha_clones_total_size_Mets)){
  sample <- alpha_clones_total_size_Mets$id[i]
  #Take all the rows relative to this samples
  split_df <- Mets_alpha_sequences_f[which(Mets_alpha_sequences_f$id==sample),]
  absolute_position <- 0
  total_count_accumulated <- NULL
  # j <- 1
  for(j in 1:nrow(split_df)){
    size <- split_df$count[j]
    total_count_accumulated <- c(total_count_accumulated,absolute_position)
    absolute_position <- absolute_position + size
  }
  split_df$absolute_position <- total_count_accumulated
  #Scale the abs position to the prop space we save for this clone
  prop_space <- alpha_clones_total_size_Mets$propsum[which(alpha_clones_total_size_Mets$id==sample)]
  tot_space <- sum(split_df$count)
  split_df$scaled_pos <- unlist(lapply(split_df$absolute_position,function(x)ceiling(x*prop_space/tot_space)))
  #Save it to the final df
  output_df <- rbind(output_df,split_df)
}
#Remove first empty line
output_df <- output_df[-1,]
Mets_alpha_sequences_f <- output_df

output_df <- data.frame(CDR3_AA=NA,count=NA,id=NA,absolute_position=NA,scaled_pos=NA)
#For each sample in Mets
i <- 1
for(i in 1:nrow(beta_clones_total_size_Mets)){
  sample <- beta_clones_total_size_Mets$id[i]
  #Take all the rows relative to this samples
  split_df <- Mets_beta_sequences_f[which(Mets_beta_sequences_f$id==sample),]
  absolute_position <- 0
  total_count_accumulated <- NULL
  # j <- 1
  for(j in 1:nrow(split_df)){
    size <- split_df$count[j]
    total_count_accumulated <- c(total_count_accumulated,absolute_position)
    absolute_position <- absolute_position + size
  }
  split_df$absolute_position <- total_count_accumulated
  #Scale the abs position to the prop space we save for this clone
  prop_space <- beta_clones_total_size_Mets$propsum[which(beta_clones_total_size_Mets$id==sample)]
  tot_space <- sum(split_df$count)
  split_df$scaled_pos <- unlist(lapply(split_df$absolute_position,function(x)ceiling(x*prop_space/tot_space)))
  #Save it to the final df
  output_df <- rbind(output_df,split_df)
}
#Remove first empty line
output_df <- output_df[-1,]
Mets_beta_sequences_f <- output_df

# Get the overlapping clones
#M1 - sc
alpha_merge1 <- merge(M1_alpha_seq_f,Sc_alpha_seq_f,by="CDR3_AA")
beta_merge1 <- merge(M1_beta_seq_f,Sc_beta_seq_f,by="CDR3_AA")
#sc - Mets
alpha_merge2 <- merge(Sc_alpha_seq_f,Mets_alpha_sequences_f,by="CDR3_AA")
beta_merge2 <- merge(Sc_beta_seq_f,Mets_beta_sequences_f,by="CDR3_AA")
#M1 - Mets
alpha_merge3 <- merge(M1_alpha_seq_f,Mets_alpha_sequences_f,by="CDR3_AA")
beta_merge3 <- merge(M1_beta_seq_f,Mets_beta_sequences_f,by="CDR3_AA")

#Create a color palette for the clones
getPalette = colorRampPalette(brewer.pal(12, "Set3"))
colors = getPalette(nrow(alpha_clones_total_size))

#Given the huge unbalance between the 3 tiempoints, assign the same space to the 3 timepoints. We will give the proportional space to each of the clones
sink("~/pycircus_commands_alpha.txt")                    
cat("gcircle = Gcircle()\n")
cat(paste0("gcircle.add_locus(\"M1\", 10000, bottom=900, height=100, facecolor=\"",colors[1],"\")\n"))
cat(paste0("gcircle.add_locus(\"sc\", 10000, bottom=900, height=100, facecolor=\"",colors[2],"\")\n"))
#Generate commands for displaying the Met clones
i <- 1
for(i in 1:nrow(alpha_clones_total_size_Mets)){
  cat(paste0("gcircle.add_locus(\"",alpha_clones_total_size_Mets$id[i],"\", ",alpha_clones_total_size_Mets$propsum[i],", bottom=900, height=100, facecolor=\"",colors[i+2],"\")\n"))
}
cat("gcircle.set_locus()\n")
cat("gcircle.save()\n")

#Create the clone links
#For each link, I have to go to the proportional positions files:
#M1: M1_alpha_seq_f, M1_beta_seq_f
#sc: Sc_alpha_seq_f, Sc_beta_seq_f
#Mets: Mets_alpha_sequences_f, Mets_beta_sequences_f 
for(i in 1:nrow(alpha_merge1)){
  clone <- alpha_merge1$CDR3_AA[i]
  #Pos inicio del clone en M1
  M1_pos <- M1_alpha_seq_f[which(M1_alpha_seq_f$CDR3_AA==clone),]
  #Pos inicio del clone en sc
  sc_pos <- Sc_alpha_seq_f[which(Sc_alpha_seq_f$CDR3_AA==clone),] 
  cat(paste0("gcircle.chord_plot([\"M1\",",M1_pos$scaled_pos,",",M1_pos$scaled_pos+M1_pos$count, "],[\"sc\",",sc_pos$scaled_pos,",",sc_pos$scaled_pos+sc_pos$count, "], bottom=900, alpha=1.0)\n"))
}

i <- 4
for(i in 1:nrow(alpha_merge2)){
  clone <- alpha_merge2$CDR3_AA[i]
  #Pos inicio del clone en Mets
  sample <- alpha_merge2$id.y[i]
  Mets_pos <- Mets_alpha_sequences_f[which(Mets_alpha_sequences_f$CDR3_AA==clone & Mets_alpha_sequences_f$id==sample),]
  #Pos inicio del clone en sc
  sc_pos <- Sc_alpha_seq_f[which(Sc_alpha_seq_f$CDR3_AA==clone),] 
  cat(paste0("gcircle.chord_plot([\"",sample,"\",",Mets_pos$scaled_pos,",",Mets_pos$scaled_pos+Mets_pos$count, "],[\"sc\",",sc_pos$scaled_pos,",",sc_pos$scaled_pos+sc_pos$count, "], color=\"#c4b000\", bottom=900, alpha=1.0)\n"))
}

i <- 1
for(i in 1:nrow(alpha_merge3)){
  clone <- alpha_merge3$CDR3_AA[i]
  #Pos inicio del clone en Mets
  sample <- alpha_merge3$id.y[i]
  Mets_pos <- Mets_alpha_sequences_f[which(Mets_alpha_sequences_f$CDR3_AA==clone & Mets_alpha_sequences_f$id==sample),]
  #Pos inicio del clone en M1
  M1_pos <- M1_alpha_seq_f[which(M1_alpha_seq_f$CDR3_AA==clone),]
  cat(paste0("gcircle.chord_plot([\"",sample,"\",",Mets_pos$scaled_pos,",",Mets_pos$scaled_pos+Mets_pos$count, "],[\"M1\",",M1_pos$scaled_pos,",",M1_pos$scaled_pos+M1_pos$count, "], color=\"#15ad3d\", bottom=900, alpha=1.0)\n"))
}
sink()


sink("~/pycircus_commands_beta.txt")                    
cat("gcircle = Gcircle()\n")
cat(paste0("gcircle.add_locus(\"M1\", 10000, bottom=900, height=100, facecolor=\"",colors[1],"\")\n"))
cat(paste0("gcircle.add_locus(\"sc\", 10000, bottom=900, height=100, facecolor=\"",colors[2],"\")\n"))
#Generate commands for displaying the Met clones
i <- 1
for(i in 1:nrow(beta_clones_total_size_Mets)){
  cat(paste0("gcircle.add_locus(\"",beta_clones_total_size_Mets$id[i],"\", ",beta_clones_total_size_Mets$propsum[i],", bottom=900, height=100, facecolor=\"",colors[i+2],"\")\n"))
}
cat("gcircle.set_locus()\n")
cat("gcircle.save()\n")

#Create the clone links
#For each link, I have to go to the proportional positions files:
#M1: M1_beta_seq_f, M1_beta_seq_f
#sc: Sc_beta_seq_f, Sc_beta_seq_f
#Mets: Mets_beta_sequences_f, Mets_beta_sequences_f 
for(i in 1:nrow(beta_merge1)){
  clone <- beta_merge1$CDR3_AA[i]
  #Pos inicio del clone en M1
  M1_pos <- M1_beta_seq_f[which(M1_beta_seq_f$CDR3_AA==clone),]
  #Pos inicio del clone en sc
  sc_pos <- Sc_beta_seq_f[which(Sc_beta_seq_f$CDR3_AA==clone),] 
  cat(paste0("gcircle.chord_plot([\"M1\",",M1_pos$scaled_pos,",",M1_pos$scaled_pos+M1_pos$count, "],[\"sc\",",sc_pos$scaled_pos,",",sc_pos$scaled_pos+sc_pos$count, "], bottom=900, alpha=1.0)\n"))
}

i <- 4
for(i in 1:nrow(beta_merge2)){
  clone <- beta_merge2$CDR3_AA[i]
  #Pos inicio del clone en Mets
  sample <- beta_merge2$id.y[i]
  Mets_pos <- Mets_beta_sequences_f[which(Mets_beta_sequences_f$CDR3_AA==clone & Mets_beta_sequences_f$id==sample),]
  #Pos inicio del clone en sc
  sc_pos <- Sc_beta_seq_f[which(Sc_beta_seq_f$CDR3_AA==clone),] 
  cat(paste0("gcircle.chord_plot([\"",sample,"\",",Mets_pos$scaled_pos,",",Mets_pos$scaled_pos+Mets_pos$count, "],[\"sc\",",sc_pos$scaled_pos,",",sc_pos$scaled_pos+sc_pos$count, "],  color=\"#c4b000\", bottom=900, alpha=1.0)\n"))
}

i <- 1
for(i in 1:nrow(beta_merge3)){
  clone <- beta_merge3$CDR3_AA[i]
  #Pos inicio del clone en Mets
  sample <- beta_merge3$id.y[i]
  Mets_pos <- Mets_beta_sequences_f[which(Mets_beta_sequences_f$CDR3_AA==clone & Mets_beta_sequences_f$id==sample),]
  #Pos inicio del clone en M1
  M1_pos <- M1_beta_seq_f[which(M1_beta_seq_f$CDR3_AA==clone),]
  cat(paste0("gcircle.chord_plot([\"",sample,"\",",Mets_pos$scaled_pos,",",Mets_pos$scaled_pos+Mets_pos$count, "],[\"M1\",",M1_pos$scaled_pos,",",M1_pos$scaled_pos+M1_pos$count, "], color=\"#15ad3d\", bottom=900, alpha=1.0)\n"))
}
sink()
