options(stringsAsFactors = F)
library(cowplot)
library(ggpubr)
library(dplyr)
library(readr)
library(tidyr)
library(ggforce)
library(pals)
library(pheatmap)
library(scales)
library(ggthemes)
library(Seurat)

## Figure 1

## Fig1b: inflammation-associated genes
inflammationScore <- function(x, cols=NULL){
  inflgenes <- unlist(strsplit("IFNG, IFNGR1, IFNGR2, IL10, IL12A, IL12B, IL12RB1, IL12RB2, IL13, IL17A, IL17F, IL18, IL18R1, IL18RAP, IL1A, IL1B, IL2, IL21, IL21R, IL22, IL23A, IL23R, IL2RG, IL4, IL4R, IL5, IL6, JUN, NFKB1, RELA, RORA, RORC, S100A8, S100A9, STAT1, STAT3, STAT4, STAT6, TGFB1, TGFB2, TGFB3, TNF", ", *"))
  Idents(x) <- x$sample
  groupExp <- AverageExpression(x, assays = 'SCT', slot = 'data')$SCT
  inflgenes <- inflgenes[inflgenes %in% rownames(groupExp)]
  gMeanExp <- groupExp[inflgenes, ]
  gMeanExp <- t(apply(gMeanExp, 1, function(x){(x - min(x)) / max(x)}))
  gMeanExp <- reshape2::melt(t(gMeanExp)) %>% 
    dplyr::rename(Sample="Var1", Gene="Var2")
  p <- ggboxplot(gMeanExp, x = "Sample", y = "value", color = "Sample", 
                 orientation = "horizontal", add = "jitter", 
                 ylab = "Scaled mean", palette = cols) 
  p
}
sampleCols <- setNames(c("#00AFBB", "#E7B800", "#E7B800", "#E7B800", "#FC4E07", "#FC4E07", "#FC4E07", "#FC4E07"),
                       c("NT-P2", "PT-P1", "PT-P2", "PT-P3", "HM-P1", "HM-P2", "HM-P3", "HM-P4"))
fig1b <- inflammationScore(seuratObj, sampleCols)
fig1b

## Fig1c: 29 clusters
clusterCols <- c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F", "#BD9E39", "#E7BA52", "#31A354", "#E41A1C", 
                          "#6BAED6", "#9ECAE1", "#AD494A", "#E7CB94", "#74C476", "#A1D99B", "#C7E9C0", "#99600F", 
                          "#E7298A", "#C3BC3F", "#D6616B", "#FF7F00", "#1B9E77", "#FDAE6B", "#B3823E", "#66A61E", 
                          "#F1788D", "#C6DBEF", "#E6550D", "#E7969C")
setNames(clusterCols) <- sprintf("%02d",0:28)
fig1c <- DimPlot(seuratObj, group.by='cluster', cols=clusterCols, pt.size=1, raster=F) 


## Fig1d 
Idents(seuratObj) <- seuratObj$cluster
avg <- AverageExpression(seuratObj, features=mmarkers, assays='SCT')$SCT 
sdat <- t(apply(avg, 1, scale))
colnames(sdat) <- colnames(avg)
phtm <- pheatmap::pheatmap(sdat, cluster_rows = F, silent = T)
cord <- phtm$tree_col$labels[phtm$tree_col$order]
cord <- data.frame(cluster=factor(PHmeta, levels = unique(PTmeta[PTorder])), 
                   id=factor(names(PHmeta),levels = cord)) %>% 
  arrange(cluster, id) %>% pull(id) %>% as.character()
geneOrder <- data.frame(cluster=factor(colnames(markerExp)[apply(markerExp, 1, which.max)], levels = cord))
geneOrder$gene <- unlist(lapply(levels(geneOrder$cluster), function(i){
  order(markerExp[which(as.character(geneOrder$cluster)==i),i], decreasing = T)
}))
markerExp <- markerExp[order(geneOrder$cluster, geneOrder$gene), ]

mxid <- factor(colnames(sdat)[apply(sdat, 1, which.max)], levels=colnames(sdat))
gord <- rownames(sdat)[order(mxid)]
gord <- intersect(rownames(markerExp), gord)
Idents(seuratObj) <- factor(seuratObj$cluster, levels = cord)
fig1d <- DotPlot(seuratObj, features = gord) + rotate() +
  scale_color_gradientn(colours=rev(brewer.rdylbu(20)), guide = "colourbar") +
  theme_minimal_grid() + theme(axis.text.x=element_text(angle=90, hjust=1))
fig1d

## Fig1e 
groupCols <- c(NT="#00AFBB", PT="#E7B800", HM="#FC4E07" )
cellstat <- table(seuratObj$sample, seuratObj$majorcluster)
cellstat <- as.data.frame(cellstat * 100 / rowSums(cellstat))
colnames(cellstat) <- c("sample", "cluster", "percentage")
cellstat$group <- factor(sampleGroups[as.character(cellstat$sample)], levels=names(groupCols))
fig1e <- ggbarplot(cellstat, x = "cluster", y = "percentage", fill = "group",
                   color = 'black', add = "mean_se", palette = groupCols,
                   position = position_dodge(0.8)) + rotate()

## Fig1g 
majorColors <- setNames(c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F", "#E41A1C", "#E7298A", "#C3BC3F", "#FF7F00", "#1B9E77", "#66A61E", "#F1788D"), 
                        c("Ductal cell", "T cell","Fibroblasts", "NKT", "MKI67+ ductal cell",
                          "Mast cell","Endocrine cell", "B cell","Plasma cell","Endothelial cell"))
fig1g1 <- DimPlot(seuratObj, group.by='majorcluster', cols=majorColors, pt.size=1, raster=F) 
fig1g1

fig1g2 <- FeaturePlot(seuratObj, features='DEgenes', reduction='umap', order=TRUE, 
                     cols=brewer.ylorrd(max(seuratObj$DEgenes)), 
                     min.cutoff='q1', pt.size=1, combine=T, raster = F)
fig1g2

