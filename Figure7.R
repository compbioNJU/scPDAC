## Figure 7
library(circlize)
library(corrplot)
library(RColorBrewer)
library(ComplexHeatmap)

pdf("Figure7.pdf", height=8.27, width=8.27)
op <- par(mfrow=c(3,3))
cellinter_stat <- NULL
plist <- list()
for(db in c("NT","PT","HM")){
  sig_means <- read.delim(sprintf("cellphonedb/out/%s/significant_means.txt",db), check.names = F)
  sub_data <- sig_means[rowSums(sig_means[,-c(1:12)], na.rm = T) > 0, c(2,13:ncol(sig_means))] %>% 
    droplevels() %>% dplyr::mutate(group=db)
  if(is.null(cellinter_stat)){
    cellinter_stat <- sub_data
  }else{
    cellinter_stat <- cellinter_stat %>% full_join(sub_data)
  }
  intersum <- sort(apply(sig_means[,-c(1:12)], 2, function(x){
    sum(!is.na(x))
  }))
  intersum[intersum < 5] <- 0
  matx <- matrix(0, nrow=length(mCols), ncol=length(mCols))
  rownames(matx) <- colnames(matx) <- names(mCols)
  intersum <- intersum / sum(intersum)
  lapply(names(intersum), function(i){
    x <- unlist(strsplit(i, split = "\\|"))
    matx[x[1],x[2]] <<- intersum[i]
  })
  ## matx[matx < 0.01] <- 0
  chordDiagram(matx, directional = 1, 
               grid.col = mCols,
               order = names(mCols), 
               transparency = 0.25,
               direction.type = c("arrows", "diffHeight"), 
               diffHeight  = -0.04,
               annotationTrack = "grid", 
               annotationTrackHeight = c(0.05, 0.1),
               ## link.arr.type = "big.arrow", 
               link.sort = TRUE, 
               link.largest.ontop = TRUE)
  # Add text and axis
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    
    circos.text(mean(xlim), ylim[1], sector.name, facing = "inside", 
                niceFacing = TRUE, adj = c(0.5, 0))
  }, bg.border = NA)
  title(db)
}
par(op)

op <- par(mfrow=c(3,3))
for(db in c("NT","PT","HM")){
  sig_means <- read.delim(sprintf("cellphonedb/out/%s/significant_means.txt",db), check.names = F)
  sum(!is.na(sig_means[,-c(1:12)]))
  intersum <- sort(apply(sig_means[,-c(1:12)], 2, function(x){
    sum(!is.na(x))
  }))
  matx <- matrix(0, nrow=length(mCols), ncol=length(mCols))
  rownames(matx) <- colnames(matx) <- names(mCols)
  intersum <- intersum / sum(intersum)
  lapply(names(intersum), function(i){
    x <- unlist(strsplit(i, split = "\\|"))
    matx[x[1],x[2]] <<- intersum[i]
  })
  ## matx[upper.tri(matx)] <- -matx[lower.tri(matx)]
  matx <- matx[rowSums(matx)>0, colSums(matx)>0]
  corrplot(matx, is.corr = FALSE, 
           method = "square", type = "full", 
           tl.col = "black", order = "original", 
           col = rev(brewer.pal(n = 11, name = "RdYlBu")))
  
  annotation_col <- data.frame(cluster=factor(colnames(matx), levels=colnames(matx)), 
                               row.names=colnames(matx))
  annotation_row <- data.frame(cluster=factor(rownames(matx), levels=rownames(matx)), 
                               row.names=rownames(matx))
  ann_colors <- list(cluster = mCols[rownames(matx)])
  
  mapal <- colorRampPalette(c(rev(brewer.blues(8))[-(1:3)],brewer.ylorrd(5)))(256)
  plist[[db]] <- ComplexHeatmap::pheatmap(log1p(matx), use_raster=F, border_color = NA, fontsize = 6, color = mapal, 
                                          labels_row = NULL, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F,
                                          annotation_colors=ann_colors, annotation_col=annotation_col, 
                                          annotation_row=annotation_row)
}
par(op)

lapply(plist, print)

cellinter_stat$group <- factor(cellinter_stat$group, levels = levels(seuratObj$group))
pheatmap::pheatmap(table(cellinter_stat$interacting_pair, cellinter_stat$group), 
                   border_color = "white", cluster_cols = F, fontsize = 5, cellwidth = 20)
for(i in colnames(cellinter_stat)[-c(1,ncol(cellinter_stat))]){
  pltd <- table((na.omit(cellinter_stat[,c('interacting_pair','group',i)]))[,1:2] %>% droplevels())
  if(nrow(pltd) > 1 & ncol(pltd) > 1 & length(unique(pltd))>1){
    pheatmap::pheatmap(pltd, main = i, border_color = "white", cluster_cols = F, 
                       fontsize = 5, cellheight=5, cellwidth = 20)
    ## print(p)
  }
}




intercells <- matrix(0, nrow = length(sCols), ncol=nrow(PH))
rownames(intercells) <- names(sCols)
colnames(intercells) <- rownames(PH)
for(db in unique(PH$group)){
  pps <- PH %>% filter(group==db) %>% pull(pair)
  xx <- cellinter_stat4[[db]][pps, -1]
  stat <- lapply(apply(xx, 1, function(x){
    unlist(strsplit(colnames(xx)[!is.na(x)], split="\\|"))
  }), function(x){
    x[x=="Endocrine cell"] <- "Acinar cell" ## Secretory cell
    table(x)
  })
  for(i in names(stat)){
    x <- stat[[i]]
    intercells[names(x),i] <- x
  }
}

annotation_row <- data.frame(cluster=names(sCols), row.names = names(sCols))
ann_colors <- list(cluster=sCols)
pltd <- apply(log1p(intercells), 2, scale, center=F)
rownames(pltd) <- rownames(intercells)
pltd[pltd>2] <- 2
p3 <- ComplexHeatmap::pheatmap(pltd, cluster_cols = F, 
                               annotation_row = annotation_row, 
                               annotation_colors = ann_colors,
                               cellwidth = 10, cellheight = 10)


print(phtm + ha)
print(ggarrange(i0+xlab(NULL)+NoLegend(),nrow=2,ncol=1))
dev.off() 


