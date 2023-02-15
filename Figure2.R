## Figure 2

## "Ductal cell", "MKI67+ ductal cell"
epiObj <- subset(x = seuratObj, subset = majorcluster %in% c("Ductal cell","MKI67+ ductal cell"))
epiObj <- RunUMAP(epiObj, verbose=TRUE) %>% 
  FindNeighbors() %>% 
  FindClusters(resolution = 0.25)
DefaultAssay(epiObj) <- "SCT"


## Fig2a 
ecols <- setNames(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"), c("0", 1:5))
mkcols <- c("#222222", "#F3C300", "#875692", "#F38400")

epiObj$flag <- ifelse(epiObj@assays$SCT@data['EPCAM',] > 0 & 
                         epiObj@assays$SCT@data['MKI67',] > 0, 'EPCAM+MKI67+', 
                       ifelse(epiObj@assays$SCT@data['EPCAM',] > 0, 'EPCAM+',
                              ifelse(epiObj@assays$SCT@data['MKI67',] > 0, 'MKI67+','Other')))
epiObj$flag <- factor(epiObj$flag, levels = c("Other",'MKI67+','EPCAM+MKI67+','EPCAM+'))

cell_eb <- data.frame(epiObj@reductions[["umap"]]@cell.embeddings[,1:2], cluster=as.character(epiObj$cluster),
                      marker=epiObj$flag, group=as.character(epiObj$group), stringsAsFactors = F)
colnames(cell_eb) <- c('x','y','cluster','marker','group')
pp1 <- ggplot(data=cell_eb, aes(x,y)) + geom_point(aes(colour = factor(marker)), size=0.25) +
  scale_color_manual(values = mkcols) + theme_pubr() +
  NoAxes() + NoLegend()
pp2 <- ggplot(data=cell_eb, aes(x,y)) + geom_point(aes(colour = factor(cluster)), size=0.25) +
  scale_color_manual(values = ecols) + theme_pubr() +
  NoAxes() + NoLegend()

cell_eb <- data.frame(epiObj@reductions[["umap"]]@cell.embeddings[,1:2], cluster=as.character(epiObj$cluster),
                      group=as.character(epiObj$group), stringsAsFactors = F)
colnames(cell_eb) <- c('x','y','cluster', 'group')
pp3 <- ggplot(data=cell_eb, aes(x,y)) + geom_point(aes(colour = factor(group)), alpha=0.5, size=0.25) +
  scale_color_manual(values = alpha(groupCols, alpha = 0.3)) + theme_pubr() +
  NoAxes() + NoLegend()
df <- epiObj@meta.data %>% dplyr::count(cluster, group) %>% dplyr::rename(cluster = cluster)
df$n <- df$n/as.numeric(table(epiObj@meta.data$group)[df$group]) # normalize 组的细胞数量
coords <- cell_eb %>% group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), .groups = "keep") %>% as.data.frame
coords <- coords %>% mutate(x=ifelse(cluster=='6',x+1,x), y=ifelse(cluster=='6',y+1,y))
piedf <- reshape2::acast(df, cluster ~ group, value.var = "n")
piedf[is.na(piedf)] <- 0
piedf <- cbind(piedf, coords[match(rownames(piedf), coords$cluster), -1])
cellnumbers <- dplyr::count(cell_eb, cluster)
sizee <- cellnumbers[match(rownames(piedf), cellnumbers[,1]), 2]
piedf$radius <- scales::rescale(sizee, c(0.4,0.9))
pp3 <- pp3 + ggnewscale::new_scale_fill() + 
  geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=piedf,
                  cols=colnames(piedf)[1:(ncol(piedf)-3)],color=NA) +
  coord_equal() + scale_fill_manual(values = groupCols)


cell_eb <- data.frame(epiObj@reductions[["umap"]]@cell.embeddings[,1:2], cluster=as.character(epiObj$cluster),
                      group=as.character(epiObj$CNV), stringsAsFactors = F)
colnames(cell_eb) <- c('x','y','cluster', 'group')
pp4 <- ggplot(data=cell_eb, aes(x,y)) + geom_point(aes(colour = factor(group)), alpha=0.5, size=0.25) +
  scale_color_manual(values = alpha(cnvCols, alpha = 0.3)) + theme_pubr() +
  NoAxes() + NoLegend()
df <- epiObj@meta.data %>% dplyr::count(cluster, CNV) %>% dplyr::rename(cluster = cluster)
df$n <- df$n/as.numeric(table(epiObj@meta.data$CNV)[df$CNV]) # normalize 组的细胞数量
coords <- cell_eb %>% group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), .groups = "keep") %>% as.data.frame
coords <- coords %>% mutate(x=ifelse(cluster=='6',x+1,x), y=ifelse(cluster=='6',y+1,y))
piedf <- reshape2::acast(df, cluster ~ CNV, value.var = "n")
piedf[is.na(piedf)] <- 0
piedf <- cbind(piedf, coords[match(rownames(piedf), coords$cluster), -1])
cellnumbers <- dplyr::count(cell_eb, cluster)
sizee <- cellnumbers[match(rownames(piedf), cellnumbers[,1]), 2]
piedf$radius <- scales::rescale(sizee, c(0.4,0.9))
pp4 <- pp4 + ggnewscale::new_scale_fill() + 
  geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=piedf,
                  cols=colnames(piedf)[1:(ncol(piedf)-3)],color=NA) +
  coord_equal() + scale_fill_manual(values = cnvCols)


print(ggarrange(as_ggplot(get_legend(p1)), as_ggplot(get_legend(p2)), as_ggplot(get_legend(p3)), as_ggplot(get_legend(p4)), ncol=2, nrow=2))
print(ggarrange(pp1+NoLegend()+NoAxes()+ggtitle(NULL), 
                pp2+NoLegend()+NoAxes()+ggtitle(NULL),
                ncol=2, nrow=2))
print(ggarrange(pp3+NoLegend()+NoAxes()+ggtitle(NULL), 
                pp4+NoLegend()+NoAxes()+ggtitle(NULL), 
                ncol=2, nrow=2))

## Fig2b: epi_split_p_heatmap

