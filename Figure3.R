## Figure 3

## Fig3a 
Idents(epiObj) <- factor(epiObj$cluster, levels = c("5",4,0,1,2,3))
fig3a <- VlnPlot(epiObj, features=c("MKI67","EPCAM","CEACAM6","CEACAM5","KLK7"), 
             stack = T, fill.by='ident', cols = ecols, pt.size=0) +
  NoLegend() 
print(fig3a)

## Fig3b
# see Velocity.ipynb as an example 

## Fig3c 
p1 <- monocle::plot_cell_trajectory(cds, cell_size=0.5, color_by="Pseudotime", show_branch_points=F) +
  gradient_color(kovesi.rainbow_bgyr_35_85_c73(100))
p2 <- monocle::plot_cell_trajectory(cds, cell_size=0.1, color_by="cluster", show_branch_points=F) + 
  scale_color_manual(values=cols)
p3 <- monocle::plot_cell_trajectory(cds, cell_size=0.1, color_by="group", show_branch_points=F) + 
  scale_color_manual(values=groupCols)
p4 <- monocle::plot_cell_trajectory(cds, cell_size=0.5, 
                                    color_by="cluster", show_branch_points=F) + 
  scale_color_manual(values=cols) + facet_wrap(~group)

fig3c <- ggarrange(p1 + NoAxes() + NoLegend(), 
                   p2 + NoAxes() + NoLegend(), 
                   p3 + NoAxes() + NoLegend(),
                   monocle::plot_cell_trajectory(cds, cell_size=0.1, color_by="State", show_branch_points=F) + 
                     scale_color_manual(values=stateColor) + NoAxes() + NoLegend(), 
                   monocle::plot_cell_trajectory(cds, cell_size=0.1, color_by="CNV", show_branch_points=F) + 
                     scale_color_manual(values=cnvCols) + NoAxes() + NoLegend(), 
                   plot_expression_trajectory('KLK7') + 
                     gradient_color(c('grey',brewer.reds(10))) + NoAxes() + NoLegend(), 
                   nrow=2, ncol=3)
fig3c


## Fig3e 
cvnCols <- setNames(c("#b9ca5d", "#00a2b3", "#f3a546", "#cf3e53"), 
                    c("Normal", "Low", "Medium", "High"))

pdf("Fig3e1.pdf", height=4, width=8.27)
op <- par(mar=rep(0.5,4), mfrow=c(4,1))
x <- table(pltd$index, pltd$cluster)
x <- x / rowSums(x)
rownames(x) <- NULL
barplot(t(x), col=ecols, border = NA, axes = F, main="Cluster")

x <- table(pltd$index, pltd$group)
x <- x / rowSums(x)
rownames(x) <- NULL
barplot(t(x), col=groupCols, border = NA, axes = F, main="Group")

pltd$CNV <- epiObj@meta.data[rownames(pData(cds)),'CNV']
x <- table(pltd$index, pltd$CNV)
x <- x / rowSums(x)
rownames(x) <- NULL
barplot(t(x), col=cnvCols, border = NA, axes = F, main="CNV")

x <- table(pltd$index, pltd$State)
x <- x / rowSums(x)
rownames(x) <- NULL
barplot(t(x), col=stateColor, border = NA, axes = F, main="State")
par(op)
dev.off()

epiObj@meta.data[,'bin'] <- pData(cds)[Cells(epiObj),'index']
Idents(epiObj) <- 'bin'
bindata <- AverageExpression(epiObj, assays = 'SCT', features=genex)$SCT %>% as.matrix()
avg <- t(apply(bindata, 1, function(x){x=log1p(x);x/max(x)}))
colnames(avg) <- NULL
pdf("Fig3e2.pdf", height=8.27, width=8.27)
pheatmap(avg, color = mapal, cluster_cols = F, cluster_rows = F, border_color = NA)
dev.off()


## Fig3f 
fig3f <- monocle::plot_cell_trajectory(cds, cell_size=0.5, color_by="group") +
  scale_color_manual(values=groupCols)
fig3f


## fig3g
# see Velocity.ipynb

