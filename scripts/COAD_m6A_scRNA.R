#01.#########
####
GSE146771_cli <- read.delim('origin_datas/GSE146771/GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt.gz')
GSE146771_count <- fread('origin_datas/GSE146771/GSE146771_CRC.Leukocyte.Smart-seq2.TPM.txt.gz',data.table = F)
class(GSE146771_count)
GSE146771_count[1:5,1:5]
rownames(GSE146771_count) <- GSE146771_count$V1
GSE146771_count <- GSE146771_count[, -1]
GSE146771_count[1:5, 1:5]
dim(GSE146771_count)
length(intersect(GSE146771_cli$CellName, colnames(GSE146771_count)))
head(GSE146771_cli)
table(GSE146771_cli$Sample,GSE146771_cli$Tissue)
table(GSE146771_cli$Tissue)
# GSE146771_cli=subset(GSE146771_cli,Tissue=='T')
dim(GSE146771_cli)
GSE146771_count=GSE146771_count[,GSE146771_cli$CellName]
dim(GSE146771_count)

datalist <- list()
for (i in unique(GSE146771_cli$Sample)) {
  ce <- GSE146771_cli$CellName[GSE146771_cli$Sample == i]
  my.count <- GSE146771_count[, ce]
  datalist[[i]] <- CreateSeuratObject(counts=my.count,project = i, min.cells = 3,min.features = 200)
  datalist[[i]]$Samples=GSE146771_cli$Sample[GSE146771_cli$Sample == i]
  datalist[[i]]$tissue=GSE146771_cli$Tissue[GSE146771_cli$Sample == i]
  datalist[[i]]$Global_Cluster=GSE146771_cli$Global_Cluster[GSE146771_cli$Sample == i]
  rm(my.count)
}
# names(datalist)=unique(GSE146771_cli$Sample)
# rm(GSE146771_count)

dir.create('results/06.scNRA')
##1.1 ######
####
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 
  datalist[[i]] <- sce
  rm(sce)
}
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
#
raw_meta=sce@meta.data
raw_count <- table(raw_meta$Samples)
raw_count
sum(raw_count)#  10468
pearplot_befor<-VlnPlot(sce,group.by ='Samples',
                        features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        pt.size = 0,
                        ncol = 3)
pearplot_befor
dev.off()
ggsave('results/06.scNRA/pearplot_befor.pdf',pearplot_befor,height = 5,width = 10)

# #
sce=subset(sce, subset=nFeature_RNA>1000 & nFeature_RNA<7000 & percent.mt<5)
# rm(datalist)
clean_meta=sce@meta.data
clean_count <- table(clean_meta$Samples)
clean_count
sum(clean_count)#10186
pearplot_after <- VlnPlot(sce,group.by ='Samples',
                          features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                          pt.size = 0,
                          ncol = 3)
pearplot_after
ggsave('results/06.scNRA/pearplot_after.pdf',pearplot_after,height = 5,width = 10)



##1.2 ####
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce, features = rownames(sce))
#
sce <- RunPCA(sce, features = VariableFeatures(sce))
colnames(sce@meta.data)
##
library(harmony)
sce = RunHarmony(sce, group.by.vars="Samples")
#
pca.plot=ElbowPlot(sce,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))
pca.plot
# ggsave('results/01.scRNA/PCA.pdf',pca.plot,height = 5,width = 5)
###
sce <- RunTSNE(sce, dims=1:20, reduction="harmony")
####
after_batch=DimPlot(sce,group.by='Samples',reduction="tsne",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
after_batch
ggsave('results/06.scNRA/after_batch.pdf',after_batch,height = 6,width = 6.5)


library(clustree)
sce <- FindNeighbors(sce, dims = 1:20, reduction="harmony")
#
sce <- FindClusters(object = sce,resolution = .1)
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
length(table(sce@meta.data$seurat_clusters))


DimPlot(sce,group.by='tissue',reduction="tsne",label = F,pt.size = 0.2,cols =pal_d3()(8))
# DimPlot(sce,group.by='Global_Cluster',reduction="umap",label = F,pt.size = 0.2,cols =pal_d3()(8))

p=DimPlot(sce,group.by='seurat_clusters',reduction="tsne",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
p=LabelClusters(p,id = 'seurat_clusters',family='Times')
p

ggsave('results/06.scNRA/cluster_umap.pdf',p,height = 5,width = 5.5)


##1.3 #####
#
Logfc = 0.5
#
Minpct = 0.25
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
Idents(sce)<-'seurat_clusters'

sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
head(sce.markers)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
table(sce.markers$cluster)
length(unique(sce.markers$gene))#762
head(sce.markers)
table(sce.markers$cluster)
write.csv(sce.markers,'results/06.scNRA/diff_marker_gene.csv')
### 
Top5 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_logFC)
length(Top5$gene)
length(unique(Top5$gene))
###
diff.marker.dotplot=DotPlot(object = sce, features = unique(Top5$gene),
                            cols=c("snow", "blue"),scale = T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()
diff.marker.dotplot

ggsave('results/06.scNRA/cluster_diffmarker.pdf',diff.marker.dotplot,height = 10,width = 10)





##1.4 #####
marker <- data.frame(cluster = 0:11,cell = 0:11)
marker[marker$cluster %in% c(0,10),2] <- 'CD8+ T cells'
marker[marker$cluster %in% c(1),2] <- 'NKT cells'
marker[marker$cluster %in% c(2),2] <- 'CD4+ cells'
marker[marker$cluster %in% c(3),2] <- 'Epithelial cells'
marker[marker$cluster %in% c(4),2] <- 'Monocyte'
marker[marker$cluster %in% c(5),2] <- 'Tregs'
marker[marker$cluster %in% c(6),2] <- 'Neutrophil'
marker[marker$cluster %in% c(7),2] <- 'B cells 1'
marker[marker$cluster %in% c(8),2] <- 'B cells 2'
marker[marker$cluster %in% c(9),2] <- 'Mast cells'
marker[marker$cluster %in% c(11),2] <- 'Fibroblasts'

marker
sce@meta.data$cell_type <- sapply(sce@meta.data$seurat_clusters,function(x){marker[x,2]})
my.cols=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#FFED6F",'#66A61E')
cell_type_umap=DimPlot(sce,group.by='cell_type',reduction="tsne",label = F,pt.size = 0.2,cols = my.cols)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
cell_type_umap=LabelClusters(cell_type_umap,id = 'cell_type',family='Times')
cell_type_umap
ggsave('results/06.scNRA/cell_type_umap.pdf',cell_type_umap,height = 5,width = 5.5)
table(sce@meta.data$cell_type)



# FeaturePlot(tumor.sce, features = c("MAL", "LEF1"))
# FeaturePlot(tumor.sce, features = c("KRT19", "EPCAM"))
# FeaturePlot(tumor.sce, features = c("IGFBP7", "COL4A1"))
# FeaturePlot(tumor.sce, features = c('CD8A','GZMK'))
# FeaturePlot(tumor.sce, features = c('S100A9','LYZ'))
# FeaturePlot(tumor.sce, features = c('TPSAB1','CPA3'))
# FeaturePlot(tumor.sce, features = c('FGFBP2','KLRF1'))
# FeaturePlot(tumor.sce, features = c('CD79A'))



marker_gene=c('CD79A','MS4A1','JCHAIN',
              'LEF1','MAL',
              'CD8A','GZMK','CD8B',
              'KRT19','KRT18','EPCAM',
              'IGFBP7','COL4A1',
              'TPSAB1','CPA3',
              'S100A9','LYZ','CD14','G0S2',
              'FGFBP2','KLRF1',
              'CTLA4')
# Idents(sce)='cell_type'
# Idents(sce)<-'seurat_clusters'
marker.dot=DotPlot(object = sce, features = marker_gene,group.by = 'cell_type',
                   cols=c("#CCEBC5", "#D95F02"),scale = T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size=12),text = element_text(family = 'Times',size = 12)) +
  xlab('')+ylab('')+coord_flip()
marker.dot
ggsave('results/06.scNRA/marker_dotplot.pdf',marker.dot,height = 6,width = 7)



saveRDS(sce,file = 'results/06.scNRA/sce.rds')



##1.5 ######
cell_freq1=data.frame(t(prop.table(table(tumor.sce$cell_type,tumor.sce$Samples),margin=2)))
cell_freq1
colnames(cell_freq1)<-c('Samples','cell_type','Freq')
cell_prop_fig1=ggplot(cell_freq1,aes(x=Samples,y=Freq,fill=cell_type))+
  scale_fill_manual(values = my.cols)+
  # facet_grid(~State,scales = 'free',space='free')+
  geom_bar(position = "fill",stat="identity")+
  xlab('')+ylab('Proportion')+theme_bw()+
  theme(text = element_text(family = 'Times',size=15),
        # axis.text.x = element_text(angle = 30,hjust = 1),
        legend.text =element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))
cell_prop_fig1

model.dotplot=DotPlot(object = sce, features =c("METTL3","YTHDC2","IGF2BP3"),group.by = 'cell_type',
        cols=c("#CCEBC5", "#D95F02"),scale = T,col.min = 0)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size=12),text = element_text(family = 'Times',size = 12)) +  xlab('')+ylab('')+coord_flip()
model.dotplot

#########
table(sce$tissue,sce$Samples)
tumor.sce=subset(sce,tissue=='T')

countexp<-tumor.sce@assays$RNA@counts
countexp<-data.frame(as.matrix(countexp))
library(AUCell)
library(GSEABase)
cells_rankings <- AUCell_buildRankings(as.matrix(countexp)) #rank

m6A.regulators=read.xlsx('origin_datas/m6A_regulators_PMID_35646049.xlsx')
m6A.regulators=m6A.regulators$Regulator
length(m6A.regulators)
geneSets=list(m6A=m6A.regulators)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #cal

m6A_AUC <- as.numeric(getAUC(cells_AUC)['m6A', ])
tumor.sce$m6A_AUC  <- m6A_AUC

sc_df=data.frame(tumor.sce@meta.data, tumor.sce@reductions$tsne@cell.embeddings)
head(sc_df)

# auc_stat <- sc_df  %>%
#   wilcox_test(m6A_AUC ~ cell_type, paired = FALSE, p.adjust.method = "fdr") %>%
#   add_significance() %>% add_xy_position(x = "cell_type")
# head(auc_stat)

# pdf('results/04.PD-L1/TCGA_PDL1_exp.pdf',height = 5,width = 6.5,onefile = F)
m6A.score.boxplot=ggplot(sc_df, aes(x = cell_type, y = m6A_AUC))+
  geom_boxplot(aes(fill = cell_type), varwidth = T, alpha = 0.9, outlier.shape = NA)+
  scale_fill_manual(values=my.cols)+
  stat_compare_means(aes(group=cell_type), label = "p.format", method = 'kruskal.test')+
  theme_bw(base_size = 20)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                 axis.text.x = element_text(angle = 30,hjust = 1),legend.position = 'none',
                                 text = element_text(family ='Times'))+
  labs(title= "AUCell" ,y="m6A score",x='')
  # stat_pvalue_manual(auc_stat, label = "p.adj.signif", hide.ns = T, tip.length = 0.008)
m6A.score.boxplot


pdf('results/06.scNRA/Fig6.pdf',height = 10,width = 12)
mg_merge_plot(cell_type_umap,marker.dot,model.dotplot,m6A.score.boxplot,nrow=2,ncol=2,heights = c(1.2,1))
dev.off()  
