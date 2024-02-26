setwd("F:/学习/数据分析/single cell sequence/测序数据/human/衰老分析NEW")
save.image(file = "衰老分析.Rdata")

library(cowplot)
library(Seurat)
library(dplyr)
library(monocle)
library(stringr)
library(pheatmap)
library(vegan)
library(biomaRt)
library(RColorBrewer)
require(scales)
require(grid)
library(plyr)
library(vegan)  
require(reshape2)
require(ggsci)
library(monocle3)

#查看调色板 display.brewer.all() display.brewer.all(type = "div")

color_plot<-c("#66C2A5","#3288BD", "#5E4FA2") #三衰老配色
color_plot<-c("#fff154", "#66C2A5","#3288BD", "#5E4FA2") #四衰老配色


#批次效应处理
allsample_combined_cca<-subset(allsample_combined_cca,features=rownames(lz005@assays$RNA@counts))
allsample_combined_cca_list <- SplitObject(allsample_combined_cca, split.by = "sample")
for (i in 1:length(allsample_combined_cca_list)) {
  allsample_combined_cca_list[[i]] <- NormalizeData(allsample_combined_cca_list[[i]], verbose = T)
  allsample_combined_cca_list[[i]] <- FindVariableFeatures(allsample_combined_cca_list[[i]], selection.method = "vst", 
                                                           nfeatures = 1000, verbose = T)
}

allsample_combined_cca_anchors <- FindIntegrationAnchors(object.list = allsample_combined_cca_list, dims = 1:30,anchor.features = 1000)
allsample_combined_cca <- IntegrateData(anchorset = allsample_combined_cca_anchors, dims = 1:30)


#质量控制
allsample_combined_cca <- subset(x = allsample_combined_cca, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 30 & nCount_RNA <100000)
FeaturePlot(allsample_combined_cca, features = c("nFeature_RNA","nCount_RNA" ,"percent.mt"), cols = c("blue", "yellow"),min.cutoff = "q9",pt.size = 0.1, reduction = "umap")
VlnPlot(allsample_combined_cca,features = c("nFeature_RNA","nCount_RNA" ,"percent.mt"),pt.size = 0.1)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )

#标准化
#DefaultAssay(allsample_combined_cca) <- "RNA"
#DefaultAssay(allsample_combined_cca) <- "integrated"
#allsample_combined_cca<-NormalizeData(allsample_combined_cca)
#allsample_combined_cca<-FindVariableFeatures(allsample_combined_cca, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
allsample_combined_cca <- ScaleData(allsample_combined_cca)
allsample_combined_cca <- RunPCA(allsample_combined_cca, npcs = 30)
#聚类
allsample_combined_cca <- FindNeighbors(object = allsample_combined_cca, dims = 1:30)
allsample_combined_cca <- FindClusters(object = allsample_combined_cca, resolution = 0.5)
#降维
for (i in 1:20) {
  allsample_combined_cca <- RunUMAP(object = allsample_combined_cca, dims = 1:30,n.neighbors = 40L,seed=i)
  a<-DimPlot(object = allsample_combined_cca, reduction = "umap",  label = T,label.size = 8,pt.size = 0.1)+labs(title = i)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())
  print(a)
}
allsample_combined_cca <- RunTSNE(object = allsample_combined_cca, dims = 1:30,n.neighbors = 35L,seed=17)

# Visualization
DimPlot(object = allsample_combined_cca, reduction = "umap",  label = T,label.size = 8,pt.size = 0.1)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())

#提取细胞高亮显示
Idents(allsample_combined_cca)<-allsample_combined_cca$sample
allsample_combined_cca_subset<-subset(x = allsample_combined_cca, idents = c("LZ007"))
a<-rownames(allsample_combined_cca_subset@meta.data)
DimPlot(object = allsample_combined_cca, reduction = "umap",cells.highlight = a,label.size = 8,pt.size = 0.1)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())


#看细胞亚型基因表达
Idents(allsample_combined_cca)<-allsample_combined_cca$seurat_10clusters
allsample_combined_cca_subset<-subset(x = allsample_combined_cca, idents = "SC")#亚型
Idents(allsample_combined_cca_subset)<-allsample_combined_cca_subset$age_4C #分组
VlnPlot(allsample_combined_cca_subset,features = c("JUN","EGR3","NR4A1"),pt.size = 0.1,group.by = "age_4C",cols = color_plot_qualitative)+NoLegend()+geom_smooth(aes(group=1),size=2,method = 'lm',formula = y~ns(x,3),color="black")+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )
FeaturePlot(object = allsample_combined_cca, reduction = "umap", features = c("SOD3"), cols = c("#eeeeee", "darkred"), ncol = 1,pt.size = 0.1)+theme(axis.title=element_blank(),axis.text = element_text(size = 0),legend.position = c(.85,.15))

VlnPlot(allsample_combined_cca_10X_KSandOA_germcell,features = c("VWF"),pt.size = 0.1,group.by = "seurat_GERMclusters",split.by = "age")+scale_fill_manual(values =mypalette)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 30),plot.title =element_text(size = 45,face="plain") )
VlnPlot(allsample_combined_cca,features = c("WISP2"),pt.size = 0.1,group.by = "seurat_12clusters",cols = cl)+geom_boxplot(fill = "black",color="gray90",size=1,width=0.1,outlier.size = 0)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )

#阳性细胞数较少时用dotplot看表达
DotPlot(allsample_combined_cca_subset,features = rev(c("JUN","EGR3","NR4A1","S100A13","ENO1","BEX1","HOPX","DEFB119","CST9L","TUBA3C","RBP7","CTSF")),group.by = "age4")+scale_size_continuous(range = c(1,8))+coord_flip()
DotPlot(allsample_combined_cca1,features = rev(c("JUN","EGR3","NR4A1","S100A13","ENO1","BEX1","HOPX","DEFB119","CST9L","TUBA3C","RBP7","CTSF")),group.by = "seurat_12clusters",split.by = "NULL",cols=c(rep("blue",12), "white"))+scale_size_continuous(range = c(1,8))+coord_flip()

#带趋势线
VlnPlot(allsample_combined_cca_LCPTM,features = c("gene_SASP1"),pt.size = 0.1)+geom_boxplot(fill = "black",color="gray90",size=1,width=0.05,outlier.size = 0)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )+stat_summary(aes(group=1),fun.y=mean, geom="smooth", shape=1,size=2,color="gray")+ scale_y_continuous(limits = c(0,2))

#计算基因集（多个基因）表达
listfori<-c("gene_SASP","gene_collagen","gene_Interleukins","gene_Chemokines","gene_Otherinflammatorymolecules","gene_Growthfactors","gene_Proteases","gene_L_R","gene_PEGNOS","gene_Insolublefactors")
for (i in listfori) {
  allsample_combined_cca<-AddModuleScore(allsample_combined_cca,features = get(i),name = i)
}


#建立表格数据
Idents(allsample_combined_cca_adult)<-allsample_combined_cca_adult$seurat_10clusters
for(j in c("SPG" ,"SPC","SPT","SC","LC","PTM","EC","VSM","MAC","LYM")){
  allsample_combined_cca_subset<-subset(x = allsample_combined_cca_adult, idents = j)
  allsample_combined_cca_subset<-NormalizeData(allsample_combined_cca_subset)
  gene_subset<-row.names(subset(allsample_combined_cca_subset@assays$RNA@meta.features,vst.mean>0.1))
  MARKERS<-t(data.frame(allsample_combined_cca_subset@assays$RNA@data))[,gene_subset]
  MARKER<-cbind(allsample_combined_cca_subset@meta.data[,c(2,3,19:35)],MARKERS) #在回归分析中加入的额外参数，需要根据meta.data修改
  MARKER$age3<-as.numeric(MARKER$age3) #将年龄等级数据转换为数值
  a<-colnames(MARKER)
  #回归分析（每次做前都建立空表格）
  c<-read.csv("空表格.csv",row.names=1) #建立空表格
  for(i in a){
    b<-lm(get(i)~age3,MARKER) #填入相应的Y~X
    b<-summary(b)
    b1<-t(data.frame(b$coefficients[2,])) #导出相关系数
    b1$R_squared<-data.frame(b$r.squared) #导出r值
    b1<-data.frame(b1)
    colnames(b1)<-c("a","Std.Error","t","p","r.squared")
    rownames(b1)<-i
    c<-rbind(c,b1)
  }
  write.csv(c,paste0("所有基因随衰老变化回归分析",j,".csv"))
}
#作图
#回归图
ggplot(data = MARKER, mapping = aes(x = age3, y = RPS26),geom="jitter")+ 
  geom_point(aes(color = nFeature_RNA,size= nCount_RNA),alpha = 1/10)+ #以drv为分组设置点的颜色
  geom_smooth(method = 'lm', formula = y ~ x,colour="ORANGE",size=2)+ scale_y_continuous(limits = c(1,5))
#


#计算比例
allsample_combined_cca[["percent.X_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= gene_X)
allsample_combined_cca[["percent.Y_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= gene_Y)

#画各个样本比例
ggplot(data = allsample_combined_cca@meta.data, mapping = aes(x = age2, fill = seurat_12clusters))+geom_bar(stat = 'count', position = 'fill')+scale_fill_manual(values=color_plot) + labs(x = '', y = 'Rate')+theme_classic()


#基因表达相关性回归分析
MARKER<-allsample_combined_cca_10X_KS_SC@meta.data
MARKERS<-allsample_combined_cca_10X_KS_SC@assays$RNA@data
MARKERS<-data.frame(MARKERS)
MARKERS<-t(MARKERS)
MARKER<-cbind(MARKER,MARKERS)

ggplot(data = MARKER, mapping = aes(x = XIST, y = percent.chr_neXi_gene))+ 
  geom_point(aes(color = percent.chr_X_gene,size= nCount_RNA))+ #以drv为分组设置点的颜色
  geom_smooth(method = 'lm', formula = y ~ x,colour="ORANGE",size=2)+labs(caption ="y = 4.092-0.698x  p-value = 8.185e-14")
MARKERS<-lm(percent.chr_eXi_gene~XIST,MARKER)
summary(MARKERS) 

#绘制火山图
MARKER<-read.csv("KS SSC 所有基因与all X回归分析.csv",header = T,stringsAsFactors=F)

ggplot(data = MARKER, aes(x = r, y = -log10(p.value),size=mean,color=UPDOWN)) +
  geom_point(alpha=0.6)+scale_size(limits=c(0,4)) + scale_color_manual(values=c("#90bff9", "grey","#f2b77c"))+
  xlim(c(-1.1, 0.7)) +
  ylim(c(0, 4.7)) +
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="gray",lwd=0.8)


#各组平均表达热图
cgene<-c("GFRA1","SYCP3","TNP1","SOX9","STAR","MYH11","VWF","NOTCH3","CD163","TPSAB1","CD3D")
MARKERS <- AverageExpression(allsample_combined_cca, return.seurat = T)
MARKERS<-ScaleData(MARKERS)
DoHeatmap(MARKERS, features = c("USP9Y","DDX3Y","UTY","HSFY1","HSFY2","RBMY1A1","RBMY1B","RBMY1C","RBMY1D","RBMY1E","RBMY1J","DAZ1","DAZ2","DAZ3","DAZ4","BPY2","CDY1B","PRY","CSPG4P1Y","DAZL","BOLL","SRY","ZFX","ZFY"), size = 3,slot = "data", disp.min = 0, disp.max = 5)+scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkorange")
DoHeatmap(MARKERS, features = c("RBMY1A1","RBMY1F","RBMY1J","RBMY1B","RBMY1D","RBMY1E"), size = 6,draw.lines = F,group.colors = cl)+scale_fill_gradient2(low = "#02197b", mid = "white", high = "#b10000")

#计算多个基因打分
IMMUNE_MARKER<-list(c("PTPRC","TPSAB1","CD163","CD3D","CD4","CD8A","MS4A1", "CD79A","GNLY","ITGAX"))
gene_androgen_synthesis<-list(c("STAR","CYP11A1","CYP17A1","HSD17B6","SCARB1","MED1","HSD3B2", "HSD17B3","SRD5A1","SRD5A2","SRD5A3"))
allsample_combined_cca<-AddModuleScore(object = allsample_combined_cca,features = gene_paternally_imprinted,name = 'gene_paternally_imprinted')

#计算差异度
#建立表格数据
Idents(allsample_combined_cca)<-allsample_combined_cca$age2
allsample_combined_cca_subset<-subset(x = allsample_combined_cca, idents = c("childbearing","agedness"))
MARKER<-AverageExpression(allsample_combined_cca_subset,add.ident = "seurat_10clusters")
MARKER<-MARKER$integrated
#MARKER<-MARKER$RNA
#MARKER<-MARKER[allsample_combined_cca_8yo@assays$RNA@var.features,]或者取前1000个差异基因
MARKER<-t(MARKER)
MARKER<-vegdist(MARKER,'bray')
MARKER<-as.matrix(MARKER)
write.csv(MARKER,"各年龄各细胞bray差异度.csv") #手动修改成作图的两列
MARKER<-read.csv("细胞差异度10cluster.csv",,header = T)
#标准化各组
MARKER$cell<-factor(MARKER$cell,levels = rev(c("SPG","SPC","SPT","SC","LC","PTM","EC","VSM","MAC","LYM")))
ggplot(data = MARKER, mapping = aes(x=Dis, y = cell))+ 
  geom_point(aes(color = Dissimilarity.Jaccard,size= Dissimilarity.Bary))+scale_color_gradient(low = "gray", high = "red")+
  scale_size_continuous(range=c(4,10))+theme(axis.title=element_blank(),axis.text = element_text(size = 20))

#小提琴图
Idents(allsample_combined_cca)<-allsample_combined_cca$seurat_10clusters
allsample_combined_cca_subset<-subset(x = allsample_combined_cca, idents = "LC")
Idents(allsample_combined_cca_subset)<-allsample_combined_cca_subset$age2
VlnPlot(allsample_combined_cca_subset,features = c("gene_triglyceride_metabolic1"),pt.size = 0.1,group.by = "age2",cols = color_plot)+NoLegend()+geom_smooth(aes(group=1),size=2,method = 'lm',formula = y~ns(x,8),color="black")+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )

allsample_combined_cca_subset <- subset(x = allsample_combined_cca_subset, subset = percent.Y_gene <0.3)

#画带透明度的Featureplot
i<-c("GFRA1")
P1<-FeaturePlot(object = allsample_combined_cca, reduction = "umap", features = c(i), ncol = 1,pt.size = 0.1)+theme(axis.title=element_blank(),axis.text = element_text(size = 0),legend.position = c(.85,.15))
ggplot(P1$data,aes(UMAP_1,UMAP_2,color=get(i)))+geom_point(size=0.1,alpha=0.7)+scale_color_gradient(low = "#f3f7f7",high = "#5c0000")+P1$theme+labs(title=i)+theme(plot.title = element_text(hjust = 0.5))

#转录因子调控关系
allsample_combined_cca_subset<-subset(x = allsample_combined_cca_sertolicells, idents = c("childbearing","agedness"))
MARKERS<-read.csv("F:/学习/数据分析/single cell sequence/测序数据/human/衰老分析/SC/OA OLD DEG.csv")
exprMatr<-as.matrix(allsample_combined_cca_subset[["RNA"]]@data)[MARKERS$gene,]
head(exprMatr[1:5,1:5])
regulators <- read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/Homo_sapiens_TF.csv",header = F,stringsAsFactors=F)[,1]
regulators<-intersect(regulators,unique(rownames(MARKERS)))
weightMat <- GENIE3(exprMatr, regulators=regulators, nCores=4, verbose=TRUE)
linkList <- getLinkList(weightMat)
linkList <- getLinkList(weightMat, reportMax=5)
linkList <- getLinkList(weightMat, threshold=0.1)

#手动气泡图
MARKER<-read.csv("F:/学习/数据分析/single cell sequence/测序数据/human/ED海绵体分析/FIGURE/FIG 2 INC/CELL PROLIFERATION AND DEATH.csv",header = T,stringsAsFactors=F)

MARKER$Diseases.or.Functions.Annotation<-factor(MARKER$Diseases.or.Functions.Annotation,levels = rev(c("Cell cycle progression","Cell proliferation of fibroblasts","Proliferation of connective tissue cells","Proliferation of epithelial cells","Proliferation of smooth muscle cells","Cell viability","Cell survival","Apoptosis","Cell death of connective tissue cells","Cell death of epithelial cells","Cell death of muscle cells")),ordered = F)

ggplot(MARKER,aes(cluster,Diseases.or.Functions.Annotation,p.Value,Bias.corrected.z.score))+ 
  geom_point(aes(size=-1*log10(p.Value),color=Bias.corrected.z.score))+
  scale_colour_gradient2(low="lightblue",mid="lightgrey",high="red")+
  scale_size_continuous(range=c(4,10))+
  scale_y_discrete(position = "right")+ theme_bw()+theme(panel.grid.major=element_line(colour=NA))

#基因表达量和占比气泡图
DotPlot(allsample_combined_cca1,features = rev(c("JUN","EGR3","NR4A1","S100A13","ENO1","BEX1","HOPX","DEFB119","CST9L","TUBA3C","RBP7","CTSF")),group.by = "seurat_12clusters",split.by = "NULL",cols=c(rep("blue",12), "white"))+scale_size_continuous(range = c(1,8))+coord_flip()+theme(axis.text.y=element_text(face=c("italic")))

#随机抽取细胞M","EC","VSM","MAC")) {
  for(j in c("childbearing", "midlife","old"))
    allsample_combined_cca_subset<-subset(x = allsample_combined_cca, idents = i)
  Idents(allsample_combined_cca_subset)<-allsample_combined_cca_subset$age2
  allsample_combined_cca_subset<-subset(x = allsample_combined_cca_subset, idents = j)
  cellname<-rownames(allsample_combined_cca_subset@meta.data)
  cellname<-sample(cellname,100)
}
allsample_combined_cca_subset<-subset(x = allsample_combined_cca, idents = c("LC"))

sample(,100)

