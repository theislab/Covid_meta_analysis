seurat()
setwd("/broad/hptmp/subraman/covid19")
library(dplyr)
library(lme4)
library(car)
library(Hmisc)

source("~/code/plot/network.r")
  
mast_gut=readRDS("/broad/hptmp/csmillie/gut.ACE2_TMPRSS2.rds")
mast_lca=readRDS("/broad/hptmp/csmillie/lca.ACE2_TMPRSS2.rds")
mast_lung=readRDS("/broad/hptmp/csmillie/lung.ACE2_TMPRSS2.rds")
mast_nose=readRDS("/broad/hptmp/csmillie/nose.ACE2_TMPRSS2.rds")

##### Figure 4c #####

a=12
genelist=NULL
Tissue="Lung"
templist=NULL
  genes_filter=mast_lung$AT2$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=data.frame(tissue=Tissue,celltype="AT2",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[c(1:a),])
  templist=NULL
  genes_filter=mast_lung$AT1$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=data.frame(tissue=Tissue,celltype="AT1",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[1:a,])
  templist=NULL
  genes_filter=mast_lung$Secretory$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=data.frame(tissue=Tissue,celltype="Secretory",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[1:a,])
  templist=NULL
  genes_filter=mast_lung$Basal$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=data.frame(tissue=Tissue,celltype="Basal",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[c(1:a,28,39),])
  templist=NULL
  genes_filter=mast_lung$Ciliated$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=data.frame(tissue=Tissue,celltype="Ciliated",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[1:a,])
  lung_genelist=genelist
  
  df=lung_genelist[,c(2:3)]
  names(df)=c("source","target")
  network=data.frame(df,weight=0.25)

  p=nice_network(edges=network,node_pal = set.colors)
  pdf("Fig4c.pdf",10,10,useDingbats=FALSE)
  print(p)
  dev.off()

  
  
##### Figure 12a #####
  
  a=12
  
  genelist=NULL
  Tissue="LCA"

   genes_filter=mast_lca$AT2$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=NULL
  templist=data.frame(tissue=Tissue,celltype="AT2",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[1:a,])
  genes_filter=mast_lca$Secretory$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=NULL
  templist=data.frame(tissue=Tissue,celltype="Secretory",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[1:a,])
  genes_filter=mast_lca$`Ciliated lineage`$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=NULL
  templist=data.frame(tissue=Tissue,celltype="Ciliated lineage",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[1:a,])
  lca_genelist=genelist
  
  df=lca_genelist[,c(2:3)]
  names(df)=c("source","target")
  network=data.frame(df,weight=0.25)
  
  p=nice_network(edges=network,node_pal = set.colors)
  pdf("FigS12a.pdf",10,10,useDingbats=FALSE)
  print(p)
  dev.off()
  
  ##### Figure 12b #####  

a=12
  genelist=NULL
  Tissue="Gut"
  genes_filter=mast_gut$Enterocytes$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=NULL
  templist=data.frame(tissue=Tissue,celltype="Enterocytes",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[c(1:a,20),])
  genes_filter=mast_gut$`Immature Enterocytes 2`$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=NULL
  templist=data.frame(tissue=Tissue,celltype="Immature Enterocytes 2",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[c(1:a,51),])
  genes_filter=mast_gut$`Cycling TA`$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=NULL
  templist=data.frame(tissue=Tissue,celltype="Cycling TA",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[c(1:8,10,12:14,35,36),])
  genes_filter=mast_gut$`TA 2`$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=NULL
  templist=data.frame(tissue=Tissue,celltype="TA 2",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[c(1:6,8:11,13,63,64),])
  genes_filter=mast_gut$`Best4+ Enterocytes`$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=NULL
  templist=data.frame(tissue=Tissue,celltype="Best4+ Enterocytes",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[c(1:5,7:13,37,66),])

  
  gut_genelist=genelist
  
  df=gut_genelist[,c(2:3)]
  names(df)=c("source","target")
  network=data.frame(df,weight=0.25)
  
  p=nice_network(edges=network,node_pal = set.colors)
  pdf("FigS12b.pdf",10,10,useDingbats=FALSE)
  print(p)
  dev.off()
  
  
  
  ##### Figure 4d ##### 
  
  a=12
  
  genelist=NULL
  Tissue="Nose"
  genes_filter=mast_nose$Ciliated$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=NULL
  templist=data.frame(tissue=Tissue,celltype="Ciliated",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[1:a,])
  genes_filter=mast_nose$Goblet$mast%>%subset(padjH<0.05)%>%subset(mastfc>0)%>%subset(alpha>0.1)%>%arrange(desc(mastfc))
  templist=NULL
  templist=data.frame(tissue=Tissue,celltype="Goblet",gene=genes_filter$gene,padjH=genes_filter$padjH,padjH=genes_filter$padjH,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
  genelist=rbind(genelist,templist[1:a,])
  
  nose_genelist=genelist
  
  df=nose_genelist[,c(2:3)]
  names(df)=c("source","target")
  network=data.frame(df,weight=0.25)
  
  p=nice_network(edges=network,node_pal = set.colors)
  pdf("Fig4d.pdf",10,10,useDingbats=FALSE)
  print(p)
  dev.off()
  
  
  