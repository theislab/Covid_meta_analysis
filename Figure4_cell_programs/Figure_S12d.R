mast_gut=readRDS("gut.ACE2_TMPRSS2.rds")
mast_lca=readRDS("lca.ACE2_TMPRSS2.rds")
mast_lung=readRDS("lung.ACE2_TMPRSS2.rds")
mast_nose=readRDS("nose.ACE2_TMPRSS2.rds")


##### Generate table for plotting ########
genelist=NULL
Tissue="Lung"
genes_filter=mast_lung$AT2$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="AT2",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
genes_filter=mast_lung$AT1$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="AT1",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
genes_filter=mast_lung$Secretory$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="Secretory",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
genes_filter=mast_lung$Basal$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="Basal",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
genes_filter=mast_lung$Ciliated$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="Ciliated",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
lung_genelist=genelist

genelist=NULL
Tissue="LCA"
genes_filter=mast_lca$AT2$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="AT2",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
genes_filter=mast_lca$Secretory$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="Secretory",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
genes_filter=mast_lca$`Ciliated lineage`$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="Ciliated lineage",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
lca_genelist=genelist

genelist=NULL
Tissue="Gut"
genes_filter=mast_gut$Enterocytes$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="Enterocytes",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
genes_filter=mast_gut$`Immature Enterocytes 2`$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="Immature Enterocytes 2",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
genes_filter=mast_gut$`Cycling TA`$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="Cycling TA",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
genes_filter=mast_gut$`TA 2`$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="TA 2",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
genes_filter=mast_gut$`Best4+ Enterocytes`$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=NULL
templist=data.frame(tissue=Tissue,celltype="Best4+ Enterocytes",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
gut_genelist=genelist


genelist=NULL
Tissue="Nose"
genes_filter=mast_nose$Ciliated$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=data.frame(tissue=Tissue,celltype="Ciliated",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)
genes_filter=mast_nose$Goblet$mast%>%subset(gene%in%c("IL6","IL6R","IL6ST"))
templist=data.frame(tissue=Tissue,celltype="Goblet",gene=genes_filter$gene,padjH=genes_filter$padjH,padjD=genes_filter$padjD,log2fc=genes_filter$log2fc,alpha=genes_filter$alpha,refalpha=genes_filter$ref_alpha)
genelist=rbind(genelist,templist)

nose_genelist=genelist

genelist=rbind(lung_genelist,gut_genelist,nose_genelist,lca_genelist)

######### PLOT ###########
library(tidyr)
library(ggpubr)


df=genelist%>%subset(alpha>0)
df$padjH[is.na(df$padjH)] = 1
df$celltype=as.character(df$celltype)
df$celltype=paste0(df$tissue,"_",df$celltype)

df$celltype=factor(df$celltype,levels=unique(df$celltype))
df=df%>%mutate(neglog_pval_adj=-log10(padjH))
df=df%>%mutate(Significant=ifelse(padjH<0.05,"True","False"))
library(ggpubr)

g <- ggballoonplot(df, 
                   x = 'gene',
                   y = 'celltype',
                   size = 'neglog_pval_adj',
                   color = "Significant",
                   fill = 'log2fc', 
                   font.xtickslab=9, 
                   size.range=c(1,5), 
                   main='ACE2+/TMPRSS2+ DP vs DN',
                   rotate.x.text=T) + 
  ggplot2::scale_fill_distiller(type='div', limit=limit, name='log2 fold-change')  + ggplot2::scale_color_manual(values=c("#808080","#990E1D"))



pdf(paste0("Fig.pdf"))
print(g)
dev.off()
