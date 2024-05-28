library(fgsea)
library(ggplot2)
library(msigdbr)
library(plyr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(stringr)
##################################################
#####------genes with significant DTU---------####
##################################################
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/scotch")
scotch_results_B_8=read.csv('scotch_BvsAll_sample8.csv')
scotch_results_Monocyte_8=read.csv('scotch_MonocytesvsAll_sample8.csv')
scotch_results_CD4_8=read.csv('scotch_CD4vsAll_sample8.csv')
scotch_results_CD8_8=read.csv('scotch_CD8vsAll_sample8.csv')
scotch_results_NK_8=read.csv('scotch_NKvsAll_sample8.csv')

scotch_results_B_7=read.csv('scotch_BvsAll_sample7.csv')
scotch_results_Monocyte_7=read.csv('scotch_MonocytesvsAll_sample7.csv')
scotch_results_CD4_7=read.csv('scotch_CD4vsAll_sample7.csv')
scotch_results_CD8_7=read.csv('scotch_CD8vsAll_sample7.csv')
scotch_results_NK_7=read.csv('scotch_NKvsAll_sample7.csv')

common_genes_8 = Reduce(intersect, list(scotch_results_B_8$genes, 
                                        scotch_results_Monocyte_8$genes, 
                                        scotch_results_CD4_8$genes, 
                                        scotch_results_CD8_8$genes, 
                                        scotch_results_NK_8$genes))
common_genes_7 = Reduce(intersect, list(scotch_results_B_7$genes, 
                                        scotch_results_Monocyte_7$genes, 
                                        scotch_results_CD4_7$genes, 
                                        scotch_results_CD8_7$genes, 
                                        scotch_results_NK_7$genes))
gene_universe=intersect(common_genes_8, common_genes_7)
#################################################################
#####------b. upset plot: genes with significant DTU---------####
#################################################################
#monocyte
a=unique(scotch_results_Monocyte_8%>%filter(p_DTU_gene_adj<=0.05)%>%pull(genes))
b=unique(scotch_results_Monocyte_7%>%filter(p_DTU_gene_adj<=0.05)%>%pull(genes))
genes_monocyte=intersect(a,b)
#B
a=unique(scotch_results_B_8%>%filter(p_DTU_gene_adj<=0.05)%>%pull(genes))
b=unique(scotch_results_B_7%>%filter(p_DTU_gene_adj<=0.05)%>%pull(genes))
genes_B=intersect(a,b)
#CD4
a=unique(scotch_results_CD4_8%>%filter(p_DTU_gene_adj<=0.05)%>%pull(genes))
b=unique(scotch_results_CD4_7%>%filter(p_DTU_gene_adj<=0.05)%>%pull(genes))
genes_CD4=intersect(a,b)
#CD8
a=unique(scotch_results_CD8_8%>%filter(p_DTU_gene_adj<=0.05)%>%pull(genes))
b=unique(scotch_results_CD8_7%>%filter(p_DTU_gene_adj<=0.05)%>%pull(genes))
genes_CD8=intersect(a,b)
#NK
a=unique(scotch_results_NK_8%>%filter(p_DTU_gene_adj<=0.05)%>%pull(genes))
b=unique(scotch_results_NK_7%>%filter(p_DTU_gene_adj<=0.05)%>%pull(genes))
genes_NK=intersect(a,b)
gene_list = list(NK=genes_NK,CD8=genes_CD8,CD4=genes_CD4,B=genes_B,Monocyte=genes_monocyte)

m = make_comb_mat(gene_list)
ss = set_size(m)
cs = comb_size(m)

##############b. upset plot #####################

ht = UpSet(m, 
           set_order = order(ss),
           comb_order = order(comb_degree(m), -cs),
           pt_size = unit(2, "mm"),
           top_annotation = HeatmapAnnotation(annotation_name_gp= gpar(fontsize = 10),
             "Intersections of\nDTU genes" = anno_barplot(cs, 
                                                  ylim = c(0, max(cs)*1.1),
                                                  border = FALSE, 
                                                  gp = gpar(fill = "black"), 
                                                  height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           right_annotation = rowAnnotation(annotation_name_gp= gpar(fontsize = 10),
             "DTU genes" = anno_barplot(ss,  baseline = 0,ylim = c(0, max(ss)*1.5),
                                               border = FALSE, 
                                               gp = gpar(fill = "black"), 
                                               width = unit(4, "cm")
             ),set_name = NULL), 
           left_annotation = NULL,row_names_gp = gpar(fontsize = 10),
           show_row_names = T);ht = draw(ht)
od = column_order(ht)
sod = row_order(ht)
decorate_annotation("Intersections of\nDTU genes", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})
decorate_annotation("DTU genes", {
  grid.text(rev(ss[sod]), y = seq_along(ss), x = unit(rev(ss[sod]), "native") + unit(10, "pt"), 
            default.units = "native", 
            gp = gpar(fontsize = 6, col = "#404040"))
})

###################################################################
#####----------a. DTU genes for DGE analysis  -----------------####
###################################################################
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/scotch")
scotch_results_B_8=read.csv('scotch_BvsAll_sample8.csv')
scotch_results_Monocyte_8=read.csv('scotch_MonocytesvsAll_sample8.csv')
scotch_results_CD4_8=read.csv('scotch_CD4vsAll_sample8.csv')
scotch_results_CD8_8=read.csv('scotch_CD8vsAll_sample8.csv')
scotch_results_NK_8=read.csv('scotch_NKvsAll_sample8.csv')

scotch_results_B_7=read.csv('scotch_BvsAll_sample7.csv')
scotch_results_Monocyte_7=read.csv('scotch_MonocytesvsAll_sample7.csv')
scotch_results_CD4_7=read.csv('scotch_CD4vsAll_sample7.csv')
scotch_results_CD8_7=read.csv('scotch_CD8vsAll_sample7.csv')
scotch_results_NK_7=read.csv('scotch_NKvsAll_sample7.csv')

#monocyte
scotch_results_Monocyte_8 = scotch_results_Monocyte_8%>%filter(p_gene_adj>0.05)%>%select(genes,p_DTU_gene_adj)%>%unique()%>%
  mutate(sig_DTU_gene=ifelse(p_DTU_gene_adj<=0.05, 1, 0))
scotch_results_Monocyte_7 = scotch_results_Monocyte_7%>%filter(p_gene_adj>0.05)%>%select(genes,p_DTU_gene_adj)%>%unique()%>%
  mutate(sig_DTU_gene=ifelse(p_DTU_gene_adj<=0.05, 1, 0))
monocyte = scotch_results_Monocyte_8%>%left_join(scotch_results_Monocyte_7, by='genes')%>%
  filter(is.na(p_DTU_gene_adj.x)==F & is.na(p_DTU_gene_adj.y)==F)%>%rowwise()%>%
  mutate(sig_DTU_gene = sig_DTU_gene.x+sig_DTU_gene.y)%>%
  mutate(sig_DTU_gene =factor(sig_DTU_gene))
monocyte$sig_DTU_gene = factor(monocyte$sig_DTU_gene,levels = c(0,1,2),labels = c('0','1','2'))

#B
scotch_results_B_8 = scotch_results_B_8%>%filter(p_gene_adj>0.05)%>%select(genes,p_DTU_gene_adj)%>%unique()%>%
  mutate(sig_DTU_gene=ifelse(p_DTU_gene_adj<=0.05, 1, 0))
scotch_results_B_7 = scotch_results_B_7%>%filter(p_gene_adj>0.05)%>%select(genes,p_DTU_gene_adj)%>%unique()%>%
  mutate(sig_DTU_gene=ifelse(p_DTU_gene_adj<=0.05, 1, 0))
B = scotch_results_B_8%>%left_join(scotch_results_B_7, by='genes')%>%
  filter(is.na(p_DTU_gene_adj.x)==F & is.na(p_DTU_gene_adj.y)==F)%>%rowwise()%>%
  mutate(sig_DTU_gene = sig_DTU_gene.x+sig_DTU_gene.y)%>%
  mutate(sig_DTU_gene =factor(sig_DTU_gene))
B$sig_DTU_gene = factor(B$sig_DTU_gene,levels = c(0,1,2),labels = c('0','1','2'))

#CD4
scotch_results_CD4_8 = scotch_results_CD4_8%>%filter(p_gene_adj>0.05)%>%select(genes,p_DTU_gene_adj)%>%unique()%>%
  mutate(sig_DTU_gene=ifelse(p_DTU_gene_adj<=0.05, 1, 0))
scotch_results_CD4_7 = scotch_results_CD4_7%>%filter(p_gene_adj>0.05)%>%select(genes,p_DTU_gene_adj)%>%unique()%>%
  mutate(sig_DTU_gene=ifelse(p_DTU_gene_adj<=0.05, 1, 0))
CD4 = scotch_results_CD4_8%>%left_join(scotch_results_CD4_7, by='genes')%>%
  filter(is.na(p_DTU_gene_adj.x)==F & is.na(p_DTU_gene_adj.y)==F)%>%rowwise()%>%
  mutate(sig_DTU_gene = sig_DTU_gene.x+sig_DTU_gene.y)%>%
  mutate(sig_DTU_gene =factor(sig_DTU_gene))
CD4$sig_DTU_gene = factor(CD4$sig_DTU_gene,levels = c(0,1,2),
                               labels = c('0','1','2'))
#CD8
scotch_results_CD8_8 = scotch_results_CD8_8%>%filter(p_gene_adj>0.05)%>%select(genes,p_DTU_gene_adj)%>%unique()%>%
  mutate(sig_DTU_gene=ifelse(p_DTU_gene_adj<=0.05, 1, 0))
scotch_results_CD8_7 = scotch_results_CD8_7%>%filter(p_gene_adj>0.05)%>%select(genes,p_DTU_gene_adj)%>%unique()%>%
  mutate(sig_DTU_gene=ifelse(p_DTU_gene_adj<=0.05, 1, 0))
CD8 = scotch_results_CD8_8%>%left_join(scotch_results_CD8_7, by='genes')%>%
  filter(is.na(p_DTU_gene_adj.x)==F & is.na(p_DTU_gene_adj.y)==F)%>%rowwise()%>%
  mutate(sig_DTU_gene = sig_DTU_gene.x+sig_DTU_gene.y)%>%
  mutate(sig_DTU_gene =factor(sig_DTU_gene))
CD8$sig_DTU_gene = factor(CD8$sig_DTU_gene,levels = c(0,1,2),
                               labels = c('0','1','2'))
#NK
scotch_results_NK_8 = scotch_results_NK_8%>%filter(p_gene_adj>0.05)%>%select(genes,p_DTU_gene_adj)%>%unique()%>%
  mutate(sig_DTU_gene=ifelse(p_DTU_gene_adj<=0.05, 1, 0))
scotch_results_NK_7 = scotch_results_NK_7%>%filter(p_gene_adj>0.05)%>%select(genes,p_DTU_gene_adj)%>%unique()%>%
  mutate(sig_DTU_gene=ifelse(p_DTU_gene_adj<=0.05, 1, 0))
NK = scotch_results_NK_8%>%left_join(scotch_results_NK_7, by='genes')%>%
  filter(is.na(p_DTU_gene_adj.x)==F & is.na(p_DTU_gene_adj.y)==F)%>%rowwise()%>%
  mutate(sig_DTU_gene = sig_DTU_gene.x+sig_DTU_gene.y)%>%
  mutate(sig_DTU_gene =factor(sig_DTU_gene))
NK$sig_DTU_gene = factor(NK$sig_DTU_gene,levels = c(0,1,2),
                               labels = c('0','1','2'))

df=rbind(monocyte%>%mutate(celltype='Monocyte'),
         B%>%mutate(celltype='B'),
         CD4%>%mutate(celltype='CD4'),
         CD8%>%mutate(celltype='CD8'),
         NK%>%mutate(celltype='NK'))%>%
  dplyr::select(sig_DTU_gene,celltype)%>%
  group_by(celltype,sig_DTU_gene)%>%tally()%>%
  mutate(sum_n=sum(n),pct=n/sum_n)

plot_a=ggplot(df,aes(x=celltype,y=pct,fill=sig_DTU_gene))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.4)+
  scale_y_continuous(sec.axis = sec_axis(~ 1 - .))+
  theme_bw()+
  scale_fill_manual(values=c('#ccece6','#66c2a4','#005824'),name = '# samples tested significant DTU')+
  labs(x='',y='gene percentage')+ggtitle('DTU analysis for non-DE genes')+
  theme(legend.position = 'bottom',plot.title = element_text(size=12));plot_a

######################################################################
#####------c. enrichment plot: genes with significant DTU---------####
######################################################################
#----pathway----#
temp=msigdbr::msigdbr_collections()
db <- msigdbr(species = "Homo sapiens", 
              category = "C2", subcategory = "CP:KEGG")
conv <- unique(db[, c("human_gene_symbol", "gs_exact_source","gs_name")])%>%
  filter(str_detect(gs_exact_source,'hsa05')==F)
pathway_list <- split(x = conv$human_gene_symbol, f = conv$gs_name)
#----gene rankes----#
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/scotch")
scotch_results_B_8=read.csv('scotch_BvsAll_sample8.csv')
scotch_results_Monocyte_8=read.csv('scotch_MonocytesvsAll_sample8.csv')
scotch_results_CD4_8=read.csv('scotch_CD4vsAll_sample8.csv')
scotch_results_CD8_8=read.csv('scotch_CD8vsAll_sample8.csv')
scotch_results_NK_8=read.csv('scotch_NKvsAll_sample8.csv')

scotch_results_B_7=read.csv('scotch_BvsAll_sample7.csv')
scotch_results_Monocyte_7=read.csv('scotch_MonocytesvsAll_sample7.csv')
scotch_results_CD4_7=read.csv('scotch_CD4vsAll_sample7.csv')
scotch_results_CD8_7=read.csv('scotch_CD8vsAll_sample7.csv')
scotch_results_NK_7=read.csv('scotch_NKvsAll_sample7.csv')

#monocyte
monocyte8=scotch_results_Monocyte_8%>%filter(is.na(p_DTU_gene_adj)==F)%>%
  dplyr::select(genes, p_DTU_gene_adj)%>%unique()%>%
  mutate(neg_log_10_padj=-log10(p_DTU_gene_adj))%>%
  arrange(-neg_log_10_padj)%>%
  mutate(neg_log_10_padj=ifelse(is.infinite(neg_log_10_padj),400,neg_log_10_padj))
monocyte7=scotch_results_Monocyte_7%>%filter(is.na(p_DTU_gene_adj)==F)%>%
  dplyr::select(genes, p_DTU_gene_adj)%>%unique()%>%
  mutate(neg_log_10_padj=-log10(p_DTU_gene_adj))%>%
  arrange(-neg_log_10_padj)%>%
  mutate(neg_log_10_padj=ifelse(is.infinite(neg_log_10_padj),400,neg_log_10_padj))
monocyte = left_join(monocyte7,monocyte8,by='genes')%>%rowwise()%>%
  mutate(statistics = (neg_log_10_padj.x +neg_log_10_padj.y)/2 )%>%
  filter(is.na(statistics)==F)%>%arrange(-statistics)%>%
  dplyr::select(genes,statistics)
monocyte_rank <- setNames(monocyte$statistics, monocyte$genes)

#b
B8=scotch_results_B_8%>%filter(is.na(p_DTU_gene_adj)==F)%>%
  dplyr::select(genes, p_DTU_gene_adj)%>%unique()%>%
  mutate(neg_log_10_padj=-log10(p_DTU_gene_adj))%>%
  arrange(-neg_log_10_padj)%>%
  mutate(neg_log_10_padj=ifelse(is.infinite(neg_log_10_padj),400,neg_log_10_padj))
B7=scotch_results_B_7%>%filter(is.na(p_DTU_gene_adj)==F)%>%
  dplyr::select(genes, p_DTU_gene_adj)%>%unique()%>%
  mutate(neg_log_10_padj=-log10(p_DTU_gene_adj))%>%
  arrange(-neg_log_10_padj)%>%
  mutate(neg_log_10_padj=ifelse(is.infinite(neg_log_10_padj),400,neg_log_10_padj))
B = left_join(B7,B8,by='genes')%>%rowwise()%>%
  mutate(statistics = (neg_log_10_padj.x +neg_log_10_padj.y)/2 )%>%
  filter(is.na(statistics)==F)%>%arrange(-statistics)%>%
  dplyr::select(genes,statistics)
B_rank <- setNames(B$statistics, B$genes)

#NK
NK8=scotch_results_NK_8%>%filter(is.na(p_DTU_gene_adj)==F)%>%
  dplyr::select(genes, p_DTU_gene_adj)%>%unique()%>%
  mutate(neg_log_10_padj=-log10(p_DTU_gene_adj))%>%
  arrange(-neg_log_10_padj)%>%
  mutate(neg_log_10_padj=ifelse(is.infinite(neg_log_10_padj),400,neg_log_10_padj))
NK7=scotch_results_NK_7%>%filter(is.na(p_DTU_gene_adj)==F)%>%
  dplyr::select(genes, p_DTU_gene_adj)%>%unique()%>%
  mutate(neg_log_10_padj=-log10(p_DTU_gene_adj))%>%
  arrange(-neg_log_10_padj)%>%
  mutate(neg_log_10_padj=ifelse(is.infinite(neg_log_10_padj),400,neg_log_10_padj))
NK = left_join(NK7,NK8,by='genes')%>%rowwise()%>%
  mutate(statistics = (neg_log_10_padj.x +neg_log_10_padj.y)/2 )%>%
  filter(is.na(statistics)==F)%>%arrange(-statistics)%>%
  dplyr::select(genes,statistics)
NK_rank <- setNames(NK$statistics, NK$genes)

#CD4
CD48=scotch_results_CD4_8%>%filter(is.na(p_DTU_gene_adj)==F)%>%
  dplyr::select(genes, p_DTU_gene_adj)%>%unique()%>%
  mutate(neg_log_10_padj=-log10(p_DTU_gene_adj))%>%
  arrange(-neg_log_10_padj)%>%
  mutate(neg_log_10_padj=ifelse(is.infinite(neg_log_10_padj),400,neg_log_10_padj))
CD47=scotch_results_CD4_7%>%filter(is.na(p_DTU_gene_adj)==F)%>%
  dplyr::select(genes, p_DTU_gene_adj)%>%unique()%>%
  mutate(neg_log_10_padj=-log10(p_DTU_gene_adj))%>%
  arrange(-neg_log_10_padj)%>%
  mutate(neg_log_10_padj=ifelse(is.infinite(neg_log_10_padj),400,neg_log_10_padj))
CD4 = left_join(CD47,CD48,by='genes')%>%rowwise()%>%
  mutate(statistics = (neg_log_10_padj.x +neg_log_10_padj.y)/2 )%>%
  filter(is.na(statistics)==F)%>%arrange(-statistics)%>%
  dplyr::select(genes,statistics)
CD4_rank <- setNames(CD4$statistics, CD4$genes)

#CD8
CD88=scotch_results_CD8_8%>%filter(is.na(p_DTU_gene_adj)==F)%>%
  dplyr::select(genes, p_DTU_gene_adj)%>%unique()%>%
  mutate(neg_log_10_padj=-log10(p_DTU_gene_adj))%>%
  arrange(-neg_log_10_padj)%>%
  mutate(neg_log_10_padj=ifelse(is.infinite(neg_log_10_padj),400,neg_log_10_padj))
CD87=scotch_results_CD8_7%>%filter(is.na(p_DTU_gene_adj)==F)%>%
  dplyr::select(genes, p_DTU_gene_adj)%>%unique()%>%
  mutate(neg_log_10_padj=-log10(p_DTU_gene_adj))%>%
  arrange(-neg_log_10_padj)%>%
  mutate(neg_log_10_padj=ifelse(is.infinite(neg_log_10_padj),400,neg_log_10_padj))
CD8 = left_join(CD87,CD88,by='genes')%>%rowwise()%>%
  mutate(statistics = (neg_log_10_padj.x +neg_log_10_padj.y)/2 )%>%
  filter(is.na(statistics)==F)%>%arrange(-statistics)%>%
  dplyr::select(genes,statistics)
CD8_rank <- setNames(CD8$statistics, CD8$genes)

#----enrichment----#
fgseaRes_monocyte <- fgsea(pathways = pathway_list, 
                  stats    = monocyte_rank,
                  minSize  = 15,
                  maxSize  = 500,
                  gseaParam = 0.5,scoreType='pos')
fgseaRes_B <- fgsea(pathways = pathway_list, 
                           stats    = B_rank,
                           minSize  = 15,
                           maxSize  = 500,
                           gseaParam = 0.5,scoreType='pos')
fgseaRes_NK <- fgsea(pathways = pathway_list, 
                    stats    = NK_rank,
                    minSize  = 15,
                    maxSize  = 500,
                    gseaParam = 0.5,scoreType='pos')
fgseaRes_CD4 <- fgsea(pathways = pathway_list, 
                     stats    = CD4_rank,
                     minSize  = 15,
                     maxSize  = 500,
                     gseaParam = 0.5,scoreType='pos')
fgseaRes_CD8 <- fgsea(pathways = pathway_list, 
                     stats    = CD8_rank,
                     minSize  = 15,
                     maxSize  = 500,
                     gseaParam = 0.5,scoreType='pos')

pathway_monocyte=fgseaRes_monocyte%>%filter(padj<=0.05)%>%arrange(-NES)
pathway_B=fgseaRes_B%>%filter(padj<=0.05)%>%arrange(-NES)
pathway_NK=fgseaRes_NK%>%filter(padj<=0.05)%>%arrange(-NES)
pathway_CD4=fgseaRes_CD4%>%filter(padj<=0.05)%>%arrange(-NES)
pathway_CD8=fgseaRes_CD8%>%filter(padj<=0.05)%>%arrange(-NES)

#pathway_monocyte=fora(pathways = pathway_list, genes = genes_monocyte,universe = gene_universe,minSize = 10, maxSize = 500)%>%
#  filter(padj<=0.05)%>%mutate(overlap_ratio=overlap/size)%>%arrange(padj)
#pathway_NK=fora(pathways = pathway_list, genes = genes_NK,universe = gene_universe,minSize = 10, maxSize = 500)%>%
#  filter(padj<=0.05)%>%mutate(overlap_ratio=overlap/size)%>%arrange(padj)
#pathway_B=fora(pathways = pathway_list, genes = genes_B,universe = gene_universe,minSize = 10, maxSize = 500)%>%
#  filter(padj<=0.05)%>%mutate(overlap_ratio=overlap/size)%>%arrange(padj)
#pathway_CD4=fora(pathways = pathway_list, genes = genes_CD4,universe = gene_universe,minSize = 10, maxSize = 500)%>%
#  filter(padj<=0.05)%>%mutate(overlap_ratio=overlap/size)%>%arrange(padj)
#pathway_CD8=fora(pathways = pathway_list, genes = genes_CD8,universe = gene_universe,minSize = 10, maxSize = 500)%>%
#  filter(padj<=0.05)%>%mutate(overlap_ratio=overlap/size)%>%arrange(padj)

top_n=15
pathway_results = rbind(pathway_monocyte%>%mutate(celltype='Monocyte')%>%head(top_n),
                        pathway_NK%>%mutate(celltype='NK')%>%head(top_n),
                        pathway_B%>%mutate(celltype='B')%>%head(top_n),
                        pathway_CD4%>%mutate(celltype='CD4')%>%head(top_n),
                        pathway_CD8%>%mutate(celltype='CD8')%>%head(top_n))
#split(f = pathway_results$celltype, x = pathway_results$pathway)

c = ggplot(pathway_results,aes(x=pathway,y=NES))+
  geom_bar(stat='identity',width = 0.5,aes(fill=-log10(padj)),color='black',alpha=0.8)+
  coord_flip()+facet_grid(~celltype)+theme_bw()+
  labs(x = '', y='NES')+
  scale_fill_gradientn(name = "-log10(p-adj)", 
                       colors = c( "white","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"),
                       values = scales::rescale(c(0,1.3,2,3,5,15)), 
                       limits = c(0, 15)) +
  theme(legend.position = 'bottom',axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size=6),
        axis.title.x = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size=6))+
  guides(fill = guide_colourbar(barwidth = 6, barheight = 0.4));c



#########################################################################
#####.  ------d. case study: nonsig gene sig DTUL IL27RA-------------####
#########################################################################
#monocyte8=scotch_results_Monocyte_8%>%filter(p_DTU_gene_adj<=0.05 & p_gene_adj>=0.05)%>%
#  pull(genes)
#monocyte7=scotch_results_Monocyte_7%>%filter(p_DTU_gene_adj<=0.05 & p_gene_adj>=0.05)%>%
#  pull(genes)
#overlap_genes = intersect(monocyte7,monocyte8)
#scotch_results_Monocyte_8%>%filter(genes%in%overlap_genes)%>%arrange(p_DTU_gene_adj)
#----IL27RA
gene="IL27RA"
df8 = scotch_results_Monocyte_8%>%filter(genes==gene)%>%
  dplyr::select(isoforms, TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
df8=rbind(df8%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'Monocyte')%>%
  dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
df8%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'Other')%>%
  dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))%>%
  mutate(Sample='Sample8')%>%
  mutate(IL27RA=str_remove(str_remove(isoforms,'IL27RA-'),'Group'))


df7 = scotch_results_Monocyte_7%>%filter(genes==gene)%>%
  dplyr::select(isoforms, TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
df7=rbind(df7%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'Monocyte')%>%
            dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
          df7%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'Other')%>%
            dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))%>%
  mutate(Sample='Sample7')%>%
  mutate(IL27RA=str_remove(str_remove(isoforms,'IL27RA-'),'Group'))
psudo = df7;psudo$IL27RA='other';psudo$`Transcript Usage`=0;psudo$max=0;psudo$min=0
df7=rbind(df7,psudo)
df = rbind(df7,df8)
df$IL27RA = str_replace(df$IL27RA,'novelIsoform','Novel-Isoform')
d=ggplot(df,aes(y=`Transcript Usage`,x=celltype,fill=IL27RA))+
  geom_bar(stat='identity',position ='dodge',width=0.4,color='black')+theme_bw()+
  geom_errorbar(aes(ymin=min,ymax=max),position =  position_dodge(width=0.4),width=0.2)+
  facet_grid(~Sample)+
  theme(legend.position = 'bottom',
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8))+
  labs(x='',y='Transcript Usage')+
  ylim(0,1)+
  scale_fill_brewer(palette = 'Set1');d


#########################################################################
#####.  ------e-g. case study: nonsig gene sig DTUL IL27RA-------------####
#########################################################################
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/coordinates")
#e gene-level
coordinates_df = read.csv('scotch_sample8_r10umap_coordinates_and_labels_with_barcodes.csv')
colnames(coordinates_df) = c('cell','UMAP1','UMAP2','celltype')
setwd("/Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/data/PBMC/scotch")
#gene files
sample8_CD4_gene=t(as.matrix(read.csv("nSample8gene_expression_TcellsCD4.csv",row.names='X')))[,gene]
sample8_CD8_gene=t(as.matrix(read.csv("nSample8gene_expression_TcellsCD8.csv",row.names='X')))[,gene]
sample8_B_gene=t(as.matrix(read.csv("nSample8gene_expression_Bcells.csv",row.names='X')))[,gene]
sample8_NK_gene=t(as.matrix(read.csv("nSample8gene_expression_NKcells.csv",row.names='X')))[,gene]
sample8_Monocytes_gene=t(as.matrix(read.csv("nSample8gene_expression_Monocytes.csv",row.names='X')))[,gene]
gene_expression = rbind(as.data.frame(sample8_CD4_gene)%>%dplyr::rename(gene_expression=sample8_CD4_gene),
                        as.data.frame(sample8_CD8_gene)%>%dplyr::rename(gene_expression=sample8_CD8_gene),
                        as.data.frame(sample8_B_gene)%>%dplyr::rename(gene_expression=sample8_B_gene),
                        as.data.frame(sample8_NK_gene)%>%dplyr::rename(gene_expression=sample8_NK_gene),
                        as.data.frame(sample8_Monocytes_gene)%>%dplyr::rename(gene_expression=sample8_Monocytes_gene))%>%
  mutate(cell = rownames(.))
#transcript files
sample8_CD4_transcript=read.csv("nSample8transcript_expression_TcellsCD4.csv",row.names = 'X')
sample8_CD4_transcript = sample8_CD4_transcript%>%as.matrix()%>%t()
sample8_CD8_transcript=read.csv("nSample8transcript_expression_TcellsCD8.csv",row.names = 'X')
sample8_CD8_transcript = sample8_CD8_transcript%>%as.matrix()%>%t()
sample8_B_transcript=read.csv("nSample8transcript_expression_Bcells.csv",row.names = 'X')
sample8_B_transcript = sample8_B_transcript%>%as.matrix()%>%t()
sample8_NK_transcript=read.csv("nSample8transcript_expression_NKcells.csv",row.names = 'X')
sample8_NK_transcript = sample8_NK_transcript%>%as.matrix()%>%t()
sample8_Monocytes_transcript=read.csv("nSample8transcript_expression_Monocytes.csv",row.names = 'X')
transcripts = data.frame(genes=str_remove(rownames(sample8_Monocytes_transcript),"-(ENST|novel|uncategorized).+"),
                                          transcripts=rownames(sample8_Monocytes_transcript))%>%
  filter(genes%in%gene)%>%pull(transcripts)
sample8_Monocytes_transcript = sample8_Monocytes_transcript%>%as.matrix()%>%t()

transcript_expression = rbind(sample8_CD4_transcript[,transcripts],
sample8_CD8_transcript[,transcripts],
sample8_B_transcript[,transcripts],
sample8_NK_transcript[,transcripts],
sample8_Monocytes_transcript[,transcripts])
transcript_expression = transcript_expression[rowSums(transcript_expression)>0,]
DTU = transcript_expression/rowSums(transcript_expression)
DTU_df = as.data.frame(DTU)%>%mutate(cell = rownames(.))

##merge
coordinates_df = coordinates_df%>%left_join(gene_expression)%>%left_join(DTU_df)
coordinates_df = coordinates_df[,-8]%>%pivot_longer(5:7)%>%
  mutate(name = case_when(name=='gene_expression'~'Gene Expression',
                          name=='IL27RA-novelIsoformGroup-3'~'NovelIsoform-3',
                          name=='IL27RA-ENST00000263379'~'ENST00000263379'))
coordinates_df_gene = coordinates_df%>%filter(name=='Gene Expression')
coordinates_df_transcript = coordinates_df%>%filter(name!='Gene Expression')%>%
  mutate(name = factor(name,levels = c('ENST00000263379','NovelIsoform-3')))


#e--- gene expression
e = ggplot(coordinates_df_gene,aes(x=UMAP1,y=UMAP2))+
  geom_point(aes(color=value),size=0.2,alpha=0.2)+
  stat_ellipse(data = subset(coordinates_df, celltype == "Monocytes"),
               aes(x = UMAP1, y = UMAP2),geom = "polygon",   
               fill = NA, color = "blue", linetype = "dashed",level = 0.99)+
  theme_light()+ggtitle('IL27RA Gene')+
  scale_color_continuous(low='blue',high='red',na.value='grey',name='Gene Expression')+
  theme(legend.key.width = unit(0.6, 'cm'),legend.key.height = unit(0.4, 'cm'),
        legend.position = 'bottom',
        plot.title = element_text(size=12))+
  guides(color = guide_colourbar(
    title.position = "left",title.vjust = 0.9,
    title.theme = element_text(margin = margin(r = 10)) 
  ));e
  
  

f = ggplot(coordinates_df_transcript%>%filter(name=='ENST00000263379'),aes(x=UMAP1,y=UMAP2))+
  geom_point(aes(color=value),size=0.2)+
  stat_ellipse(data = subset(coordinates_df_transcript, celltype == "Monocytes"),
               aes(x = UMAP1, y = UMAP2),geom = "polygon",   
               fill = NA, color = "blue", linetype = "dashed",level = 0.99)+
  theme_bw()+ggtitle('ENST00000263379')+
  scale_color_continuous(low='blue',high='red',na.value='grey',name='Transcript Usage')+
  theme(legend.key.width = unit(0.8, 'cm'),legend.key.height = unit(0.4, 'cm'),
        legend.position = 'bottom',
        plot.title = element_text(size=12))+
  guides(color = guide_colourbar(
    title.position = "left",title.vjust = 0.9,
    title.theme = element_text(margin = margin(r = 10)) 
  ));f

g = ggplot(coordinates_df_transcript%>%filter(name=='NovelIsoform-3'),aes(x=UMAP1,y=UMAP2))+
  geom_point(aes(color=value),size=0.2)+
  stat_ellipse(data = subset(coordinates_df_transcript, celltype == "Monocytes"),
               aes(x = UMAP1, y = UMAP2),geom = "polygon",   
               fill = NA, color = "blue", linetype = "dashed",level = 0.99)+
  theme_bw()+ggtitle('NovelIsoform-3')+
  scale_color_continuous(low='blue',high='red',na.value='grey',name='Transcript Usage')+
  theme(legend.key.width = unit(0.8, 'cm'),legend.key.height = unit(0.4, 'cm'),
        legend.position = 'bottom',
        plot.title = element_text(size=12))+
  guides(color = guide_colourbar(
    title.position = "left",title.vjust = 0.9,
    title.theme = element_text(margin = margin(r = 10)) 
  ));g

fg=ggarrange(f,g,ncol=2,common.legend = T,legend = 'bottom');fg
ef = ggarrange(e,fg,ncol = 2, widths = c(0.32,0.68),labels = c('e','f'));ef


###########################################################
#####.  ------sup.  novel isoform of  IL27RA-------------####
###########################################################

#conda activate trackplot
#cd /Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/sashimi

#trackplot \
#-e chr19:14031761-14053218 \
#-r genes_IL27RA.gtf \
#--density bamfile.tsv \
#-o IL27RA_sashimi.pdf \
#--dpi 300 \
#--width 12 \
#--height 1 --show-junction-num -t 10 \
#--intron-scale 1 



#######################################
#####.  ------Figure 3-------------####
#######################################
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/figure")
pdf('Figure3.pdf',width = 12,height =11)
#png('Figure3.png',width = 360*4,height =330*4)
grid.newpage()
#DTU DGE
print(plot_a,vp = viewport(y=0.84, x=0,height = 0.32, width=0.32, just = c("left")))
#subset plot
pushViewport(viewport(y=0.84,x=0.32,height = 0.32,width=0.68,just = c("left")))
ht = draw(ht,newpage=F)
od = column_order(ht)
sod = row_order(ht)
decorate_annotation("Intersections of\nDTU genes", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})
decorate_annotation("DTU genes", {
  grid.text(rev(ss[sod]), y = seq_along(ss), x = unit(rev(ss[sod]), "native") + unit(10, "pt"), 
            default.units = "native", 
            gp = gpar(fontsize = 6, col = "#404040"))
})
popViewport()
#pathway plot
print(c,vp = viewport(y=0.68, height = 0.35, width=1, just = c("top")))
#gene
print(d,vp = viewport(y=0.165,x=0,height = 0.33,width = 0.37,  just = c("left")))
print(ef,vp = viewport(y=0.165,x=0.37,height = 0.33,width = 0.61,  just = c("left")))
grid.text('a',x = 0.01,y=0.985,gp=gpar(fontsize=14, col="black",face='bold'))
grid.text('a',x = 0.01,y=0.985,gp=gpar(fontsize=14, col="black",face='bold'))
grid.text('a',x = 0.01,y=0.985,gp=gpar(fontsize=14, col="black",face='bold'))
grid.text('b',x = 0.32,y=0.985,gp=gpar(fontsize=14, col="black",face='bold'))
grid.text('b',x = 0.32,y=0.985,gp=gpar(fontsize=14, col="black",face='bold'))
grid.text('b',x = 0.32,y=0.985,gp=gpar(fontsize=14, col="black",face='bold'))
grid.text('c',x = 0.01,y=0.68,gp=gpar(fontsize=14, col="black",face='bold'))
grid.text('c',x = 0.01,y=0.68,gp=gpar(fontsize=14, col="black",face='bold'))
grid.text('c',x = 0.01,y=0.68,gp=gpar(fontsize=14, col="black",face='bold'))
grid.text('d',x = 0.01,y=0.31,gp=gpar(fontsize=14, col="black",face='bold'))
grid.text('d',x = 0.01,y=0.31,gp=gpar(fontsize=14, col="black",face='bold'))
grid.text('d',x = 0.01,y=0.31,gp=gpar(fontsize=14, col="black",face='bold'))
dev.off()
##############################################################
