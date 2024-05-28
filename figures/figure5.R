setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/coordinates")
#sample7
X_illumina_7 = read.csv("n7illumina1875019891umap_coordinates_and_labels_with_barcodes.csv")
X_R9_ED1_7 = read.csv("n7_R9_ED1umap_coordinates_and_labels_with_barcodes.csv")
X_R9_ED2_7 = read.csv("n7_R9_ED2umap_coordinates_and_labels_with_barcodes.csv")
X_R10_ED1_7 = read.csv("n7_R10_ED1umap_coordinates_and_labels_with_barcodes.csv")
X_R10_ED2_7 = read.csv("n7_R10_ED2umap_coordinates_and_labels_with_barcodes.csv")
X_R10_ED1_7_SCOTCH=read.csv("scotch_sample7_r10umap_coordinates_and_labels_with_barcodes.csv")
parse_illumina_7 = read.csv("Sample7ParseIlluminaumap_coordinates_and_labels_with_barcodes.csv")
parse_R10_ED1_7 = read.csv("nSample7ParseNanoporeumap_coordinates_and_labels_with_barcodes.csv")
parse_R10_ED1_7_SCOTCH=read.csv("scotch_sample7_parseumap_coordinates_and_labels_with_barcodes.csv" )

#sample8
X_illumina_8 = read.csv("n8illumina5554309213umap_coordinates_and_labels_with_barcodes.csv")
X_R9_ED1_8 = read.csv("n8_R9_ED1umap_coordinates_and_labels_with_barcodes.csv")
X_R9_ED2_8 = read.csv("n8_R9_ED2umap_coordinates_and_labels_with_barcodes.csv")
X_R10_ED1_8 = read.csv("n8_R10_ED1umap_coordinates_and_labels_with_barcodes.csv")
X_R10_ED2_8 = read.csv("n8_R10_ED2umap_coordinates_and_labels_with_barcodes.csv")
X_R10_ED1_8_SCOTCH=read.csv("scotch_sample8_r10umap_coordinates_and_labels_with_barcodes.csv")
parse_illumina_8 = read.csv("Sample8ParseIlluminaumap_coordinates_and_labels_with_barcodes.csv")
parse_R10_ED1_8 = read.csv("nSample8ParseNanoporeumap_coordinates_and_labels_with_barcodes.csv")
parse_R10_ED1_8_SCOTCH=read.csv("scotch_sample8_parseumap_coordinates_and_labels_with_barcodes.csv" )


df7 = rbind(X_illumina_7%>%mutate(label='10X-illumina'),
            X_R9_ED1_7%>%mutate(label='10X-R9ED1-vendor'),
            X_R9_ED2_7%>%mutate(label='10X-R9ED2-vendor'),
            X_R10_ED1_7%>%mutate(label='10X-R10ED1-vendor'),
            X_R10_ED2_7%>%mutate(label='10X-R10ED2-vendor'),
            X_R10_ED1_7_SCOTCH%>%mutate(label='10X-R10ED1-SCOTCH'),
            parse_illumina_7%>%mutate(label='Parse-illumina'),
            parse_R10_ED1_7%>%mutate(label='Parse-R10ED1-vendor'),
            parse_R10_ED1_7_SCOTCH%>%mutate(label='Parse-R10ED1-SCOTCH'))%>%
  mutate(label=factor(label,levels=unique(label)))


df8 = rbind(X_illumina_8%>%mutate(label='10X-illumina'),
            X_R9_ED1_8%>%mutate(label='10X-R9ED1-vendor'),
            X_R9_ED2_8%>%mutate(label='10X-R9ED2-vendor'),
            X_R10_ED1_8%>%mutate(label='10X-R10ED1-vendor'),
            X_R10_ED2_8%>%mutate(label='10X-R10ED2-vendor'),
            X_R10_ED1_8_SCOTCH%>%mutate(label='10X-R10ED1-SCOTCH'),
            parse_illumina_8%>%mutate(label='Parse-illumina'),
            parse_R10_ED1_8%>%mutate(label='Parse-R10ED1-vendor'),
            parse_R10_ED1_8_SCOTCH%>%mutate(label='Parse-R10ED1-SCOTCH'))%>%
  mutate(label=factor(label,levels=unique(label)))

a = ggplot(df7,aes(x=umap_1,y=umap_2,color=singleR.labels))+
  geom_point(size=0.2)+theme_bw()+
  scale_color_brewer(palette = 'Set1',name='')+
  facet_wrap(~label,nrow=3)+
  labs(x='UMAP1',y='UMAP2')+
  guides(color=guide_legend(override.aes=list(size=1.5)))+
  theme(legend.position = 'bottom',
        strip.text = element_text(size = 8));a

b = ggplot(df8,aes(x=umap_1,y=umap_2,color=singleR.labels))+
  geom_point(size=0.2)+theme_bw()+
  scale_color_brewer(palette = 'Set1',name='')+
  facet_wrap(~label,nrow=3)+
  labs(x='UMAP1',y='UMAP2')+
  guides(color=guide_legend(override.aes=list(size=1.5)))+
  theme(legend.position = 'bottom',
  strip.text = element_text(size = 8));b


df7_summary = df7%>%group_by(label,singleR.labels)%>%tally()%>%left_join(df7%>%group_by(label)%>%tally(),by='label')%>%
  ungroup()%>%rowwise()%>%mutate(percentage = n.x/n.y)%>%mutate(sample='Sample 7')
df8_summary = df8%>%group_by(label,singleR.labels)%>%tally()%>%left_join(df8%>%group_by(label)%>%tally(),by='label')%>%
  ungroup()%>%rowwise()%>%mutate(percentage = n.x/n.y)%>%mutate(sample='Sample 8')
clustering_summary = rbind(df7_summary,df8_summary)

c = ggplot(clustering_summary,aes(x=label,y=percentage,fill=singleR.labels))+
  geom_bar(stat='identity',width=0.4)+theme_bw()+
  scale_fill_brewer(palette = 'Set1',name='')+
  theme(axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        strip.text = element_text(size = 8),
        axis.title.x = element_text(size=8),
        legend.position = 'right')+
  facet_grid(~sample)+labs(x='',y='cell type composition');c


##########################-------number of DTU--------------###################
#########-------scotch
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/scotch")
scotch_results_B_8=read.csv('scotch_BvsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(scotch_results_B_8)
scotch_results_Monocyte_8=read.csv('scotch_MonocytesvsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(scotch_results_Monocyte_8)
scotch_results_CD4_8=read.csv('scotch_CD4vsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(scotch_results_CD4_8)
scotch_results_CD8_8=read.csv('scotch_CD8vsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(scotch_results_CD8_8)
scotch_results_NK_8=read.csv('scotch_NKvsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(scotch_results_NK_8)

scotch_results_B_7=read.csv('scotch_BvsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(scotch_results_B_7)
scotch_results_Monocyte_7=read.csv('scotch_MonocytesvsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(scotch_results_Monocyte_7)
scotch_results_CD4_7=read.csv('scotch_CD4vsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(scotch_results_CD4_7)
scotch_results_CD8_7=read.csv('scotch_CD8vsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(scotch_results_CD8_7)
scotch_results_NK_7=read.csv('scotch_NKvsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(scotch_results_NK_7)

scotch_B = intersect(scotch_results_B_8,scotch_results_B_7);length(scotch_B)
scotch_Monocyte = intersect(scotch_results_Monocyte_8,scotch_results_Monocyte_7);length(scotch_Monocyte)
scotch_NK = intersect(scotch_results_NK_8,scotch_results_NK_7);length(scotch_NK)
scotch_CD4 = intersect(scotch_results_CD4_8,scotch_results_CD4_7);length(scotch_CD4)
scotch_CD8 = intersect(scotch_results_CD8_8,scotch_results_CD8_7);length(scotch_CD8)


read.csv('scotch_BvsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%filter(genes=='CD74')
#######-----nanopore
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/nanopore")
nanopore_results_B_8=read.csv('nanopore_BvsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(nanopore_results_B_8)
nanopore_results_Monocyte_8=read.csv('nanopore_MonocytesvsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(nanopore_results_Monocyte_8)
nanopore_results_CD4_8=read.csv('nanopore_CD4vsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(nanopore_results_CD4_8)
nanopore_results_CD8_8=read.csv('nanopore_CD8vsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(nanopore_results_CD8_8)
nanopore_results_NK_8=read.csv('nanopore_NKvsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(nanopore_results_NK_8)

nanopore_results_B_7=read.csv('nanopore_BvsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(nanopore_results_B_7)
nanopore_results_Monocyte_7=read.csv('nanopore_MonocytesvsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(nanopore_results_Monocyte_7)
nanopore_results_CD4_7=read.csv('nanopore_CD4vsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(nanopore_results_CD4_7)
nanopore_results_CD8_7=read.csv('nanopore_CD8vsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(nanopore_results_CD8_7)
nanopore_results_NK_7=read.csv('nanopore_NKvsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(nanopore_results_NK_7)

nanopore_B = intersect(nanopore_results_B_8,nanopore_results_B_7);length(nanopore_B)
nanopore_Monocyte = intersect(nanopore_results_Monocyte_8,nanopore_results_Monocyte_7);length(nanopore_Monocyte)
nanopore_NK = intersect(nanopore_results_NK_8,nanopore_results_NK_7);length(nanopore_NK)
nanopore_CD4 = intersect(nanopore_results_CD4_8,nanopore_results_CD4_7);length(nanopore_CD4)
nanopore_CD8 = intersect(nanopore_results_CD8_8,nanopore_results_CD8_7);length(nanopore_CD8)


#######-----parse
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/parse")
parse_results_B_8=read.csv('parse_BvsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(parse_results_B_8)
parse_results_Monocyte_8=read.csv('parse_MonocytesvsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(parse_results_Monocyte_8)
parse_results_CD4_8=read.csv('parse_CD4vsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(parse_results_CD4_8)
parse_results_CD8_8=read.csv('parse_CD8vsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(parse_results_CD8_8)
parse_results_NK_8=read.csv('parse_NKvsAll_sample8.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(parse_results_NK_8)

parse_results_B_7=read.csv('parse_BvsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(parse_results_B_7)
parse_results_Monocyte_7=read.csv('parse_MonocytesvsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(parse_results_Monocyte_7)
parse_results_CD4_7=read.csv('parse_CD4vsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(parse_results_CD4_7)
parse_results_CD8_7=read.csv('parse_CD8vsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(parse_results_CD8_7)
parse_results_NK_7=read.csv('parse_NKvsAll_sample7.csv')%>%filter(p_DTU_gene_adj<=0.05)%>%
  pull(genes)%>%unique();length(parse_results_NK_7)

parse_B = intersect(parse_results_B_8,parse_results_B_7);length(parse_B)
parse_Monocyte = intersect(parse_results_Monocyte_8,parse_results_Monocyte_7);length(parse_Monocyte)
parse_NK = intersect(parse_results_NK_8,parse_results_NK_7);length(parse_NK)
parse_CD4 = intersect(parse_results_CD4_8,parse_results_CD4_7);length(parse_CD4)
parse_CD8 = intersect(parse_results_CD8_8,parse_results_CD8_7);length(parse_CD8)


df_scotch = data.frame(celltype = c('B','Monocyte','NK','CD4','CD8'),
           sample7=c(3190,4479,959,1033,475),sample8=c(3693,5477,1883,1019,260),
           both=c(2309,3766,636,468,121))%>%
  mutate(Sample7=both/sample7,Sample8=both/sample8)%>%
  dplyr::select(-both,-sample7,-sample8)%>%pivot_longer(2:3)
df_nanopore = data.frame(celltype = c('B','Monocyte','NK','CD4','CD8'),
                       sample7=c(2237,3228,684,549,126),
                       sample8=c(2610,3862,1383,474,93),
                       both=c(1266,2160,380,179,22))%>%
  mutate(Sample7=both/sample7,Sample8=both/sample8)%>%
  dplyr::select(-both,-sample7,-sample8)%>%pivot_longer(2:3)
df = rbind(df_scotch%>%mutate(pipeline='SCOTCH'),
           df_nanopore%>%mutate(pipeline='wf-single-cell'))
df$pipeline = factor(df$pipeline,levels =unique(df$pipeline))
df$pipeline = factor(df$pipeline,labels = c("10X-SCOTCH","10X-wf-single-cell"))


d = ggplot(df,aes(x=celltype,y=value,fill=pipeline))+
  geom_bar(stat='identity',position='dodge',width=0.4,color='black',alpha=0.8,show.legend=F)+
  facet_grid(~name)+theme_bw()+labs(x='',y='% reproducible DTU genes')+
  scale_fill_manual(values = c('#fc4e2a','#2171b5'))+ylim(0,1)+
  theme(#legend.position = 'n',
        strip.text = element_text(size = 8),
          axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8,angle = 45, hjust = 1));d


overlap_B = intersect(scotch_B,nanopore_B)
overlap_monocyte = intersect(scotch_Monocyte,nanopore_Monocyte)
overlap_NK = intersect(scotch_NK,nanopore_NK)
overlap_CD4 = intersect(scotch_CD4,nanopore_CD4)
overlap_CD8 = intersect(scotch_CD8,nanopore_CD8)

overlap_B_ = intersect(scotch_B,parse_results_B_7)
overlap_monocyte_ = intersect(scotch_Monocyte,parse_results_Monocyte_7)
overlap_NK_ = intersect(scotch_NK,parse_results_NK_7)
overlap_CD4_ = intersect(scotch_CD4,parse_results_CD4_7)
overlap_CD8_ = intersect(scotch_CD8,parse_results_CD8_7)

df=data.frame(celltype=c('B','Monocyte','NK','CD4','CD8'),
           `10X-SCOTCH` = c(length(overlap_B)/length(scotch_B),length(overlap_monocyte)/length(scotch_Monocyte),
                      length(overlap_NK)/length(scotch_NK),length(overlap_CD4)/length(scotch_CD4),
                      length(overlap_CD8)/length(scotch_CD8)),
           `10X-wf-single-cell` = c(length(overlap_B)/length(nanopore_B),length(overlap_monocyte)/length(nanopore_Monocyte),
                      length(overlap_NK)/length(nanopore_NK),length(overlap_CD4)/length(nanopore_CD4),
                      length(overlap_CD8)/length(nanopore_CD8)),
           `Parse-SCOTCH` = c(length(overlap_B_)/length(parse_results_B_7),length(overlap_monocyte_)/length(parse_results_Monocyte_7),
                              length(overlap_NK_)/length(parse_results_NK_7),length(overlap_CD4_)/length(parse_results_CD4_7),
                              length(overlap_CD8_)/length(parse_results_CD8_7)))%>%
  pivot_longer(2:4)

df$name = factor(df$name,levels = unique(df$name))
df$name = factor(df$name,labels = c("10X-SCOTCH","10X-wf-single-cell","Parse-SCOTCH"))

e = ggplot(df,aes(x=celltype,y=value,fill=name))+
  geom_bar(stat='identity',position='dodge',width=0.4,color='black',alpha=0.8)+
  theme_bw()+labs(x='',y='% consensus DTU genes')+
  scale_fill_manual(values = c('#fc4e2a','#2171b5','#8c6bb1'),name='')+ylim(0,1)+
  theme(legend.position = 'bottom',
        strip.text = element_text(size = 8),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8,angle = 45, hjust = 1));e

setwd("~/Documents/ResearchProject/SingleCellLongReads/data/figure")
ab=ggarrange(a,b,ncol=2,common.legend = T,legend = 'bottom',labels = c('a','b'))
de = ggarrange(d,e,ncol=2,common.legend = T,legend = 'bottom',labels = c('d','e'),widths = c(0.45,0.27));de
cde = ggarrange(c,de,ncol=2,labels = c('c',''),widths = c(0.55,0.45+0.27));cde
P = ggarrange(ab,cde,nrow=2,heights = c(0.63,0.37));P
ggsave('Figure5.png',width = 12,height =9)


######################################
####IFITM2 in scotch, nanopore, parse#cd37
parse_results_Monocyte_7[parse_results_Monocyte_7%in%overlap_monocyte]
gene='IFITM2'
###########function##########
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/scotch")
scotch_Monocyte_sample8=read.csv('scotch_MonocytesvsAll_Sample8.csv')
scotch_Monocyte_sample7=read.csv('scotch_MonocytesvsAll_Sample7.csv')
isoforms = unique(c(scotch_Monocyte_sample8%>%filter(genes==gene)%>%pull(isoforms),
                    scotch_Monocyte_sample7%>%filter(genes==gene)%>%pull(isoforms)));isoforms
results_Monocyte8 = scotch_Monocyte_sample8%>%filter(genes==gene)%>%dplyr::select(isoforms,TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
results_Monocyte8 = rbind(results_Monocyte8%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'Monocyte')%>%
                           dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
                          results_Monocyte8%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'other')%>%
                           dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))
results_Monocyte7 = scotch_Monocyte_sample7%>%filter(genes==gene)%>%dplyr::select(isoforms,TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
results_Monocyte7 = rbind(results_Monocyte7%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'Monocyte')%>%
                            dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
                          results_Monocyte7%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'other')%>%
                            dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))
scotch_results_Monocyte = rbind(results_Monocyte7%>%mutate(sample='Sample7'),results_Monocyte8%>%mutate(sample='Sample8'))
scotch_results_Monocyte$isoforms = str_remove(scotch_results_Monocyte$isoforms,paste0(gene,'-'));scotch_results_Monocyte$isoforms

setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/nanopore")
nanopore_Monocyte_sample8=read.csv('nanopore_MonocytesvsAll_Sample8.csv')
nanopore_Monocyte_sample7=read.csv('nanopore_MonocytesvsAll_Sample7.csv')
isoforms = unique(c(nanopore_Monocyte_sample8%>%filter(genes==gene)%>%pull(isoforms),
                    nanopore_Monocyte_sample7%>%filter(genes==gene)%>%pull(isoforms)))
results_Monocyte8 = nanopore_Monocyte_sample8%>%filter(genes==gene)%>%dplyr::select(isoforms,TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
results_Monocyte8 = rbind(results_Monocyte8%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'Monocyte')%>%
                            dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
                          results_Monocyte8%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'other')%>%
                            dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))
results_Monocyte7 = nanopore_Monocyte_sample7%>%filter(genes==gene)%>%dplyr::select(isoforms,TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
results_Monocyte7 = rbind(results_Monocyte7%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'Monocyte')%>%
                            dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
                          results_Monocyte7%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'other')%>%
                            dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))
nanopore_results_Monocyte = rbind(results_Monocyte7%>%mutate(sample='Sample7'),results_Monocyte8%>%mutate(sample='Sample8'))

setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/parse")
parse_Monocyte_sample8=read.csv('parse_MonocytesvsAll_Sample8.csv')
parse_Monocyte_sample7=read.csv('parse_MonocytesvsAll_Sample7.csv')
isoforms = unique(c(parse_Monocyte_sample8%>%filter(genes==gene)%>%pull(isoforms),
                    parse_Monocyte_sample7%>%filter(genes==gene)%>%pull(isoforms)))
results_Monocyte8 = parse_Monocyte_sample8%>%filter(genes==gene)%>%dplyr::select(isoforms,TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
results_Monocyte8 = rbind(results_Monocyte8%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'Monocyte')%>%
                            dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
                          results_Monocyte8%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'other')%>%
                            dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))
results_Monocyte7 = parse_Monocyte_sample7%>%filter(genes==gene)%>%dplyr::select(isoforms,TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
results_Monocyte7 = rbind(results_Monocyte7%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'Monocyte')%>%
                            dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
                          results_Monocyte7%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'other')%>%
                            dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))
parse_results_Monocyte = rbind(results_Monocyte7%>%mutate(sample='Sample7'),results_Monocyte8%>%mutate(sample='Sample8'))
parse_results_Monocyte$isoforms = str_remove(parse_results_Monocyte$isoforms,paste0(gene,'-'))

results_Monocyte = rbind(scotch_results_Monocyte%>%mutate(Pipeline='10X-SCOTCH'),
                         nanopore_results_Monocyte%>%mutate(Pipeline='10X-wf-single-cell'),
                         parse_results_Monocyte%>%mutate(Pipeline='Parse-SCOTCH'))

ggplot(results_Monocyte,aes(x=celltype,y=`Transcript Usage`,fill=isoforms))+
  geom_bar(stat='identity',position='dodge',alpha=0.9,color='black',width = 0.4)+
  geom_errorbar(aes(ymin=min,ymax=max),position =  position_dodge(width=0.4),width=0.2)+
  labs(y='transcript usage',x='')+theme_bw()+
  scale_fill_brewer(palette = 'Set1')+
  theme(legend.position = 'bottom')+
  facet_grid(sample~Pipeline,scales = 'free')



#conda activate trackplot
#cd /Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/sashimi

#trackplot \
#-e chr11:307631-315272 \
#-r genes_IFITM2.gtf \
#--density bamfile_IFITM2.tsv \
#-o IFITM2_sashimi.pdf \
#--dpi 300 \
#--width 12 \
#--height 1.5 --show-junction-num -t 100 \
#--intron-scale 1
