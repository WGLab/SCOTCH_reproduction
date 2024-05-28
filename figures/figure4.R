
###############################AIF1###################################
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rtracklayer)
library(ggcoverage)
library(ggpattern)
library(ggtern)

gene='AIF1'
#gene='PTPRC'
###########function##########
process_mat = function(mat,isoforms){
  df=mat%>%as.data.frame()%>%filter(rowSums(.)>0)
  df = diag(1/rowSums(df))%*%as.matrix(df)
  df = df[,colnames(df)%in%isoforms]
  df = df%>%as.data.frame()%>%mutate(other = 1-rowSums(.))
  return(df)
}
##############################
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/scotch")
Sample8_CD4_transcript=read.csv("nSample8transcript_expression_TcellsCD4.csv",row.names = 'X')
Sample8_CD4_transcript = Sample8_CD4_transcript%>%as.matrix()%>%t()
Sample8_CD8_transcript=read.csv("nSample8transcript_expression_TcellsCD8.csv",row.names = 'X')
Sample8_CD8_transcript = Sample8_CD8_transcript%>%as.matrix()%>%t()
Sample8_B_transcript=read.csv("nSample8transcript_expression_Bcells.csv",row.names = 'X')
Sample8_B_transcript = Sample8_B_transcript%>%as.matrix()%>%t()
Sample8_NK_transcript=read.csv("nSample8transcript_expression_NKcells.csv",row.names = 'X')
Sample8_NK_transcript = Sample8_NK_transcript%>%as.matrix()%>%t()
Sample8_Monocytes_transcript=read.csv("nSample8transcript_expression_Monocytes.csv",row.names = 'X')
Sample8_Monocytes_transcript = Sample8_Monocytes_transcript%>%as.matrix()%>%t()
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/scotch")
results_B=read.csv('scotch_BvsAll_Sample8.csv')
results_Monocyte=read.csv('scotch_MonocytesvsAll_Sample8.csv')
results_CD4=read.csv('scotch_CD4vsAll_Sample8.csv')
results_CD8=read.csv('scotch_CD8vsAll_Sample8.csv')
results_NK=read.csv('scotch_NKvsAll_Sample8.csv')

isoforms = unique(c(results_Monocyte%>%filter(genes==gene)%>%pull(isoforms),
                    results_B%>%filter(genes==gene)%>%pull(isoforms),
                    results_CD4%>%filter(genes==gene)%>%pull(isoforms),
                    results_CD8%>%filter(genes==gene)%>%pull(isoforms),
                    results_NK%>%filter(genes==gene)%>%pull(isoforms)))

mat_CD4=Sample8_CD4_transcript[,str_starts(colnames(Sample8_CD4_transcript),gene)]
mat_CD8=Sample8_CD8_transcript[,str_starts(colnames(Sample8_CD8_transcript),gene)]
mat_B=Sample8_B_transcript[,str_starts(colnames(Sample8_B_transcript),gene)]
mat_NK=Sample8_NK_transcript[,str_starts(colnames(Sample8_NK_transcript),gene)]
mat_Monocyte=Sample8_Monocytes_transcript[,str_starts(colnames(Sample8_Monocytes_transcript),gene)]

##########################################################################################
library(forcats)
#sample 8 
results_B = results_B%>%filter(genes==gene)%>%dplyr::select(isoforms,TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
results_B = rbind(results_B%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'B')%>%
        dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
      results_B%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'otherB')%>%
        dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))

results_NK = results_NK%>%filter(genes==gene)%>%dplyr::select(isoforms,TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
results_NK = rbind(results_NK%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'NK')%>%
                    dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
                   results_NK%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'otherNK')%>%
                    dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))

results_CD4 = results_CD4%>%filter(genes==gene)%>%dplyr::select(isoforms,TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
results_CD4 = rbind(results_CD4%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'CD4')%>%
                     dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
                   results_CD4%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'otherCD4')%>%
                     dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))

results_CD8 = results_CD8%>%filter(genes==gene)%>%dplyr::select(isoforms,TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
results_CD8 = rbind(results_CD8%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'CD8')%>%
                      dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
                    results_CD8%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'otherCD8')%>%
                      dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))

results_Monocyte = results_Monocyte%>%filter(genes==gene)%>%dplyr::select(isoforms,TU1,TU2,TU_var1,TU_var2)%>%
  mutate(ymin1=TU1-TU_var1,ymin2 = TU2-TU_var2,ymax1=TU1+TU_var1,ymax2 = TU2+TU_var2)%>%
  dplyr::select(-contains("var"))
results_Monocyte = rbind(results_Monocyte%>%dplyr::select(isoforms,contains('1'))%>%mutate(celltype = 'Monocyte')%>%
                      dplyr::rename("Transcript Usage"="TU1","min"="ymin1","max"="ymax1"),
                    results_Monocyte%>%dplyr::select(isoforms,contains('2'))%>%mutate(celltype = 'otherMonocyte')%>%
                      dplyr::rename("Transcript Usage"="TU2","min"="ymin2","max"="ymax2"))

df_temp = expand.grid(isoforms = unique(c(results_B$isoforms,results_NK$isoforms,results_CD4$isoforms,
                               results_CD8$isoforms,results_Monocyte$isoforms)),TranscriptUsage = 0,min=0,max=0,
                     celltype=c('B','NK','CD4','CD8','Monocyte','otherB','otherNK','otherCD4','otherCD8','otherMonocyte'))
colnames(df_temp)[2] = c("Transcript Usage")


df_data=rbind(rbind(results_B,df_temp%>%filter(celltype%in%results_B$celltype & isoforms%in%results_B$isoforms==F))%>%
  mutate(`cell type` = 'B'),
rbind(results_NK,df_temp%>%filter(celltype%in%results_NK$celltype & isoforms%in%results_NK$isoforms==F))%>%
  mutate(`cell type` = 'NK'),
rbind(results_CD4,df_temp%>%filter(celltype%in%results_CD4$celltype & isoforms%in%results_CD4$isoforms==F))%>%
  mutate(`cell type` = 'CD4'),
rbind(results_CD8,df_temp%>%filter(celltype%in%results_CD8$celltype & isoforms%in%results_CD8$isoforms==F))%>%
  mutate(`cell type` = 'CD8'),
rbind(results_Monocyte,df_temp%>%filter(celltype%in%results_Monocyte$celltype & isoforms%in%results_Monocyte$isoforms==F))%>%
  mutate(`cell type` = 'Monocyte'))%>%
  mutate(isoforms=str_remove(isoforms,paste0(gene,'-')))%>%dplyr::rename(AIF1=isoforms)


df_data <- df_data %>%
  mutate(celltype = fct_recode(celltype,
                                  "other" = "otherB",
                                  "other" = "otherNK",
                                  "other" = "otherCD4",
                                  "other" = "otherCD8",
                                  "other" = "otherMonocyte"))
df_data$AIF1 = str_replace(df_data$AIF1,'novelIsoformGroup','Novel-Isoform')

setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/scotch")
results_B=read.csv('scotch_BvsAll_Sample8.csv')%>%filter(genes==gene)%>%mutate(celltype='B')
results_Monocyte=read.csv('scotch_MonocytesvsAll_Sample8.csv')%>%filter(genes==gene)%>%mutate(celltype='Monocyte')
results_CD4=read.csv('scotch_CD4vsAll_Sample8.csv')%>%filter(genes==gene)%>%mutate(celltype='CD4')
results_CD8=read.csv('scotch_CD8vsAll_Sample8.csv')%>%filter(genes==gene)%>%mutate(celltype='CD8')
results_NK=read.csv('scotch_NKvsAll_Sample8.csv')%>%filter(genes==gene)%>%mutate(celltype='NK')
results = rbind(results_B, results_Monocyte, results_CD4,results_CD8,results_NK)

a = ggplot(df_data,aes(x=celltype,y=`Transcript Usage`,fill=AIF1))+
  geom_bar(stat='identity',position='dodge',alpha=0.9,color='black',width = 0.4)+
  geom_errorbar(aes(ymin=min,ymax=max),position =  position_dodge(width=0.4),width=0.2)+
  labs(y='transcript usage',x='')+theme_bw()+
  scale_fill_brewer(palette = 'Set1')+
  theme(legend.position = 'bottom')+
  facet_grid(~`cell type`,scales = 'free');a


##########################################################################################
df_CD4=process_mat(mat_CD4,isoforms = isoforms)
df_CD8=process_mat(mat_CD8,isoforms = isoforms)
df_B=process_mat(mat_B,isoforms = isoforms)
df_NK=process_mat(mat_NK,isoforms = isoforms)
df_Monocyte=process_mat(mat_Monocyte,isoforms = isoforms)
df_data = rbind(df_CD4%>%mutate(celltype='CD4',cell=1:nrow(.)),
                df_CD8%>%mutate(celltype='CD8',cell=1:nrow(.)),
                df_B%>%mutate(celltype='B',cell=1:nrow(.)),
                df_NK%>%mutate(celltype='NK',cell=1:nrow(.)),
                df_Monocyte%>%mutate(celltype='Monocyte',cell=1:nrow(.)))%>%
  pivot_longer(cols=1:length(isoforms),names_to = 'Isoform',values_to = 'Transcript Usage')%>%
  arrange(celltype,Isoform,`Transcript Usage`)%>%
  mutate(Isoform = str_remove(Isoform,paste0(gene,'-')))

b = ggplot(df_data,aes(x=`Transcript Usage`))+
  geom_histogram(bins=50)+labs(y='cell count')+
  facet_grid(celltype~Isoform,scales = 'free')+
  labs(x='Transcript Usage')+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        text = element_text(size = 8))+scale_y_sqrt();b

##########################################################################################
df_data = rbind(df_Monocyte%>%mutate(celltype='Monocyte'),
                rbind(df_CD4,df_CD8,df_B,df_NK)%>%
                  mutate(celltype='Other Cells'))
set.seed(42)
df_data = rbind(df_data%>%filter(celltype=='Other Cells')%>%mutate(cell=1:nrow(.)),
                df_data%>%filter(celltype=='Monocyte')%>%sample_n(159)%>%mutate(cell=1:nrow(.)))

df_data = df_data%>%
  pivot_longer(cols=1:length(isoforms),names_to = 'Isoform',values_to = 'Transcript Usage')%>%
  arrange(celltype,Isoform,`Transcript Usage`)%>%
  mutate(Isoform=str_remove(Isoform,paste0(gene,'-')))%>%dplyr::rename(AIF1=Isoform)

df_data$AIF1 = str_replace(df_data$AIF1,'novelIsoformGroup','Novel-Isoform')
c = ggplot(df_data,aes(y=`Transcript Usage`,x=cell,fill=AIF1))+
  geom_bar(stat='identity',position = 'stack')+
  facet_wrap(celltype~.,scales = 'free',nrow=2)+theme_bw()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        legend.position = 'bottom')+labs(x='',y='Transcript Usage')+
  scale_fill_brewer(palette = 'Set1');c


ggplot(df_data,aes(y=`Transcript Usage`,x=cell,fill=AIF1))+
  geom_bar(stat='identity',position = 'stack')+
  facet_wrap(celltype~.,scales = 'free',ncol=2)+theme_bw()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        legend.position = 'bottom')+labs(x='',y='Transcript Usage')+
  scale_fill_brewer(palette = 'Set1')
##########################################################################################
#conda activate trackplot
#cd /Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/sashimi

#trackplot \
#-e chr6:31615217-31617021 \
#-r genes_AIF1.gtf \
#--density bamfile.tsv \
#-o AIF1_sashimi_new.pdf \
#--dpi 300 \
#--width 12 \
#--height 1.5 --show-junction-num -t 100 \
#--intron-scale 1 


####################################plot######################################################
setwd("~/Documents/ResearchProject/SingleCellLongReads/data/figure")
ab=ggarrange(a,c,ncol=2,common.legend = T,legend = 'bottom',widths = c(0.7,0.37),labels = c('a','b'));ab
P = ggarrange(ab,b,nrow = 2,heights = c(0.4,0.6),labels = c('','c'));P
ggsave('Figure4_abc.pdf',width = 12,height =9)

