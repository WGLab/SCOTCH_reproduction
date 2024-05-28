################################################################################################
#########################scotch statistics analysis in scotch count matrix#################### 
################################################################################################
library(ggpubr)
library(stringr)
library(plyr)
library(dplyr)
source('/Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/Project_single-cell-long-read-transcriptomic-analysis/src/cluster.R')

####################################################################
######################--------sample8----------#####################
####################################################################
setwd("/Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/data/PBMC/parse")

#gene files
sample8_CD4_gene=t(as.matrix(read.csv("nSample8gene_expression_TcellsCD4.csv",row.names='X')))
sample8_CD8_gene=t(as.matrix(read.csv("nSample8gene_expression_TcellsCD8.csv",row.names='X')))
sample8_B_gene=t(as.matrix(read.csv("nSample8gene_expression_Bcells.csv",row.names='X')))
sample8_NK_gene=t(as.matrix(read.csv("nSample8gene_expression_NKcells.csv",row.names='X')))
sample8_Monocytes_gene=t(as.matrix(read.csv("nSample8gene_expression_Monocytes.csv",row.names='X')))

#transcript files
sample8_CD4_transcript=read.csv("nSample8transcript_expression_TcellsCD4.csv",row.names = 'X')
gene_transcript_CD4_df = data.frame(genes=str_remove(rownames(sample8_CD4_transcript),"-(ENST|novel|uncategorized).+"),
                                    transcripts=rownames(sample8_CD4_transcript))
sample8_CD4_transcript = sample8_CD4_transcript%>%as.matrix()%>%t()

sample8_CD8_transcript=read.csv("nSample8transcript_expression_TcellsCD8.csv",row.names = 'X')
gene_transcript_CD8_df = data.frame(genes=str_remove(rownames(sample8_CD8_transcript),"-(ENST|novel|uncategorized).+"),
                                    transcripts=rownames(sample8_CD8_transcript))
sample8_CD8_transcript = sample8_CD8_transcript%>%as.matrix()%>%t()

sample8_B_transcript=read.csv("nSample8transcript_expression_Bcells.csv",row.names = 'X')
gene_transcript_B_df = data.frame(genes=str_remove(rownames(sample8_B_transcript),"-(ENST|novel|uncategorized).+"),
                                  transcripts=rownames(sample8_B_transcript))
sample8_B_transcript = sample8_B_transcript%>%as.matrix()%>%t()

sample8_NK_transcript=read.csv("nSample8transcript_expression_NKcells.csv",row.names = 'X')
gene_transcript_NK_df = data.frame(genes=str_remove(rownames(sample8_NK_transcript),"-(ENST|novel|uncategorized).+"),
                                   transcripts=rownames(sample8_NK_transcript))
sample8_NK_transcript = sample8_NK_transcript%>%as.matrix()%>%t()

sample8_Monocytes_transcript=read.csv("nSample8transcript_expression_Monocytes.csv",row.names = 'X')
gene_transcript_Monocytes_df = data.frame(genes=str_remove(rownames(sample8_Monocytes_transcript),"-(ENST|novel|uncategorized).+"),
                                          transcripts=rownames(sample8_Monocytes_transcript))
sample8_Monocytes_transcript = sample8_Monocytes_transcript%>%as.matrix()%>%t()

setwd("/Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/parse")
###################CD4 VS all##################
#human PBMC samples----gene level: CD4 VS all#
df_gene = scotch_gene(sample8_CD4_gene, rbind(sample8_B_gene,sample8_CD8_gene,sample8_NK_gene,sample8_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD4 VS CD8#
df_transcript = scotch_transcript(gene_transcript_CD4_df,gene_transcript_CD4_df, 
                                  sample8_CD4_transcript, 
                                  rbind(sample8_B_transcript,sample8_CD8_transcript,sample8_NK_transcript,sample8_Monocytes_transcript),
                                  ncores=10)
df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
write.csv(df_scotch,'parse_CD4vsAll_sample8.csv')

###################CD8 VS all##################
#human PBMC samples----gene level: CD8 VS all#
df_gene = scotch_gene(sample8_CD8_gene, rbind(sample8_B_gene,sample8_CD4_gene,sample8_NK_gene,sample8_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD8 VS all#
df_transcript = scotch_transcript(gene_transcript_CD8_df,gene_transcript_B_df, sample8_CD8_transcript, 
                                  rbind(sample8_B_transcript,sample8_CD4_transcript,sample8_NK_transcript,sample8_Monocytes_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
write.csv(df_scotch,'parse_CD8vsAll_sample8.csv')


###################B VS all##################
#human PBMC samples----gene level: CD8 VS all#
df_gene = scotch_gene(sample8_B_gene, 
                      rbind(sample8_CD8_gene,sample8_CD4_gene,sample8_NK_gene,sample8_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD8 VS all#
df_transcript = scotch_transcript(gene_transcript_B_df,gene_transcript_CD8_df, 
                                  sample8_B_transcript, 
                                  rbind(sample8_CD8_transcript,sample8_CD4_transcript,sample8_NK_transcript,sample8_Monocytes_transcript),
                                  ncores=10)
df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
write.csv(df_scotch,'parse_BvsAll_sample8.csv')

###################NK VS all##################
#human PBMC samples----gene level: CD8 VS all#
df_gene = scotch_gene(sample8_NK_gene,
                      rbind(sample8_CD8_gene,sample8_CD4_gene,sample8_B_gene,sample8_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD8 VS all#
df_transcript = scotch_transcript(gene_transcript_NK_df,gene_transcript_CD8_df, 
                                  sample8_NK_transcript, 
                                  rbind(sample8_CD8_transcript,sample8_CD4_transcript,sample8_B_transcript,sample8_Monocytes_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
write.csv(df_scotch,'parse_NKvsAll_sample8.csv')

###################Monocytes VS all##################
#human PBMC samples----gene level: CD8 VS all#
df_gene = scotch_gene(sample8_Monocytes_gene,
                      rbind(sample8_CD8_gene,sample8_CD4_gene,sample8_B_gene,sample8_NK_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD8 VS all#
df_transcript = scotch_transcript(gene_transcript_Monocytes_df,gene_transcript_CD8_df, 
                                  sample8_Monocytes_transcript, 
                                  rbind(sample8_CD8_transcript,sample8_CD4_transcript,sample8_B_transcript,sample8_NK_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
write.csv(df_scotch,'parse_MonocytesvsAll_sample8.csv')



####################################################################
######################--------sample7----------#####################
####################################################################
setwd("/Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/data/PBMC/parse")

#gene files
sample7_CD4_gene=t(as.matrix(read.csv("nSample7gene_expression_TcellsCD4.csv",row.names='X')))
sample7_CD8_gene=t(as.matrix(read.csv("nSample7gene_expression_TcellsCD8.csv",row.names='X')))
sample7_B_gene=t(as.matrix(read.csv("nSample7gene_expression_Bcells.csv",row.names='X')))
sample7_NK_gene=t(as.matrix(read.csv("nSample7gene_expression_NKcells.csv",row.names='X')))
sample7_Monocytes_gene=t(as.matrix(read.csv("nSample7gene_expression_Monocytes.csv",row.names='X')))

#transcript files
sample7_CD4_transcript=read.csv("nSample7transcript_expression_TcellsCD4.csv",row.names = 'X')
gene_transcript_CD4_df = data.frame(genes=str_remove(rownames(sample7_CD4_transcript),"-(ENST|novel|uncategorized).+"),
                                    transcripts=rownames(sample7_CD4_transcript))
sample7_CD4_transcript = sample7_CD4_transcript%>%as.matrix()%>%t()

sample7_CD8_transcript=read.csv("nSample7transcript_expression_TcellsCD8.csv",row.names = 'X')
gene_transcript_CD8_df = data.frame(genes=str_remove(rownames(sample7_CD8_transcript),"-(ENST|novel|uncategorized).+"),
                                    transcripts=rownames(sample7_CD8_transcript))
sample7_CD8_transcript = sample7_CD8_transcript%>%as.matrix()%>%t()

sample7_B_transcript=read.csv("nSample7transcript_expression_Bcells.csv",row.names = 'X')
gene_transcript_B_df = data.frame(genes=str_remove(rownames(sample7_B_transcript),"-(ENST|novel|uncategorized).+"),
                                  transcripts=rownames(sample7_B_transcript))
sample7_B_transcript = sample7_B_transcript%>%as.matrix()%>%t()


sample7_NK_transcript=read.csv("nsample7transcript_expression_NKcells.csv",row.names = 'X')
gene_transcript_NK_df = data.frame(genes=str_remove(rownames(sample7_NK_transcript),"-(ENST|novel|uncategorized).+"),
                                   transcripts=rownames(sample7_NK_transcript))
sample7_NK_transcript = sample7_NK_transcript%>%as.matrix()%>%t()

sample7_Monocytes_transcript=read.csv("nsample7transcript_expression_Monocytes.csv",row.names = 'X')
gene_transcript_Monocytes_df = data.frame(genes=str_remove(rownames(sample7_Monocytes_transcript),"-(ENST|novel|uncategorized).+"),
                                          transcripts=rownames(sample7_Monocytes_transcript))
sample7_Monocytes_transcript = sample7_Monocytes_transcript%>%as.matrix()%>%t()

setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/parse")
###################CD4 VS all##################
#human PBMC samples----gene level: CD4 VS all#
df_gene = scotch_gene(sample7_CD4_gene, rbind(sample7_B_gene,sample7_CD8_gene,sample7_NK_gene,sample7_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD4 VS CD8#
df_transcript = scotch_transcript(gene_transcript_CD4_df,gene_transcript_B_df, 
                                  sample7_CD4_transcript, rbind(sample7_B_transcript,sample7_CD8_transcript,sample7_NK_transcript,sample7_Monocytes_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
write.csv(df_scotch,'parse_CD4vsAll_sample7.csv')

###################CD8 VS all##################
#human PBMC samples----gene level: CD8 VS all#
df_gene = scotch_gene(sample7_CD8_gene, rbind(sample7_B_gene,sample7_CD4_gene,sample7_NK_gene,sample7_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD8 VS all#
df_transcript = scotch_transcript(gene_transcript_CD8_df,gene_transcript_B_df, sample7_CD8_transcript, 
                                  rbind(sample7_B_transcript,sample7_CD4_transcript,sample7_NK_transcript,sample7_Monocytes_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
write.csv(df_scotch,'parse_CD8vsAll_sample7.csv')


###################B VS all##################
#human PBMC samples----gene level: CD8 VS all#
df_gene = scotch_gene(sample7_B_gene, 
                      rbind(sample7_CD8_gene,sample7_CD4_gene,sample7_NK_gene,sample7_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD8 VS all#
df_transcript = scotch_transcript(gene_transcript_B_df,gene_transcript_CD8_df, 
                                  sample7_B_transcript, 
                                  rbind(sample7_CD8_transcript,sample7_CD4_transcript,sample7_NK_transcript,sample7_Monocytes_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
write.csv(df_scotch,'parse_BvsAll_sample7.csv')

###################NK VS all##################
#human PBMC samples----gene level: CD8 VS all#
df_gene = scotch_gene(sample7_NK_gene,
                      rbind(sample7_CD8_gene,sample7_CD4_gene,sample7_B_gene,sample7_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD8 VS all#
df_transcript = scotch_transcript(gene_transcript_NK_df,gene_transcript_CD8_df, 
                                  sample7_NK_transcript, 
                                  rbind(sample7_CD8_transcript,sample7_CD4_transcript,sample7_B_transcript,sample7_Monocytes_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
write.csv(df_scotch,'parse_NKvsAll_sample7.csv')

###################Monocytes VS all##################
#human PBMC samples----gene level: CD8 VS all#
df_gene = scotch_gene(sample7_Monocytes_gene,
                      rbind(sample7_CD8_gene,sample7_CD4_gene,sample7_B_gene,sample7_NK_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD8 VS all#
df_transcript = scotch_transcript(gene_transcript_Monocytes_df,gene_transcript_CD8_df, 
                                  sample7_Monocytes_transcript, 
                                  rbind(sample7_CD8_transcript,sample7_CD4_transcript,sample7_B_transcript,sample7_NK_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
write.csv(df_scotch,'parse_MonocytesvsAll_sample7.csv')


