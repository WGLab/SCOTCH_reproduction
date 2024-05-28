################################################################################################
#########################scotch statistics analysis in nanopore count matrix#################### 
################################################################################################
library(ggpubr)
library(stringr)
library(plyr)
library(dplyr)
source('/Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/Project_single-cell-long-read-transcriptomic-analysis/src/cluster.R')

####################################################################
######################--------sample8----------#####################
####################################################################
setwd("/Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/data/PBMC/nanopore")

#gene files
sample8_CD4_gene=t(as.matrix(read.csv("Sample8gene_expression_TcellsCD4.csv",row.names='X')))
sample8_CD8_gene=t(as.matrix(read.csv("Sample8gene_expression_TcellsCD8.csv",row.names='X')))
sample8_B_gene=t(as.matrix(read.csv("Sample8gene_expression_Bcells.csv",row.names='X')))
sample8_NK_gene=t(as.matrix(read.csv("Sample8gene_expression_NKcells.csv",row.names='X')))
sample8_Monocytes_gene=t(as.matrix(read.csv("Sample8gene_expression_Monocytes.csv",row.names='X')))

#transcript files
sample8_CD4_transcript=read.csv("Sample8transcript_expression_TcellsCD4_genename.csv",row.names = 'X')
gene_transcript_CD4_df = data.frame(genes=sample8_CD4_transcript$Gene,transcripts=sample8_CD4_transcript$ID)
sample8_CD4_transcript = sample8_CD4_transcript%>%dplyr::select(-ID,-Gene)%>%as.matrix()%>%t()

sample8_CD8_transcript=read.csv("Sample8transcript_expression_TcellsCD8_genename.csv",row.names = 'X')
gene_transcript_CD8_df = data.frame(genes=sample8_CD8_transcript$Gene,transcripts=sample8_CD8_transcript$ID)
sample8_CD8_transcript = sample8_CD8_transcript%>%dplyr::select(-ID,-Gene)%>%as.matrix()%>%t()

sample8_B_transcript=read.csv("Sample8transcript_expression_Bcells_genename.csv",row.names = 'X')
gene_transcript_B_df = data.frame(genes=sample8_B_transcript$Gene,transcripts=sample8_B_transcript$ID)
sample8_B_transcript = sample8_B_transcript%>%dplyr::select(-ID,-Gene)%>%as.matrix()%>%t()

sample8_NK_transcript=read.csv("Sample8transcript_expression_NKcells_genename.csv",row.names = 'X')
gene_transcript_NK_df = data.frame(genes=sample8_NK_transcript$Gene,transcripts=sample8_NK_transcript$ID)
sample8_NK_transcript = sample8_NK_transcript%>%dplyr::select(-ID,-Gene)%>%as.matrix()%>%t()

sample8_Monocytes_transcript=read.csv("Sample8transcript_expression_Monocytes_genename.csv",row.names = 'X')
gene_transcript_Monocytes_df = data.frame(genes=sample8_Monocytes_transcript$Gene,transcripts=sample8_Monocytes_transcript$ID)
sample8_Monocytes_transcript = sample8_Monocytes_transcript%>%dplyr::select(-ID,-Gene)%>%as.matrix()%>%t()

setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/nanopore")
###################CD4 VS all##################
#human PBMC samples----gene level: CD4 VS all#
df_gene = scotch_gene(sample8_CD4_gene, rbind(sample8_B_gene,sample8_CD8_gene,sample8_NK_gene,sample8_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD4 VS CD8#
df_transcript = scotch_transcript(gene_transcript_CD4_df,gene_transcript_B_df, 
                                  sample8_CD4_transcript, rbind(sample8_B_transcript,sample8_CD8_transcript,sample8_NK_transcript,sample8_Monocytes_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
#write.csv(df_scotch,'nanopore_CD4vsAll_sample8.csv')

###################CD8 VS all##################
#human PBMC samples----gene level: CD8 VS all#
df_gene = scotch_gene(sample8_CD8_gene, rbind(sample8_B_gene,sample8_CD4_gene,sample8_NK_gene,sample8_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD8 VS all#
df_transcript = scotch_transcript(gene_transcript_CD8_df,gene_transcript_B_df, sample8_CD8_transcript, 
                                  rbind(sample8_B_transcript,sample8_CD4_transcript,sample8_NK_transcript,sample8_Monocytes_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
#write.csv(df_scotch,'nanopore_CD8vsAll_sample8.csv')


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
#write.csv(df_scotch,'nanopore_BvsAll_sample8.csv')

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
#write.csv(df_scotch,'nanopore_NKvsAll_sample8.csv')

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
#write.csv(df_scotch,'nanopore_MonocytesvsAll_sample8.csv')



####################################################################
######################--------sample7----------#####################
####################################################################
setwd("/Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/data/PBMC/nanopore")

#gene files
sample7_CD4_gene=t(as.matrix(read.csv("sample7gene_expression_TcellsCD4.csv",row.names='X')))
sample7_CD8_gene=t(as.matrix(read.csv("sample7gene_expression_TcellsCD8.csv",row.names='X')))
sample7_B_gene=t(as.matrix(read.csv("sample7gene_expression_Bcells.csv",row.names='X')))
sample7_NK_gene=t(as.matrix(read.csv("sample7gene_expression_NKcells.csv",row.names='X')))
sample7_Monocytes_gene=t(as.matrix(read.csv("sample7gene_expression_Monocytes.csv",row.names='X')))

#transcript files
sample7_CD4_transcript=read.csv("sample7transcript_expression_TcellsCD4_genename.csv",row.names = 'X')
gene_transcript_CD4_df = data.frame(genes=sample7_CD4_transcript$Gene,transcripts=sample7_CD4_transcript$ID)
sample7_CD4_transcript = sample7_CD4_transcript%>%dplyr::select(-ID,-Gene)%>%as.matrix()%>%t()

sample7_CD8_transcript=read.csv("sample7transcript_expression_TcellsCD8_genename.csv",row.names = 'X')
gene_transcript_CD8_df = data.frame(genes=sample7_CD8_transcript$Gene,transcripts=sample7_CD8_transcript$ID)
sample7_CD8_transcript = sample7_CD8_transcript%>%dplyr::select(-ID,-Gene)%>%as.matrix()%>%t()

sample7_B_transcript=read.csv("sample7transcript_expression_Bcells_genename.csv",row.names = 'X')
gene_transcript_B_df = data.frame(genes=sample7_B_transcript$Gene,transcripts=sample7_B_transcript$ID)
sample7_B_transcript = sample7_B_transcript%>%dplyr::select(-ID,-Gene)%>%as.matrix()%>%t()

sample7_NK_transcript=read.csv("sample7transcript_expression_NKcells_genename.csv",row.names = 'X')
gene_transcript_NK_df = data.frame(genes=sample7_NK_transcript$Gene,transcripts=sample7_NK_transcript$ID)
sample7_NK_transcript = sample7_NK_transcript%>%dplyr::select(-ID,-Gene)%>%as.matrix()%>%t()

sample7_Monocytes_transcript=read.csv("sample7transcript_expression_Monocytes_genename.csv",row.names = 'X')
gene_transcript_Monocytes_df = data.frame(genes=sample7_Monocytes_transcript$Gene,transcripts=sample7_Monocytes_transcript$ID)
sample7_Monocytes_transcript = sample7_Monocytes_transcript%>%dplyr::select(-ID,-Gene)%>%as.matrix()%>%t()

setwd("~/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/nanopore")
###################CD4 VS all##################
#human PBMC samples----gene level: CD4 VS all#
df_gene = scotch_gene(sample7_CD4_gene, rbind(sample7_B_gene,sample7_CD8_gene,sample7_NK_gene,sample7_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD4 VS CD8#
df_transcript = scotch_transcript(gene_transcript_CD4_df,gene_transcript_B_df, 
                                  sample7_CD4_transcript, rbind(sample7_B_transcript,sample7_CD8_transcript,sample7_NK_transcript,sample7_Monocytes_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
#write.csv(df_scotch,'nanopore_CD4vsAll_sample7.csv')

###################CD8 VS all##################
#human PBMC samples----gene level: CD8 VS all#
df_gene = scotch_gene(sample7_CD8_gene, rbind(sample7_B_gene,sample7_CD4_gene,sample7_NK_gene,sample7_Monocytes_gene), epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#human PBMC samples----transcript level: CD8 VS all#
df_transcript = scotch_transcript(gene_transcript_CD8_df,gene_transcript_B_df, sample7_CD8_transcript, 
                                  rbind(sample7_B_transcript,sample7_CD4_transcript,sample7_NK_transcript,sample7_Monocytes_transcript),
                                  ncores=10)

df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
#write.csv(df_scotch,'nanopore_CD8vsAll_sample7.csv')


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
#write.csv(df_scotch,'nanopore_BvsAll_sample7.csv')

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
#write.csv(df_scotch,'nanopore_NKvsAll_sample7.csv')

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
#write.csv(df_scotch,'nanopore_MonocytesvsAll_sample7.csv')




####################################################################
##################-----sample7 vs sample8----------#################
####################################################################
setwd("/Users/zhuoranx/Documents/ResearchProject/SingleCellLongReads/data/PBMC/results/nanopore")
###B 
sample7_B = read.csv('nanopore_BvsAll_sample7.csv')[,-1]
sample8_B = read.csv('nanopore_BvsAll_sample8.csv')[,-1]

genes = intersect(sample7_B$genes,sample8_B$genes)
B_sample7vs8 = full_join(sample7_B,sample8_B,by=c('genes','isoforms'))%>%
  arrange(genes)
write.csv(B_sample7vs8,"B_sample7vs8.csv")


###CD4
sample7_CD4 = read.csv('nanopore_CD4vsAll_sample7.csv')[,-1]
sample8_CD4 = read.csv('nanopore_CD4vsAll_sample8.csv')[,-1]
genes = intersect(sample7_CD4$genes,sample8_CD4$genes)
CD4_sample7vs8 = full_join(sample7_CD4,sample8_CD4,by=c('genes','isoforms'))%>%
  arrange(genes)
write.csv(CD4_sample7vs8,"CD4_sample7vs8.csv")

###CD8
sample7_CD8 = read.csv('nanopore_CD8vsAll_sample7.csv')[,-1]
sample8_CD8 = read.csv('nanopore_CD8vsAll_sample8.csv')[,-1]
genes = intersect(sample7_CD8$genes,sample8_CD8$genes)
CD8_sample7vs8 = full_join(sample7_CD8,sample8_CD8,by=c('genes','isoforms'))%>%
  arrange(genes)
write.csv(CD8_sample7vs8,"CD8_sample7vs8.csv")


###NK
sample7_NK = read.csv('nanopore_NKvsAll_sample7.csv')[,-1]
sample8_NK = read.csv('nanopore_NKvsAll_sample8.csv')[,-1]
genes = intersect(sample7_NK$genes,sample8_NK$genes)
NK_sample7vs8 = full_join(sample7_NK,sample8_NK,by=c('genes','isoforms'))%>%
  arrange(genes)
write.csv(NK_sample7vs8,"NK_sample7vs8.csv")

###Monocytes
sample7_Monocytes = read.csv('nanopore_MonocytesvsAll_sample7.csv')[,-1]
sample8_Monocytes = read.csv('nanopore_MonocytesvsAll_sample8.csv')[,-1]
genes = intersect(sample7_Monocytes$genes,sample8_Monocytes$genes)
Monocytes_sample7vs8 = full_join(sample7_Monocytes,sample8_Monocytes,by=c('genes','isoforms'))%>%
  arrange(genes)
write.csv(Monocytes_sample7vs8,"Monocytes_sample7vs8.csv")


