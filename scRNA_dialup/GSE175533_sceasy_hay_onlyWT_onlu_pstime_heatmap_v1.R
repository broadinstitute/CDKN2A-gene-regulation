setwd('/Volumes/etorlait/p16/Published_Data/Chan_et_al_2022_eLife/scRNAseq/')

a=read.delim('GSE175533_sceasy_hay_onlyWT_redoUMAP_onlu_pstime_table_for_heatmap_newSignScores.txt', sep='\t')

library(ComplexHeatmap)
library(circlize)
library(magick)

a2 <- a[order(a$pstime),]

#PDL, Mki67 for cell cycle, p21?CDKN1A, TGFbeta, p14ARF signature, p16INK4A signature, p15INK4B signature 

#PDL and pstime in annotation, 

#Mki67, p14 sign, p16 sign, p15 sign, CKDN1A, TGFb


set.seed(3)

try=a2[1:5,]

col_fun = colorRamp2(c(-3, 0, 5,18), c('grey', 'white', '#cb181d', '#99000d'))

col_fun2 = colorRamp2(c(0, max(a2$pstime)), c('#deebf7','#0c2c84'))
#col_fun2 = magma(256)

column_ha = HeatmapAnnotation(pstime=a2$pstime,
                              cellcycle=a2$Phase,
                              PDL = a2$PDL, 
                              col=list(pstime=col_fun2,
                                       cellcycle=c('G1'='#66c2a5',
                                             'S'='#8da0cb', 
                                             'G2M'='#fc8d62'),
                                       PDL=c('PDL_25'='#C15D80',
                                             'PDL_29'='#D191A7', 
                                             'PDL_33'='#E1CAD1', 
                                             'PDL_37'='#D5D2E9', 
                                             'PDL_46'='#A9A9D8', 
                                             'PDL_50'='#747DBF'), 
                                       na_col="white"
                                       ))

to_plot=as.matrix(subset(a2, select=c('MKI67', 'high_corr_p14_norm_genes','high_corr_p16_norm_genes', 'high_corr_CDKN2B_genes', 'CDKN1A', 'senmayo', 'sasp', 'TGFB1','CDKN2B', 'p16.norm', 'p14.norm', 'CDKN2A')))

to_plot_scaled = scale(to_plot)

pdf('GSE175533_sceasy_hay_onlyWT_redoUMAP_onlu_pstime_heatmap_wCDKN2A_newSignScores_noRaster.pdf', width = 12, height =3 )
Heatmap(t(to_plot_scaled), 
        top_annotation = column_ha, cluster_rows = F, cluster_columns = F, show_column_names = F, 
        col = col_fun, name = 'Z-score', use_raster = F) 
dev.off()

