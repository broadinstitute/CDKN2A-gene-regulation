library('ComplexHeatmap')
library('UpSetR')

a=read.delim('GTEx_heatmap_data.csv', sep=',')

a$p14_expr=ifelse(a$CDKN2A.p14 > 1, 1, 0)
a$p16_expr=ifelse(a$CDKN2A.p16 > 1, 1, 0)
a$p15_expr=ifelse(a$CDKN2B > 1, 1, 0)
a$ANRIL_expr=ifelse(a$CDKN2B.AS1 > 1, 1, 0)

ups=subset(a, select=c(p14_expr, p16_expr, p15_expr, ANRIL_expr))

m1=make_comb_mat(ups)

#pdf('GTEx_alltissues_upsetRplot.pdf', width = 7, height = 5.5)

UpSet(m1,  comb_order = order(-comb_size(m1)), top_annotation = upset_top_annotation(m1, height = unit(10, "cm"), add_numbers = TRUE))

#dev.off()