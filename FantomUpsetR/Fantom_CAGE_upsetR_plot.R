library('ComplexHeatmap')
library('UpSetR')

a=read.delim('hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.averaged_tech_rep.averaging_donors.txt', sep='\t')

a$ANRIL_expr=ifelse(a$ANRIL>1, 1, 0)
a$p16_expr=ifelse(a$p16>1, 1, 0)
a$p14_expr=ifelse(a$p14>1, 1, 0)
a$p15_expr=ifelse(a$p15>1, 1, 0)

ups=subset(a, select=c(p14_expr, p16_expr, p15_expr, ANRIL_expr))

m1=make_comb_mat(ups)

pdf('Fantom_CAGE_onlyAdult_avgtechrepl_avgdonors_upsetRplot.pdf', width = 7, height = 5.5)

UpSet(m1,  comb_order = order(-comb_size(m1)), top_annotation = upset_top_annotation(m1, height = unit(10, "cm"), add_numbers = TRUE))

dev.off()