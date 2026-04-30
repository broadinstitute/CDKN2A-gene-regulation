
a=read.delim('heatmap_data_Brugge.csv', sep=',')
rownames(a)=a$Tissue

a$ratio=a$CDKN2A.p16/a$CDKN2A.p14

a_ordered <- a[order(a$ratio, decreasing = T),]

a1=subset(a_ordered, select=c('CDKN2A.p14', 'CDKN2A.p16'))

a2=subset(a_ordered, select=c('CDKN2A', 'CDKN2B', 'CDKN2B.AS1', 'MTAP'))

library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(0, 0.95,1, 65), c("#af8dc3",'white',"#fdbb84", "#e6550d"))

pdf('GTEx_heatmap_CDKN2A_CDKN2B_CDKN2BAS1_GTEx_RawExpression_rev_order_all_tissue_v2.pdf', width = 5, height = 20)

#change order 

Heatmap(a2, cluster_rows = F, cluster_columns = F, col=col_fun,name = "tpm",
        cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", a2[i, j]), x, y, gp = gpar(fontsize = 10), rot=90)
})

dev.off()

a1_scaled = t(scale(t(a1)))

col_fun2 = colorRamp2(c(0, 0.95, 1, 10), c("#a6bddb", 'white',"#fee0d2", "#de2d26"))

pdf('GTEx_heatmap_p14_p16_GTEx_RawExpression_rev_order_all_tissue_v2.pdf', width = 4, height = 20)

Heatmap(a1, cluster_rows = F, cluster_columns = F, col=col_fun2,name = "tpm", 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", a1[i, j]), x, y, gp = gpar(fontsize = 10), rot=90)
        })

dev.off()

## discuss which tissue to add -- this is what we had before
selected_tissues=c("Adrenal Gland", "Artery - Tibial", "Artery - Coronary", "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)", "Artery - Aorta", "Adipose - Subcutaneous", "Nerve - Tibial", "Lung", "Pituitary", "Cells - Cultured fibroblasts", "Colon - Transverse", "Cells - EBV-transformed lymphocytes", "Esophagus - Mucosa", "Vagina", "Esophagus - Muscularis", "Spleen", "Small Intestine - Terminal Ileum", "Whole Blood", "Testis")

a_sel=subset(a, Tissue %in% selected_tissues)

a1_sel=subset(a_sel, select=c('CDKN2A.p14', 'CDKN2A.p16'))

a2_sel=subset(a_sel, select=c('CDKN2A', 'CDKN2B', 'CDKN2B.AS1', 'MTAP'))

