library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

set.seed(123)

## load cell type annotations (see Methods for how they are generated)

cells=read.delim('1E_Cell_Counts_Tissues_Simplified_MajorTissue_withCellType.txt',  sep='\t', header = T)
cells=subset(cells, select=c('Cell_ID2', 'Label_Fig1E', 'tissue_label' ))
colnames(cells)=c("cell.name", 'cell_group', 'main_tissue')

cell_ids=read.delim('Zhang_celltypes_withIDs.txt',  sep='\t',header = F)
colnames(cell_ids)=c("cell.name", 'cell.id')
cell_ids=subset(cell_ids, cell.id <= 111)

cells=merge(cells, cell_ids, by='cell.name')
cells=cells[order(cells$cell.id),]

### matrix.mtx.gz: matrix in matrix market format, each line represents a cCRE-cell type pair. For example, "2 22" means cCRE No. 2 is accessibile in cell type No. 22.
### from Zhang,Hocker et al, only regions overlapping 9p21 promoters and CREs 

a=read.delim('matrix_only_overlapping_9p21_f01_regions.txt', sep='\t', header = F)
colnames(a)=c('cre.id', "cell.id")

##keep only "adult" cell types
a=subset(a, cell.id <=111)

##mark those as open region
a$status=1

## get the cell types that are in the a matrix

selected_cell_ids=sort(unique(a$cell.id))
length(selected_cell_ids)
#all 111 cell types are present at least once

wide_a= a %>%
  pivot_wider(names_from = cell.id, values_from = status, values_fill = 0)

row.names(wide_a)=wide_a$cre.id

matrix_a=as.matrix(wide_a)
matrix_a<-matrix_a[,-1]
matrix_a1=matrix_a[order(as.numeric(row.names(matrix_a))),order(as.numeric(colnames(matrix_a)))]

colnames(matrix_a1)=cells$cell.name

col_fun <- c("0" = "#f0f0f0", "1" = "#3182bd")

ha = columnAnnotation(group=cells$cell_group, organ=cells$main_tissue, 
                    col=list(group=c('AdrCor'='#A6CEE3',
                                    'BrainGlia'='#C51B7D',
                                    'BrainNeur'='#7D54A5',
                                    'Endoth'='#EAD27A',
                                    'GI_Epith'='#005b4d',
                                    'Other_Epith'='#4DAF4A',
                                    'ImmuneMyelo'='red4',
                                    'ImmuneLymp'='#FB6A4A',
                                    'Islet'='#b8e186',
                                    'Mural'='#B9A499',
                                    'Myocyte'='thistle2',
                                    'SmMusc'='#2b3990',
                                   'Stromal'='#B15928'),
                             
                             organ=c('AdipOm'='#fcf75e',
                                     'Adrenal'='#33a02c',
                                     'Brain'='#900c3c',
                                     'Breast'='#bcbddc',
                                     'Cardiovasc'='#e31a1c',
                                     'FemRepr'='#54278f',
                                     'GI'='#fecdc8',
                                     'Liver'='#8d4505',
                                     'Lung'='#b2df8a',
                                     'Multiple'='#dadada',
                                     'Nerve'='steelblue',
                                     'Pancreas'='#66c2a5',
                                     'Skin'='#fdb863',
                                     'SkMuscl'='#016c59',
                                     'Thyroid'='#c98a3a')), 
                    annotation_name_side= 'left')

## Load up annotation on which region regulates which transcript 

valid=read.delim('Gene_Region_Validation_Matrix_withZhangIDs.txt')

row_groups <- c('g1', 
                       'g2','g2','g2','g2','g2','g2',
                       'g3','g3','g3','g3','g3','g3','g3','g3')


ha_tp = rowAnnotation(
  p16=valid$p16_on,
  p14=valid$p14_on,
  ANRIL=valid$ANRIL_on,
  p15=valid$p15_on,
  col = list(
             p15 = c('1'='red', 
                            '0'='white'), 
             ANRIL = c('1'='purple', 
                              '0'='white'),
             p14 = c('1'='green', 
                     '0'='white'), 
             p16 = c('1'='blue', 
                     '0'='white')),
  annotation_name_gp= gpar(fontsize = 9),
  annotation_name_side= 'top')

# get dendrogram separately for customisation
col_dend <- as.dendrogram(hclust(dist(t(matrix_a1)), method = 'complete'))

pdf('7D_CRE_heatmap.pdf', width = 20, height = 7)
Heatmap(matrix_a1, cluster_columns = col_dend, cluster_rows = F, col=col_fun, row_title='CRE annotation (WI-38)', 
        column_title = 'Cell types', name = '0:closed\n1:open', right_annotation = ha_tp, 
        row_split=row_groups,
        bottom_annotation = ha, 
        show_row_names = T, row_names_side = "left", row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 9), column_names_side = 'bottom')
dev.off()
    

