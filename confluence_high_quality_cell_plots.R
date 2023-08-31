for (row in 1:nrow(newtable)) {
  if (!is.na(newtable[row, "both"])) newcol[row] <- "pass both"
  else if (!is.na(newtable[row, "cr"]) && newtable[row, "mine"]==FALSE) {
      if (newtable[row, "percent.mt"] >15 && newtable[row, "nUMI"]>=999 && newtable[row, "nGene"]>=500 && newtable[row, "nGene"]<6000 && newtable[row, "nFrags"] > 300 && newtable[row, "TSS"]>4) newcol[row]<-"-mito"
      else if  (newtable[row, "nFrags"] > 300 && newtable[row, "TSS"]>4 &&  (newtable[row, "nUMI"] < 1000 || newtable[row, "nGene"]<500 || newtable[row, "nGene"]>6000)) newcol[row] <- "-other_rna"
      else if (newtable[row, "nFrags"] < 300 && newtable[row, "TSS"]>4 && newtable[row, "percent.mt"] <=15 && newtable[row, "nGene"]>=500 && newtable[row, "nGene"]<6000 && newtable[row, "nUMI"]>=999) newcol[row]<-"-frags"
      else if (newtable[row, "TSS"]<4  && newtable[row, "percent.mt"] <=15 && newtable[row, "nGene"]>=500 && newtable[row, "nGene"]<6000 && newtable[row, "nUMI"]>=999) newcol[row]<-"-TSS"
      else newcol[row]<-"-fail multi"
   }
   else newcol[row] <- "-fail all"
}
newtable$classification<-newcol
ggplot(newtable, aes(x = nFrags, y = nUMI, color = classification, order=nFrags)) + geom_point() + scale_y_log10() + scale_x_log10() + scale_color_viridis(discrete = TRUE, option = "D")

ggplot(newtable, aes(x = nFrags, y = nUMI, color = classification)) + geom_point() + scale_y_log10() + scale_x_log10() + scale_color_manual(values = c("#000000", "#ffe119", "#9A6324", "purple", "#469990", "orange","#f032e6"))

threshs<-c(0,0.0002,0.0005)
inds<-which(Multiome.rna2$p16.norm==threshs[1])
df1<-data.frame(p16.0=Multiome.rna2@meta.data$p16.norm[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p16.norm>threshs[1] & Multiome.rna2$p16.norm<=threshs[2])
df2<-data.frame(p16.0=Multiome.rna2@meta.data$p16.norm[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p16.norm>threshs[2] & Multiome.rna2$p16.norm<=threshs[3])
df3<-data.frame(p16.0=Multiome.rna2@meta.data$p16.norm[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p16.norm>threshs[3])
df4<-data.frame(p16.0=Multiome.rna2@meta.data$p16.norm[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
boxplot(df1$score,df2$score,df3$score,df4$score,notch=TRUE,names=c("p16 0", "p16 <0.0002", "p16 <0.0005", "p16 >0.0005"), ylab="Senescence score (from bulk Hayflick)", main="Senescence score by grouped p16 norm")



threshs<-c(0,3,10)
inds<-which(Multiome.rna2$p14.norm==threshs[1])
df1<-data.frame(p16.0=Multiome.rna2@meta.data$p14.norm[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p14.norm>threshs[1] & Multiome.rna2$p14.norm<=threshs[2])
df2<-data.frame(p16.0=Multiome.rna2@meta.data$p14.norm[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p14.norm>threshs[2] & Multiome.rna2$p14.norm<=threshs[3])
df3<-data.frame(p16.0=Multiome.rna2@meta.data$p14.norm[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p14.norm>threshs[3])
df4<-data.frame(p16.0=Multiome.rna2@meta.data$p14.norm[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
boxplot(df1$score,df2$score,df3$score,df4$score,notch=TRUE,names=c("p14 0", "p14 1-3", "p14 4-10", "p14 >10"), ylab="Senescence score (from bulk Hayflick)", main="Senescence score by grouped p14 dialout UMIs")




threshs<-c(0,3,6)
inds<-which(Multiome.rna2$p16.dialout==threshs[1])
df1<-data.frame(p16.0=Multiome.rna2@meta.data$p16.dialout[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p16.dialout>threshs[1] & Multiome.rna2$p16.dialout<=threshs[2])
df2<-data.frame(p16.0=Multiome.rna2@meta.data$p16.dialout[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p16.dialout>threshs[2] & Multiome.rna2$p16.dialout<=threshs[3])
df3<-data.frame(p16.0=Multiome.rna2@meta.data$p16.dialout[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p16.dialout>threshs[3])
df4<-data.frame(p16.0=Multiome.rna2@meta.data$p16.dialout[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
boxplot(df1$score,df2$score,df3$score,df4$score,notch=TRUE,names=c("p16 0", "p16 1-3", "p16 4-6", "p16 >6"), ylab="Senescence score (from bulk Hayflick)", main="Senescence score by grouped p16 dialout UMIs")



threshs<-c(0,3,10)
inds<-which(Multiome.rna2$p14.dialout==threshs[1])
df1<-data.frame(p16.0=Multiome.rna2@meta.data$p14.dialout[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p14.dialout>threshs[1] & Multiome.rna2$p14.dialout<=threshs[2])
df2<-data.frame(p16.0=Multiome.rna2@meta.data$p14.dialout[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p14.dialout>threshs[2] & Multiome.rna2$p14.dialout<=threshs[3])
df3<-data.frame(p16.0=Multiome.rna2@meta.data$p14.dialout[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
inds<-which(Multiome.rna2$p14.dialout>threshs[3])
df4<-data.frame(p16.0=Multiome.rna2@meta.data$p14.dialout[inds], score=Multiome.rna2@meta.data$Senescence.upregulated.score1[inds])
boxplot(df1$score,df2$score,df3$score,df4$score,notch=TRUE,names=c("p14 0", "p14 1-3", "p14 4-10", "p14 >10"), ylab="Senescence score (from bulk Hayflick)", main="Senescence score by grouped p14 dialout UMIs")
