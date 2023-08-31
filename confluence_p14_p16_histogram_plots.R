df = data.frame(p14=Multiome.rna$p14,p16=Multiome.rna$p16,CDKN2A=tmp)
tmp <- FindMarkers(Multiome.rna, features = "CDKN2A",ident.1 = "CDKN2A")
t <- (Multiome.rna$p14 >= 5 & Multiome.rna$p16 < 5 )
u <- (Multiome.rna$p16 >= 5 & Multiome.rna$p14 < 5 )
w <- (Multiome.rna$p16 >= 5 & Multiome.rna$p14 >= 5)
sum(t|u)*100/(sum(t|u)+sum(w))

x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
for (val in x) {
t <- (p14.umis >= val & p16.umis < val )
u <- (p16.umis >= val & p14.umis < val )
w <- (p16.umis >= val & p14.umis >= val)
tot <- (p14.umis >= val | p16.umis >= val)
print(c(val,sum(t|u)*100/(sum(t|u)+sum(w)), sum(tot)))
}

p16.umis <- Multiome.rna@meta.data$p16.dialout
p14.umis <- Multiome.rna@meta.data$p14.dialout
df <- data.frame(p16.umis, p14.umis)

ggplot(df,aes(p16.umis))+geom_histogram(aes(y = ..density..))+facet_grid(~p14.umis)+xlim(0,8)
ggplot(df,aes(p14.umis))+geom_histogram()+facet_grid(~p16.umis,rows = 12,cols=1)
ggplot(df,aes(p14.umis))+geom_histogram()+facet_grid(~p16.umis,rows = p16.umis)

ggplot(df,aes(p14))+geom_histogram()+facet_grid(rows=vars(p16),scales="free_y")

ggplot(df,aes(p16.umis))+geom_histogram()+facet_grid(rows=vars(p14.umis)) + ylim(0,10)
