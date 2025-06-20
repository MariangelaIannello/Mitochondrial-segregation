#!/usr/bin/env Rscript

tab <- read.table(file="Results_summary.tsv",header=F,col.names=c("OG","lnL_null","lnL_alt","Delta_w"))

for (i in tab[,1]) {
	tab[tab$OG==i,5] <- -2*(as.numeric(tab[tab$OG==i,2])-as.numeric(tab[tab$OG==i,3]))
	tab[tab$OG==i,6] <- pchisq(tab[tab$OG==i,5],df=1,lower.tail=FALSE)
}

names(tab)[5] <- "LRT"
names(tab)[6] <- "pvalue"
tab$BH_fdr <- p.adjust(tab$pvalue,method="BH")

write.table(tab,file="Results_summary_LRT.tsv",quote=F,row.names=F,col.names=T,sep="\t")


