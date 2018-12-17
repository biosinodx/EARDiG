### GSED: Gene set enrichment in disease categories
# Copyright (C) 2018  Dong, Xiao
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.

# Xiao Dong, 2018.12.17, biosinodx@gmail.com, xiao.dong@einstein.yu.edu
# 2018.12.17 - v1.0.0 First publish version
# 2016.08.16 - v0.1.3 (Optional changes) trim off "Others" in disease categories from network output
# 2016.08.14 - v0.1.2 replace v0.1.2 overlaping by igraph output
# 2016.07.20 - v0.1.1 output overlapping genes with disease categories; correct naming of disease categories; add disease output
# 2016.06.09 - v0.1.1 add command line options
# 2016.03.16 - v0.1.0 initial version

# Usage
# Rscript EARDiG.R workingdir projectname inputgenelist.txt backgroundgenelist.txt repeattime diseaseclass
# workingdir: e.g. ./
# projectname: anyname, e.g. test
# inputgenelist.txt: file name and its root of an input gene list, should be tab limited file with header, the first colunmn is the gene names of interest.
# backgroundgenelist.txt: file name and its root of an input background gene list, should be tab limited file with header, the first colunmn is the gene names of interest, 2rd colunmn gene start, 3rd colunmn gene end; (e.g. By Xiao Dong, for protein coding genes, as in the depository "./background_genelist.txt")
# repeattime: recommend 2000
# diseaseclass: file name and its root of a predefined disease - gene classes; (e.g. by Simon C Johnson as in the depository "./diseasecat_simon_agingcell_2015.RData")

verion='v1.0.0'
print(paste("EARDiG: Enrichment in ARD-associated GWAS: , verion", verion))
Args <- commandArgs(TRUE)

# project name
#PN='iis'
PN=Args[2]
# input gene list, should be tab limited file with header, the first colunmn is the gene names of interest.
#genelist_candidate='iis_simon.txt.dedup'
genelist_candidate=Args[3]
# input background gene list, should be tab limited file with header, the first colunmn is the gene names of interest, 2rd colunmn gene start, 3rd colunmn gene end
#genelist_background='background_genelist.txt'
genelist_background=Args[4]
# input number of random gene set should be generated
#numrandom=2
numrandom=as.numeric(Args[5])
# load gene class assignment
# gdclass='diseasecat_simon_agingcell_2015.RData'
gdclass=Args[6]
# setworking dir
#wd='~/Desktop/2015-diseaseenrichment/workdir4'
wd=Args[1]

### analysis
setwd(wd)
dir.create(PN)
dir.create(paste(PN,'/randomset',sep=''))

genelist_background=read.table(genelist_background, header=T)
genelist_candidate=read.table(genelist_candidate, header=T)

# generate percentiles of gene length in background gene list, may take a few minutes
x=genelist_background
x$l=x[,3]- x[,2]
a=levels(x[,1])
b=vector()
for(i in 1:length(a)){
	b[i]=mean(x$l[x$Associated_Gene_Name==a[i]])
}
blist=data.frame(a,b)
colnames(blist)=c('gene','length')

q=quantile(blist$length, 1:100/100)
x=vector()
for(i in 1:nrow(blist)){
	x[i]=which(blist$length[i]<=q)[1]
}
blist$quantile = x
genelist=genelist_candidate
g = merge(genelist[,1], blist, by=1)

# generate random gene sets which saves to project dir, this may take an hour for 1000 random gene sets
tmp=as.factor(g$quantile)
tmp=summary(tmp, maxsum=nrow(g))
print('Generating random gene sets')

for(k in 1:numrandom){
	print(paste('Random set #', k, sep=''))
	randomdat=vector()
	for(i in 1:length(tmp)){
		x=blist[names(tmp)[i]==blist$quantile,]
		x=x[sample(nrow(x),tmp[i], replace=T),]
		randomdat=rbind(randomdat,x)
	}
	write.table(randomdat, paste('./',PN,'/randomset/randomgene_set', k, '.txt' , sep=''), col.names=T, row.names=F, sep='\t', quote=F)
}


# 
counting=function(c_genelist,c_categorylist){
	c_out=vector()
	for(i in 1:length(c_categorylist)){
		c_out[i] <- nrow(merge(c_genelist, c_categorylist[[i]], by=1))
	}
	names(c_out)=names(c_categorylist)
	return(c_out)
}

#overlap=function(c_genelist,c_categorylist){
#	c_out=vector()
#	for(i in 1:length(c_categorylist)){
#		a=as.vector(merge(c_genelist, c_categorylist[[i]], by=1)[,1])
#		if(length(a)>0)	c_out <- rbind(c_out, cbind(names(c_categorylist)[i], a))
#	}
#	colnames(c_out)=c('DiseaseCat','Gene')
#	return(c_out)
#}

load(gdclass)

#tmp=overlap(g, diseasecat_gene)
#tmp2=vector()
#for(i in 1:length(disease_gene)){
#	tmp2=rbind(tmp2, cbind(names(disease_gene)[i], disease_gene[[i]]) )
#}
#tmp3=unique(merge(tmp, tmp2, by=2))
#colnames(tmp3)=c('QueryGene','DiseaseCat','Disease')
#write.table(tmp3, paste('./',PN,'/overlaps.txt',sep=''), col.names=T, row.names=F, sep='\t', quote=F)

datout=vector()
datout=rbind(datout, counting(g, diseasecat_gene))
rownames(datout)='observed'
randomout=vector()
for(i in 1:numrandom){
	x=read.table(paste('./',PN,'/randomset/randomgene_set', i, '.txt' , sep=''), header=T)
	randomout=rbind(randomout, counting(as.vector(x[,1]), diseasecat_gene))
}
rownames(randomout)=paste('randomset',1:numrandom,sep='')
pvalue=vector()
for(i in 1:length(diseasecat_gene)){
	pvalue[i]=(length(which(datout[1,i] <= randomout[,i]))+1)/(numrandom+1)
}

fileout=rbind(datout, pvalue, randomout)
write.table(fileout,paste('./',PN,'/results.txt',sep=''), col.names=T, row.names=T, sep='\t', quote=F)

for(i in 1:length(diseasecat_gene)){
	pdf(file=paste('./',PN,'/DensityPlot_',names(diseasecat_gene[i]),'.pdf', sep=''))
	plot(density(randomout[,i]), main=paste(names(diseasecat_gene[i]),', pvalue=', round(pvalue[i],4),sep=''),xlab='# genes (distribution - null hypothesis; red dashed line - observed)', xlim=c(0,max(c(randomout[,i], datout[,i]))))
	abline(v=datout[,i], col='red', cex=100, lty=2)
	dev.off()
}

require(igraph)

tmp=vector()
namecat=c('Others','Cancer','Cardiovascular','Other age-related','Metabolic','Neurodegenerative')
for(i in 1:length(diseasecat_disease)) tmp=c(tmp, rep(namecat[i],length(diseasecat_disease[[i]])))
t1=cbind(tmp, as.vector(unlist(diseasecat_disease)))
#
t1=t1[t1[,1]!='Others',]
#

tmp=vector()
namecat=names(disease_gene)
for(i in 1:length(disease_gene)) tmp=c(tmp, rep(namecat[i],length(disease_gene[[i]])))
t2=cbind(tmp, as.vector(unlist(disease_gene)))
#
pick=vector()
for(i in 1:nrow(t2)) if(length(which(as.vector(t2[i,1])==t1[,2]))>=1) pick = c(pick, i)
#
t2=t2[pick,]

t2tmp=merge(genelist_candidate, t2, by.x=1, by.y=2)
t2tmp=unique(t2tmp[,c(2,1)])
colnames(t2tmp)=c('a','b')

t1tmp=merge(as.vector(t2tmp[,1]), t1, by.x=1, by.y=2)
t1tmp=unique(t1tmp[,c(2,1)])
colnames(t1tmp)=c('a','b')

t0tmp=cbind('DiseaseCat',unique(as.vector(t1tmp[,1])))

tmp=rbind(t0tmp, as.matrix(t1tmp, ncol=2), as.matrix(t2tmp, ncol=2))

g=graph.data.frame(tmp,directed=TRUE)

g_font=vector()
g_cex=vector()
g_color=vector()
tmp=V(g)$name
for(i in 1:length(tmp)){
	if(tmp[i]=='DiseaseCat'){
		g_font[i]=2
		g_cex[i]=1
		g_color[i]='black'
	}else if(length(which(tmp[i]==t0tmp[,2]))>=1){
		g_font[i]=2
		g_cex[i]=0.8
		g_color[i]='black'
	}else if(length(which(tmp[i]==t1tmp[,2]))>=1){
		g_font[i]=1
		g_cex[i]=0.6
		g_color[i]='black'
	}else if(length(which(tmp[i]==t2tmp[,2]))>=1){
		g_font[i]=3
		g_cex[i]=0.4
		g_color[i]='black'
	}
}
g_font=g_font[!is.na(g_font)]
g_cex=g_cex[!is.na(g_cex)]
g_color=g_color[!is.na(g_color)]

V(g)$size=0.01
#plot(g, vertex.label.color= "black",vertex.label.font=3,vertex.label.cex=1, edge.arrow.size=1.2, vertex.shape="circle",edge.color="light grey", vertex.color="red", asp=0)
pdf(paste('./',PN,'/overlaps.pdf',sep=''),width=8,height=8)
plot(g, vertex.label.color= g_color,vertex.label.font=g_font,vertex.label.cex=g_cex, edge.arrow.size=0.5, vertex.shape="circle",edge.color="light grey", vertex.color="red", asp=0)
dev.off()

#plot(g, layout =  layout.kamada.kawai,  vertex.label = V(g)$name,  vertex.label.color= "white",edge.arrow.size=0.5,  edge.curved=T, edge.label=E(g)$Freq, edge.label.color="pink", edge.label.font=5,vertex.shape="circle",edge.color="white", vertex.color="red", asp=0)
#title("This is my first igraph",cex.main=3,col.main="green")

#dtree <- dominator.tree(g, root='DiseaseCat')
#layout <- layout_as_tree(dtree$domtree, root='DiseaseCat')
#layout[,2] <- -layout[,2]
#plot(dtree$domtree, layout=layout, vertex.label=V(dtree$domtree)$name)
