rm(list=ls())
gc()
library("edgeR")
library("gplots")
library("RColorBrewer")
library(ggplot2)
library(matrixStats)

# This file contains the code to produce the replicate specific gene expression analysis
# The beginning of the code is strictly identical to the first analysis

#Read the merged count table per replicate
counts_merged=read.table("Counts_merged.txt",header=TRUE,sep='\t')

Pop_merged=substring(names(counts_merged),1,1)
Repl_merged=substring(names(counts_merged),6,7)
Temp_merged=substring(names(counts_merged),3,4)

c_B23_m=c(1:dim(counts_merged)[2])[Pop_merged =="B" & Temp_merged =="23"][1]
c_C23_m=c(1:dim(counts_merged)[2])[Pop_merged =="C" & Temp_merged =="23"][1]
c_H23_m=c(1:dim(counts_merged)[2])[Pop_merged =="H" & Temp_merged =="23"][1]

c_B15_m=c(1:dim(counts_merged)[2])[Pop_merged =="B" & Temp_merged =="15"][1]
c_C15_m=c(1:dim(counts_merged)[2])[Pop_merged =="C" & Temp_merged =="15"][1]
c_H15_m=c(1:dim(counts_merged)[2])[Pop_merged =="H" & Temp_merged =="15"][1]

## We retain genes with more than 1 counts per million reads on average > 0

keep=log2(rowSums(cpm(counts_merged,normalized.lib.sizes=TRUE))/dim(counts_merged)[2])>0
summary(keep)
expressed_genes <- row.names(counts_merged)[keep]

## Number of analyzed genes
# summary(keep)
#   Mode   FALSE    TRUE 
#logical    2062   11200

#############################################################
#### Now we will analyse the replicate change separately ####
#############################################################

# We restric to the ancestral population ("B") and the Hot evolved ones ("H")
counts_BH <- counts_merged[keep,Pop_merged%in%c("B","H")]

DGEList=DGEList(counts_BH)
DGEList=calcNormFactors(DGEList, method=c("TMM"))

Colors_temp=rep(c("chartreuse3","cornflowerblue","firebrick3"),each=2)

Pop_sp <- substring(names(counts_BH),1,1)
Temp_sp <- substring(names(counts_BH),3,4)
Repl_sp <- substring(names(counts_merged),6,7)

Pop_sp[Pop_sp=="H"] <- paste0(Pop_sp, Repl_sp)[Pop_sp=="H"]

design_plasticity=model.matrix(~ Pop_sp* Temp_sp)

####### Now we run the model

BCH_DGEList_GLM <- estimateGLMRobustDisp(DGEList, design_plasticity,verbose=TRUE,maxit=10)

BCH_fit_DGEList_GLM=glmFit(BCH_DGEList_GLM, design=design_plasticity)

### Build the contrasts and save the results in table
design <- design_plasticity
L1 = design[Pop_sp=="B" & Temp_sp=="23",][1,]-design[Pop_sp=="B" & Temp_sp=="15",][1,]

L2 = design[Pop_sp=="H11" & Temp_sp=="23",]-design[Pop_sp=="H11" & Temp_sp=="15",]
L3 = design[Pop_sp=="H12" & Temp_sp=="23",]-design[Pop_sp=="H12" & Temp_sp=="15",]
L4 = design[Pop_sp=="H13" & Temp_sp=="23",]-design[Pop_sp=="H13" & Temp_sp=="15",]
L5 = design[Pop_sp=="H14" & Temp_sp=="23",]-design[Pop_sp=="H14" & Temp_sp=="15",]
L6 = design[Pop_sp=="H15" & Temp_sp=="23",]-design[Pop_sp=="H15" & Temp_sp=="15",]

L7 <- L2-L1
L8 <- L3-L1
L9 <- L4-L1
L10 <- L5-L1
L11 <- L6-L1

L12 = design[Pop_sp=="H11" & Temp_sp=="23",]-design[Pop_sp=="B" & Temp_sp=="23",][1,]
L13 = design[Pop_sp=="H12" & Temp_sp=="23",]-design[Pop_sp=="B" & Temp_sp=="23",][1,]
L14 = design[Pop_sp=="H13" & Temp_sp=="23",]-design[Pop_sp=="B" & Temp_sp=="23",][1,]
L15 = design[Pop_sp=="H14" & Temp_sp=="23",]-design[Pop_sp=="B" & Temp_sp=="23",][1,]
L16 = design[Pop_sp=="H15" & Temp_sp=="23",]-design[Pop_sp=="B" & Temp_sp=="23",][1,]

L17 = design[Pop_sp=="H11" & Temp_sp=="15",]-design[Pop_sp=="B" & Temp_sp=="15",][1,]
L18 = design[Pop_sp=="H12" & Temp_sp=="15",]-design[Pop_sp=="B" & Temp_sp=="15",][1,]
L19 = design[Pop_sp=="H13" & Temp_sp=="15",]-design[Pop_sp=="B" & Temp_sp=="15",][1,]
L20 = design[Pop_sp=="H14" & Temp_sp=="15",]-design[Pop_sp=="B" & Temp_sp=="15",][1,]
L21 = design[Pop_sp=="H15" & Temp_sp=="15",]-design[Pop_sp=="B" & Temp_sp=="15",][1,]

L <- cbind(L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15,L16,L17,L18,L19,L20,L21)

LRT_fits <- list()
LRT_fits_DE <- list()
LRT_fits_DE_low <- list()
for(k in 1: ncol(L)){
	LRT_fits[[k]] <- glmLRT(BCH_fit_DGEList_GLM, contrast=L[,k])	
} 

for(k in 1:ncol(L)) LRT_fits_DE[[k]] <- decideTestsDGE(LRT_fits[[k]], p=0.05, adjust="BH")
for(k in 1:ncol(L)) LRT_fits_DE_low[[k]] <- decideTestsDGE(LRT_fits[[k]], p=0.1, adjust="BH")

## Plastic genes in the base :
table(LRT_fits_DE[[1]])
table(LRT_fits_DE[[2]]) 
table(LRT_fits_DE[[3]]) 

###
library(VennDiagram)
vec_Venn=cbind(LRT_fits_DE_low[[7]],
LRT_fits_DE_low[[8]],
LRT_fits_DE_low[[9]],
LRT_fits_DE_low[[10]],
LRT_fits_DE_low[[11]])

vec_Venn_stringent=cbind(LRT_fits_DE[[7]],
LRT_fits_DE[[8]],
LRT_fits_DE[[9]],
LRT_fits_DE[[10]],
LRT_fits_DE[[11]])

a=vennCounts(vec_Venn)
vennDiagram(a,names=paste("R",1:5))
table(rowSums(vec_Venn))

sign_all5 <-(rowSums(abs(vec_Venn)))==5
data.frame(expressed_genes[sign_all5])

sign_all5 <-(rowSums(abs(vec_Venn_stringent)))==5
data.frame(expressed_genes[sign_all5])

sign_all1 <-(rowSums(abs(vec_Venn_stringent)))==1
sum(sign_all1) # 272

summary(colSums(vec_Venn_stringent!=0))
sum(rowSums(vec_Venn_stringent!=0)>=1) # 409 genes

###############################################
#Inc/dec plasticity for each of the replicate #
###############################################

inc_plast <- NULL
for(i in 1:5){
inc_plast <- cbind(inc_plast,LRT_fits_DE_low[[6+i]]!=0 & (LRT_fits_DE[[11+i]]!=0 | LRT_fits_DE[[16+i]]!=0) & abs(LRT_fits[[1]]$table[,1])<abs(LRT_fits[[1+i]]$table[,1]) )
}

dec_plast <- NULL
for(i in 1:5){
dec_plast <- cbind(dec_plast,LRT_fits_DE_low[[6+i]]!=0 & (LRT_fits_DE[[11+i]]!=0 | LRT_fits_DE[[16+i]]!=0) & abs(LRT_fits[[1]]$table[,1])>abs(LRT_fits[[1+i]]$table[,1]) )
}

sort(colSums(inc_plast)) 
sort(colSums(dec_plast))

sum(rowSums(inc_plast)>=1) # 318 genes
sum(rowSums(dec_plast)>=1) # 101 genes

table(rowSums(inc_plast))
table(rowSums(dec_plast))

#### Proportion of shared genes evolving plasticity

shared_in=data.frame(rep=1:5,nb_sign=as.numeric(colSums(inc_plast)),at_least2=c(
sum(inc_plast[,1] & rowSums(inc_plast[,1:5])>=3),
sum(inc_plast[,2] & rowSums(inc_plast[,1:5])>=3),
sum(inc_plast[,3] & rowSums(inc_plast[,1:5])>=3),
sum(inc_plast[,4] & rowSums(inc_plast[,1:5])>=3),
sum(inc_plast[,5] & rowSums(inc_plast[,1:5])>=3)
))
shared_in$ratio = shared_in$at_least2/shared_in$nb_sign
summary(shared_in$ratio)

shared_de=data.frame(rep=1:5,nb_sign=as.numeric(colSums(dec_plast)),at_least2=c(
sum(dec_plast[,1] & rowSums(dec_plast[,1:5])>=3),
sum(dec_plast[,2] & rowSums(dec_plast[,1:5])>=3),
sum(dec_plast[,3] & rowSums(dec_plast[,1:5])>=3),
sum(dec_plast[,4] & rowSums(dec_plast[,1:5])>=3),
sum(dec_plast[,5] & rowSums(dec_plast[,1:5])>=3)
))
shared_de$ratio = shared_de$at_least2/shared_de$nb_sign
summary(shared_de$ratio)

## Changes in reaction norms

list_inc_plast_23 <- list()
for(i in 1:5){
vect_table=NULL
test23 <-  inc_plast[,i] & LRT_fits_DE[[11+i]]!=0
print(sum(test23))

for(j in 1:5){
	if(j!=i){
		vect_table <- cbind(vect_table,LRT_fits[[6+j]]$table[test23,1])
}}
list_inc_plast_23[[i]] <- vect_table

}

rsq <- function(x, y) summary(lm(y~x))$r.squared

list_cor <- NULL
list_rsq <- NULL
for(j in 1:5){
list_cor_t <- NULL
list_rsq_t <- NULL

test23 <-  inc_plast[,j] & LRT_fits_DE[[11+j]]!=0
for(i in 1:4) list_cor_t <- c(list_cor_t ,cor(list_inc_plast_23[[j]][,i],LRT_fits[[1+j]]$table[test23,1]))
for(i in 1:4) list_rsq_t <- c(list_rsq_t ,rsq(list_inc_plast_23[[j]][,i],LRT_fits[[1+j]]$table[test23,1]))

list_cor=rbind(list_cor, list_cor_t)
list_rsq=rbind(list_cor, list_rsq_t)
}
mean(list_cor) # .81
mean(list_rsq) # .77

##############
list_dec_plast_23 <- list()
for(i in 1:5){
vect_table=NULL
test23 <-  dec_plast[,i] & LRT_fits_DE[[11+i]]!=0
print(sum(test23))

for(j in 1:5){
	if(j!=i){
		vect_table <- cbind(vect_table,LRT_fits[[6+j]]$table[test23,1])
}}
list_dec_plast_23[[i]] <- vect_table

}

list_cor2 <- NULL
list_rsq2 <- NULL

for(j in 1:5){
list_cor_t <- NULL
list_rsq_t <- NULL
test23 <-  dec_plast[,j] & LRT_fits_DE[[11+j]]!=0
for(i in 1:4) list_cor_t <- c(list_cor_t ,cor(list_dec_plast_23[[j]][,i],LRT_fits[[1+j]]$table[test23,1]))
for(i in 1:4) list_rsq_t <- c(list_rsq_t ,rsq(list_dec_plast_23[[j]][,i],LRT_fits[[1+j]]$table[test23,1]))
list_cor2 =rbind(list_cor2, list_cor_t)
list_rsq2 =rbind(list_rsq2, list_rsq_t)
}

mean(list_cor)
mean(list_cor2)
mean(list_rsq)
mean(list_rsq2)

###################################
##### Replicate specific genes ####
#####   for GO enrichment      ####
###################################

sign_all_sp <- (rowSums(vec_Venn)) <= c(-1)
sum(sign_all_sp) # 295 genes

# All replicates
write.table(expressed_genes[sign_all_sp],file="Rep_sp_Down.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(expressed_genes[sign_all_sp],file="Rep_sp_Down_stringent.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Replicate specific

for(i in 1:5) write.table(expressed_genes[inc_plast[,i]],file=paste0("Rep_sp_Inc_Plast",i,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)

for(i in 1:5) write.table(expressed_genes[dec_plast[,i]],file=paste0("Rep_sp_Dec_Plast",i,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)

## Load the table with the concatenated GO
GO_by_rep <- read.table("GO_by_rep_Inc.txt")
GO_df<-data.frame(GO=unique(GO_by_rep$V2))
table(table(GO_by_rep$V2))
(table(GO_by_rep$V1))

for(i in 1:5){	
	GO_df <- cbind(GO_df,GO_df$GO%in% subset(GO_by_rep,V1==paste0("R",i))$V2)
}
names(GO_df)[2:6] <- paste0("R",1:5)
head(GO_df)
data.frame(GO_df$GO[rowSums(GO_df[2:6])==5])

############################
########  Figures   ########
############################

# Figure 5B

ox_pho_genes=unique(read.table("GO_Oxidative_phosphorylation.txt",sep="\t")$V1)
length(ox_pho_genes) # 188

ox_pho=cpm(counts_BH,normalized.lib.sizes=TRUE)[expressed_genes%in% ox_pho_genes,]
for(i in 1:dim(ox_pho)[1]) ox_pho[i,]=(ox_pho[i,]-mean(ox_pho[i,Pop_sp=="B"])) / sd(ox_pho[i,Pop_sp=="B"])

Base_plast <-rowMeans(ox_pho[,Pop_sp=="B" & Temp_sp=="23"])-rowMeans(ox_pho[,Pop_sp=="B" & Temp_sp=="15"])
H1_plast <-(ox_pho[,Pop_sp=="H11" & Temp_sp=="23"])-(ox_pho[,Pop_sp=="H11" & Temp_sp=="15"])
H2_plast <-(ox_pho[,Pop_sp=="H12" & Temp_sp=="23"])-(ox_pho[,Pop_sp=="H12" & Temp_sp=="15"])
H3_plast <-(ox_pho[,Pop_sp=="H13" & Temp_sp=="23"])-(ox_pho[,Pop_sp=="H13" & Temp_sp=="15"])
H4_plast <-(ox_pho[,Pop_sp=="H14" & Temp_sp=="23"])-(ox_pho[,Pop_sp=="H14" & Temp_sp=="15"])
H5_plast <-(ox_pho[,Pop_sp=="H15" & Temp_sp=="23"])-(ox_pho[,Pop_sp=="H15" & Temp_sp=="15"])

df.violin <- data.frame(plast=c(Base_plast[Base_plast<0], H1_plast[Base_plast<0], H2_plast[Base_plast<0], H3_plast[Base_plast<0], H4_plast[Base_plast<0], H5_plast[Base_plast<0]),Pop=rep(c("Ancestral","Hot-1","Hot-2","Hot-3","Hot-4","Hot-5"),each=sum(Base_plast<0)))

ggplot(df.violin, aes(x=Pop, y=plast,fill=Pop)) + 
  geom_violin(trim=TRUE)+scale_fill_manual(values=rep(c("chartreuse3","firebrick3"),c(1,5)))+geom_boxplot(width=0.1)



##############################
######### Figure 5A ##########
##############################

genes=read.table("genesID2.txt",sep="\t")

glyc=c("Hex-A","Hex-C","Pgi","Pfk","fbp","Tpi","Ald","Gapdh1","Gapdh2","Pgk","Pglym78","Eno","PyK","ImpL3","CG11876","CG5261","CG7430","l(1)G0334")

gg=subset(genes,V4%in%glyc)
gg=gg[order(gg$V4),]
gg=gg[order(order(glyc)),]

gly_exp=cpm(counts_BH,normalized.lib.sizes=TRUE)[row.names(counts_BH)%in%gg$V1,]
for(i in 1:dim(gly_exp)[1]) gly_exp[i,]=(gly_exp[i,]-mean(gly_exp[i,Pop_sp=="B"])) / sd(gly_exp[i,Pop_sp=="B"])
dim(gly_exp)

Base_plast <- rowMeans(gly_exp[,Pop_sp=="B" & Temp_sp=="23"])-rowMeans(gly_exp[,Pop_sp=="B" & Temp_sp=="15"])
H1_plast <-(gly_exp[,Pop_sp=="H11" & Temp_sp=="23"])-(gly_exp[,Pop_sp=="H11" & Temp_sp=="15"])
H2_plast <-(gly_exp[,Pop_sp=="H12" & Temp_sp=="23"])-(gly_exp[,Pop_sp=="H12" & Temp_sp=="15"])
H3_plast <-(gly_exp[,Pop_sp=="H13" & Temp_sp=="23"])-(gly_exp[,Pop_sp=="H13" & Temp_sp=="15"])
H4_plast <-(gly_exp[,Pop_sp=="H14" & Temp_sp=="23"])-(gly_exp[,Pop_sp=="H14" & Temp_sp=="15"])
H5_plast <-(gly_exp[,Pop_sp=="H15" & Temp_sp=="23"])-(gly_exp[,Pop_sp=="H15" & Temp_sp=="15"])

df.violin <- data.frame(plast=c(Base_plast[Base_plast<0], H1_plast[Base_plast<0], H2_plast[Base_plast<0], H3_plast[Base_plast<0], H4_plast[Base_plast<0], H5_plast[Base_plast<0]),Pop=rep(c("Ancestral","Hot-1","Hot-2","Hot-3","Hot-4","Hot-5"),each=sum(Base_plast<0)))

ggplot(df.violin, aes(x=Pop, y=plast,fill=Pop)) + 
  geom_violin(trim=TRUE)+scale_fill_manual(values=rep(c("chartreuse3","firebrick3"),c(1,5)))+geom_boxplot(width=0.1)

# A function to plot replicate variation of specific genes

do_plot<-function(gene_id,save_pdf=FALSE){

FBgn = as.character(gg$V1[gg$V4==gene_id])
name_pdf=paste(getwd(),"/",gene_id,".pdf",sep="")
if(save_pdf) pdf(name_pdf)
t=as.numeric(cpm(counts_BH,normalized.lib.sizes=TRUE)[expressed_genes==FBgn,])

plot(t ~ I(as.numeric(Temp_sp) + as.numeric(as.factor(substring(Pop_sp,1,1)))-1),pch=21,las=1,bty="n",xaxt="n",bg=rep(c("chartreuse3","firebrick3"),c(1,5))[as.numeric(as.factor(Pop_sp))],xlab="Temperature (°C)",ylab="cpm",main=gene_id,xlim=c(14,25))
axis(side=1,at=c(15.5,23.5),labels=c(15,23))

vect_BH <- as.numeric(as.factor(substring(Pop_sp,1,1)))
mean_cpm <- tapply(t,paste0(Temp_sp, vect_BH),mean)

points(c(15,23), mean_cpm[c(1,3)],type="l",lwd=2,col="chartreuse3")
points(c(16,24), mean_cpm[c(2,4)],type="l",lwd=2,col="firebrick3")

if(save_pdf) dev.off()

}

# All the glycolysis genes
gg=subset(genes,V4%in%glyc)
gg=gg[order(gg$V4),]
gg=gg[order(order(glyc)),]

pdf("Sup_file_Glycolysis_all.pdf")
par(mfrow=c(3,3))
for(i in gg$V4) do_plot(i, save_pdf=FALSE)
dev.off()

# All the ox pho genes
length(ox_pho_genes) #188
gg=subset(genes,V1%in%row.names(ox_pho))

pdf("Sup_file_Ox_pho_all.pdf")
par(mfrow=c(3,3))
for(i in gg$V4) do_plot(i, save_pdf=FALSE)
dev.off()

##############################
######### Figure S4 ##########
##############################

pdf("Fig_S4.pdf", width=20, height =20)
	layout.mat <- matrix(1:25,5,5)
	layout.mat <- t(layout.mat)
	layout(layout.mat, widths=rep(1,5) , heights=rep(1,5))
	par(mar=c(0,0,0,0))
	k=0
	for(i in 1:5){
	}
for(j in 1:5){
	for(i in 1:5){
	if(i==j){
		plot(0,0,type="n",xaxt="n",yaxt="n",ylab="",xlab="",bty="n")

	}else{	
	if(i>j){
	test_R1 <-  inc_plast[,j] & LRT_fits_DE[[11+j]]!=0
	test_R2 <-  inc_plast[,i] & LRT_fits_DE[[11+i]]!=0		
	vcol= "violetred1"		
	}else{
	test_R1 <-  dec_plast[,j] & LRT_fits_DE[[11+j]]!=0
	test_R2 <-  dec_plast[,i] & LRT_fits_DE[[11+i]]!=0		
	vcol="darkgoldenrod1"			
	}
	test_all <- test_R1 | test_R2

plot(LRT_fits[[6+i]]$table[test_all,1]~LRT_fits[[6+j]]$table[test_all,1],pch=16,asp=1,xlab=j,ylab=i,bty="n",xlim=c(-9,6),ylim=c(-9,6),bty="n",xaxt="n",yaxt="n",col= vcol)
print(k)
k=k+1
if(!k %in% c(1,5,9,13,17)){
	axis(side=2,at=c(-5,0,5),labels=c("","",""),pos=c(-6),cex=2,las=1)
	}else{
	axis(side=2,at=c(-5,0,5),pos=c(-6),cex=2,las=1)		
		}
if(!k %in% c(16:20)){
	axis(side=1,at=c(-5,0,5),labels=c("","",""),pos=c(-6),cex=2,las=1)	
}else{
	axis(side=1,at=c(-5,0,5),pos=c(-6),cex=2,las=1)	
	}

lines(c(-4,4),c(-4,4),lty=2)

}}
}	
	
dev.off()
