rm(list=ls())
gc()
library("edgeR")
library("gplots")
library("RColorBrewer")
library(ggplot2)
library(matrixStats)

#This code contains the code to produce the main analysis (all replicate together).

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

## Number of analyzed genes
# summary(keep)
#   Mode   FALSE    TRUE 
#logical    2062   11200

DGEList=DGEList(counts_merged[keep,])
DGEList=calcNormFactors(DGEList, method=c("TMM"))
Colors_temp=rep(c("chartreuse3","cornflowerblue","firebrick3"),each=2)

design_plasticity=model.matrix(~Pop_merged*Temp_merged)

####### A general model to see

DGEList=calcNormFactors(DGEList, method=c("TMM"))
BCH_DGEList_GLM <- estimateGLMRobustDisp(DGEList, design_plasticity,verbose=TRUE,maxit=10)
BCH_fit_DGEList_GLM=glmFit(BCH_DGEList_GLM, design=design_plasticity)

### Build the contrasts and save the results in a table

L1 = design_plasticity[c_B23_m,]-design_plasticity[c_B15_m,]
L2 = design_plasticity[c_C23_m,]-design_plasticity[c_C15_m,]
L3 = design_plasticity[c_H23_m,]-design_plasticity[c_H15_m,]

L4 = design_plasticity[c_C15_m,]-design_plasticity[c_B15_m,]
L5 = design_plasticity[c_H15_m,]-design_plasticity[c_B15_m,]

L6 = design_plasticity[c_C23_m,]-design_plasticity[c_B23_m,]
L7 = design_plasticity[c_H23_m,]-design_plasticity[c_B23_m,]

L8 = L3 - L1
L9 = L2 - L1

L <- cbind(L1,L2,L3,L4,L5,L6,L7,L8,L9)

# Fits LRT on each contrast
LRT_fits <- list()
LRT_fits_DE <- list()
LRT_fits_DE_low <- list()
for(k in 1: ncol(L)){
	LRT_fits[[k]] <- glmLRT(BCH_fit_DGEList_GLM, contrast=L[,k])	
} 

for(k in 1:ncol(L)) LRT_fits_DE[[k]] <- decideTestsDGE(LRT_fits[[k]], p=0.05, adjust="BH")
for(k in 1:ncol(L)) LRT_fits_DE_low[[k]] <- decideTestsDGE(LRT_fits[[k]], p=0.1, adjust="BH")

save(list=ls(),file="EdgeR_merged.RData")

#################################################################
###### Alternalively, load processed data from here	#############
#################################################################

load("EdgeR_merged.RData")

########################################
## Plastic genes in the base :
table(LRT_fits_DE[[1]])
# 2172 + 2180 = 4352
table(LRT_fits_DE[[2]]) #Cold
# 1896 + 1706 = 3602
table(LRT_fits_DE[[3]]) #Hot
# 2451 + 2458 = 4909

sum(LRT_fits_DE[[1]]!=0)/sum(keep) # 0.3885714
sum(LRT_fits_DE[[2]]!=0)/sum(keep) # 0.3216071
sum(LRT_fits_DE[[3]]!=0)/sum(keep) # 0.4383036

(sum(LRT_fits_DE[[3]]!=0)-sum(LRT_fits_DE[[1]]!=0))/sum(LRT_fits_DE[[1]]!=0) #0.1279871

expressed_genes <- row.names(counts_merged)[keep]
length(expressed_genes)
# Export table for GO analysis
write.table(expressed_genes,file="all_expressed.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(expressed_genes[LRT_fits_DE[[1]]==-1],file="Base_down.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(expressed_genes[LRT_fits_DE[[1]]==1],file="Base_up.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

Base_plast=cbind(data.frame(V1=c(expressed_genes[LRT_fits_DE[[1]]==-1], expressed_genes[LRT_fits_DE[[1]]==1])),rbind(LRT_fits[[1]]$table[LRT_fits_DE[[1]]==-1,c(1,4)],LRT_fits[[1]]$table[LRT_fits_DE[[1]]==1,c(1,4)]))
head(Base_plast)

# The FBgn symbols were translated in Flybase using the Tool "Upload/Convert IDs" on 14th April 2020 by uploading the tables "all_expressed.txt" generated above
trans=read.table("FBgn_translated.txt",h=F)
dim(trans)
head(trans)

B2=merge(Base_plast,trans[trans$V1%in%Base_plast$V1,c(1,4)])
# Three genes had multiple hits. They won't be translated
mult_hits <- subset(B2, V1%in% names(table(B2[,1]))[table(B2[,1])>1])
B2 <- subset(B2, !V1%in% names(table(B2[,1]))[table(B2[,1])>1])
mult_hits$V4 <- "unknown_ID"; mult_hits <- unique(mult_hits)
B2 <- rbind(B2, mult_hits)

dim(Base_plast); dim(B2) # OK
#Order
B2=B2[order(B2$logFC),]
head(B2)

write.table(B2[,c(1,4,2,3)],file="Base_plastic_translated.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

### GO comparison

GO_Chen=read.table("GO_class_Chen.txt")
head(GO_Chen)
GO_Port=read.table("GO_class_Port.txt")
head(GO_Port)

dim(subset(GO_Port,V2=="inc"))
dim(subset(GO_Chen,V2=="inc"))
table(GO_Chen[GO_Chen$V1%in%subset(GO_Port,V2=="inc")$V1,]$V2)
subset(GO_Chen[GO_Chen$V1%in%subset(GO_Port,V2=="inc")$V1,],V2=="inc")
sum(table(GO_Chen[GO_Chen$V1%in%subset(GO_Port,V2=="inc")$V1,]$V2)) # 107 present

dim(subset(GO_Port,V2=="dec"))
dim(subset(GO_Chen,V2=="dec"))
table(GO_Chen[GO_Chen$V1%in%subset(GO_Port,V2=="dec")$V1,]$V2)
sum(table(GO_Chen[GO_Chen$V1%in%subset(GO_Port,V2=="dec")$V1,]$V2)) # 107 present

names(GO_Port)=c("GO","Portugal")
names(GO_Chen)=c("GO","Chen")
dim(merge(GO_Chen,GO_Port)) #128

head(GO_Chen)
head(GO_Port)
GO_merged=merge(GO_Chen,GO_Port)
table(GO_merged$Port,GO_merged$Chen)

write.table(GO_merged,file="TableS2.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

### Correlation between ancestral plasticity and
### plasticity in the evolved populations

cor(LRT_fits[[1]]$table[,1],LRT_fits[[2]]$table[,1]) # 0.91
cor(LRT_fits[[1]]$table[,1],LRT_fits[[3]]$table[,1]) # 0.89

## Evolution of gene expression in the Cold evolved experiment

table(LRT_fits_DE[[4]]) # 18 + 24 -> 42 genes
table(LRT_fits_DE[[6]]) # 21 + 54 -> 75 genes

#table(LRT_fits_DE[[4]]!=0 & LRT_fits_DE[[6]]!=0) # 10 genes are differentially expressed in both environments
#data.frame(expressed_genes[(LRT_fits_DE_low[[4]]!=0 & LRT_fits_DE_low[[6]]!=0)])

# Export these genes for the GO analysis
# up : 95 genes
data.frame(expressed_genes[LRT_fits_DE_low[[6]]==1])
data.frame(expressed_genes[LRT_fits_DE_low[[6]]==-1])

###Change plasticity
data.frame(expressed_genes[LRT_fits_DE_low[[9]]!=0 & (LRT_fits_DE[[4]]!=0 | LRT_fits_DE[[6]]!=0) ])
data.frame(expressed_genes[LRT_fits_DE_low[[9]]!=0 & LRT_fits_DE[[4]]!=0 ])
data.frame(expressed_genes[LRT_fits_DE_low[[9]]!=0 & LRT_fits_DE[[6]]!=0 ])

###Export tables of evolved gene expression at 15°C / Cold
dim(data.frame(expressed_genes[LRT_fits_DE_low[[4]]%in%c(-1,1)])) # 58
BC_to_add=LRT_fits[[4]]$table[LRT_fits_DE_low[[4]]%in%c(-1,1),]
BC_to_add$V1=row.names(BC_to_add)
BC_trans_15=merge(trans, BC_to_add);dim(BC_trans_15) #58
write.table(BC_trans_15[,c(1,4:8)],file="BC15_Final.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


###Export tables of evolved gene expression at 23°C / Cold
dim(data.frame(expressed_genes[LRT_fits_DE_low[[6]]%in%c(-1,1)])) # 132
BC_to_add=LRT_fits[[6]]$table[LRT_fits_DE_low[[6]]%in%c(-1,1),]
BC_to_add$V1=row.names(BC_to_add)
BC_trans_23=merge(trans, BC_to_add);dim(BC_trans_23) #132
write.table(BC_trans_23[,c(1,4:8)],file="BC23_Final.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
 subset(BC_trans_23,V4=="sro")

## Evolution of gene expression in the Cold evolved experiment

table(LRT_fits_DE[[5]]) # 54 + 215 -> 269 genes
table(LRT_fits_DE[[7]]) # 298 + 427 -> 725 genes

## Ev. Plasticity
sum(LRT_fits_DE_low[[8]]!=0 & LRT_fits_DE[[5]]!=0) # 77
sum(LRT_fits_DE_low[[8]]!=0 & LRT_fits_DE[[7]]!=0) # 263

## Constitutive
sum(LRT_fits_DE[[5]]!=0 & LRT_fits_DE[[7]]!=0 & LRT_fits_DE[[7]]==LRT_fits_DE[[5]]) # 50

# Export
BH_const <- LRT_fits_DE[[5]]!=0 & LRT_fits_DE[[7]]!=0 & LRT_fits_DE[[7]]==LRT_fits_DE[[5]]
BH_const <- data.frame(V1=expressed_genes[BH_const])
BH_const =merge(trans, BH_const);dim(BH_const) # 50
write.table(BH_const[,c(1,4)],file="BH_const.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

#### BOTH Cold and Hot evolved at 23°C (ox. reduction) for Table S3

data.frame(expressed_genes[(LRT_fits_DE_low[[6]]==-1 & LRT_fits_DE_low[[7]]==-1)])
names_shared <- expressed_genes[(LRT_fits_DE_low[[6]]==-1 & LRT_fits_DE_low[[7]]==-1)]
length(names_shared)
B_shared <- merge(data.frame(V1=names_shared),trans)
write.table(B_shared,file="Shared_23_evolved.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

##############################
## Export data for Table S4 ##
##############################

dim(data.frame(expressed_genes[LRT_fits_DE[[5]]%in%c(-1,1)])) # 269
BH_to_add=LRT_fits[[5]]$table[LRT_fits_DE[[5]]%in%c(-1,1),]
BH_to_add$V1=row.names(BH_to_add);dim(BH_to_add)
BH_trans=merge(trans, BH_to_add);dim(BH_trans) # 269
write.table(BH_trans[,c(1,4:8)],file="BH15_Final.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


dim(data.frame(expressed_genes[LRT_fits_DE[[7]]%in%c(-1,1)])) # 725
BH_to_add=LRT_fits[[7]]$table[LRT_fits_DE[[7]]%in%c(-1,1),]
BH_to_add$V1=row.names(BH_to_add)
BH_trans=merge(trans, BH_to_add);dim(BH_trans) #725
write.table(BH_trans[,c(1,4:8)],file="BH23_Final.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

##############################
########## Figure 1 ##########
##############################

# Plot of hot evolved plasticity
v1=LRT_fits_DE[[3]]==0 & LRT_fits_DE[[1]]!=0
v2=LRT_fits_DE[[3]]!=0 & LRT_fits_DE[[1]]==0
v3=LRT_fits_DE[[3]]!=0 & LRT_fits_DE[[1]]!=0
color_code=rep("Not plastic",length(v1))
color_code[v1]="Ancestral only"
color_code[v2]="Hot ev. only"
color_code[v3]="Shared"

df_plot=as.data.frame(cbind(LRT_fits[[1]]$table[,1],LRT_fits[[2]]$table[,1],LRT_fits[[3]]$table[,1]))

names(df_plot)=c("Anc","Cold","Hot")

df_plot$color_code_Hot="Not plastic"
df_plot$color_code_Hot[LRT_fits_DE[[3]]==0 & LRT_fits_DE[[1]]!=0]="Ancestral only"
df_plot$color_code_Hot[LRT_fits_DE[[3]]!=0 & LRT_fits_DE[[1]]==0]="Hot ev. only"
df_plot$color_code_Hot[LRT_fits_DE[[3]]!=0 & LRT_fits_DE[[1]]!=0]="Shared"

df_plot$color_code_Cold="Not plastic"
df_plot$color_code_Cold[LRT_fits_DE[[2]]==0 & LRT_fits_DE[[1]]!=0]="Ancestral only"
df_plot$color_code_Cold[LRT_fits_DE[[2]]!=0 & LRT_fits_DE[[1]]==0]="Cold ev. only"
df_plot$color_code_Cold[LRT_fits_DE[[2]]!=0 & LRT_fits_DE[[1]]!=0]="Shared"

df_plot$color_code_Hot =factor(df_plot$color_code_Hot,levels=c("Hot ev. only","Cold ev. only","Ancestral only","Not plastic","Shared"))

df_plot$color_code_Cold=factor(df_plot$color_code_Cold,levels=c("Hot ev. only","Cold ev. only","Ancestral only","Not plastic","Shared"))

vec_color=c("firebrick3","cornflowerblue","forestgreen","grey","blueviolet")

par(mfrow=c(1,2))

plot(df_plot[,1],df_plot[,2],las=1,bty="n",type="n",asp=1,xlab="Ancestral plasticity",ylab="Cold evolved plasticity",xlim=c(-6.7,6.7))
abline(h=0,lty=2);abline(v=0,lty=2)
abline(a=0,b=1)
points(df_plot[,1],df_plot[,2],col=vec_color[as.numeric(df_plot$color_code_Cold)],pch=16,las=1,bty="n",cex=.7)
legend(1,-3,c("Hot evolved only","Cold evolved only","Ancestral only","Not plastic","Shared"),pch=16,col= vec_color,bty="n",cex=.8)

plot(df_plot[,1],df_plot[,3],las=1,bty="n",type="n",asp=1,xlab="Ancestral plasticity",ylab="Hot evolved plasticity",xlim=c(-6.7,6.7))
abline(h=0,lty=2);abline(v=0,lty=2)
abline(a=0,b=1)
points(df_plot[,1],df_plot[,3],col=vec_color[as.numeric(df_plot$color_code_Hot)],pch=16,las=1,bty="n",cex=.7)


#### End of Figure 1

###############################
### Evolution of plasticity ###
###############################

# Evolved plasticity in the Hot experiment
change_rn_BH=LRT_fits_DE_low[[8]]!=0 & (LRT_fits_DE[[7]]!=0 | LRT_fits_DE[[5]]!=0)
table(change_rn_BH) # 325 genes are plastic


#########

BH_trans=merge(trans, data.frame(V1=expressed_genes[change_rn_BH]));dim(BH_trans)

#Final export
#write.table(BH_trans[,c(1,4:8)],file="Change_rn_FlyBase_Final.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

evolved_Hot=(LRT_fits_DE[[7]]!=0 | LRT_fits_DE[[5]]!=0)
table(evolved_Hot)# 930

################################
# Genes that increase plasticity
#We compare the absolute values of the reaction norms

v15= LRT_fits_DE[[5]]!=0 & LRT_fits_DE_low[[8]]!=0 & abs(LRT_fits[[1]]$table[,1])<abs(LRT_fits[[3]]$table[,1]) 
v23= LRT_fits_DE[[7]]!=0 & LRT_fits_DE_low[[8]]!=0 & abs(LRT_fits[[1]]$table[,1])<abs(LRT_fits[[3]]$table[,1]) 
sum(v15)
sum(v23)
sum(v15 | v23)
sum(v15)+sum(v23)-sum(v15 | v23)  # 12 that evolved both

inc_plast=v15 | v23
sum(inc_plast) # 242

# Genes that decrease plasticity

v15= LRT_fits_DE[[5]]!=0 & LRT_fits_DE_low[[8]]!=0 & abs(LRT_fits[[1]]$table[,1])>abs(LRT_fits[[3]]$table[,1]) 
v23= LRT_fits_DE[[7]]!=0 & LRT_fits_DE_low[[8]]!=0 & abs(LRT_fits[[1]]$table[,1])>abs(LRT_fits[[3]]$table[,1])

sum(v15)
sum(v23)
sum(v15 | v23)
sum(v15)+sum(v23)-sum(v15 | v23)  # 3 that evolved both

dec_plast=v15 | v23
table(dec_plast) # 84
table(inc_plast) # 241 genes

# we can sample 84 genes multiple times in order to see if it depends on the number of genes ?
for(i in 1:100) write.table(file=paste0("sample_inc_plast_genes_",i,".txt"),expressed_genes[inc_plast][sample(1:242,85)],row.names=FALSE,quote=FALSE,col.names=FALSE)

data.frame(expressed_genes[dec_plast])

write.table(merge(trans,cbind(data.frame(V1=expressed_genes[dec_plast]),LRT_fits[[8]]$table[dec_plast,]))[,c(1,4:8)],file="BH_dec_plast.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(merge(trans,cbind(data.frame(V1=expressed_genes[inc_plast]),LRT_fits[[8]]$table[inc_plast,]))[,c(1,4:8)],file="BH_inc_plast.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

# Restricted to ancestrally plastic genes only
table(dec_plast & LRT_fits_DE[[1]]!=0 & abs(LRT_fits[[1]]$table[,1]) > 1) # 20
table(inc_plast & LRT_fits_DE[[1]]!=0 & abs(LRT_fits[[1]]$table[,1]) > 1) # 62 genes

################################

summary(glm(rep(c(0,1),c(84,241))~1,family="binomial"))
summary(glm(rep(c(0,1),c(21,63))~1,family="binomial"))

################################
## Is there an overlap with the genes showing a constitutive change ?

table(inc_plast & expressed_genes%in%BH_const$V1)
table(dec_plast & expressed_genes%in%BH_const$V1)
# Only 1 gene
counts_cpm=cpm(counts_merged,normalized.lib.sizes=TRUE)[keep,]
vect_temp= counts_cpm[expressed_genes[inc_plast & expressed_genes%in%BH_const$V1],]


#######


## EXPORT FOR TABLE
#BH_plast_dn=read.table("BH_plast_down_translated.txt",sep="\t")[,c(1,2,4)]


#BH_trans_const=merge(BH_trans23,BH_trans15)
#dim(BH_trans_const) # 151 - 117 + 34 in diff direction
#BH_trans_const= BH_trans_const[BH_trans_const$V1%in%row.names(LRT_fits_DE[[7]])[constitutive_BH],]
#dim(BH_trans_const) # 117 OK
#write.table(BH_trans_const,file="Constitutive_FlyBase_Final.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

## Overlap between evolution at 15°C and 23°C 

table(LRT_fits_DE[[7]]!=0 & LRT_fits_DE[[5]]!=0) # 64 genes

####### Genes with increased plasticity
v15= LRT_fits_DE[[5]]!=0 & LRT_fits_DE_low[[8]]!=0 & abs(LRT_fits[[1]]$table[,1])<abs(LRT_fits[[3]]$table[,1]) 
v23= LRT_fits_DE[[7]]!=0 & LRT_fits_DE_low[[8]]!=0 & abs(LRT_fits[[1]]$table[,1])<abs(LRT_fits[[3]]$table[,1]) 

v15l= LRT_fits_DE_low[[5]]!=0 & LRT_fits_DE_low[[8]]!=0 & abs(LRT_fits[[1]]$table[,1])<abs(LRT_fits[[3]]$table[,1]) 

v23l= LRT_fits_DE_low[[7]]!=0 & LRT_fits_DE_low[[8]]!=0 & abs(LRT_fits[[1]]$table[,1])<abs(LRT_fits[[3]]$table[,1]) 

#################
### Figure 2A ###
#################

par(mfrow=c(1,2))

plot(LRT_fits[[5]]$table[,1],LRT_fits[[7]]$table[,1],type="n",xlab="Change in expression at 15°C (log2FC)",ylab="Change in expression at 23°C (log2FC)",las=1,bty="n",xlim=c(-1,2),ylim=c(-4,4))
abline(h=0);abline(v=0)

points(LRT_fits[[5]]$table[v15 & !v23,1],LRT_fits[[7]]$table[v15 & !v23,1],col="cornflowerblue",pch=16)

points(LRT_fits[[5]]$table[!v15 & v23,1],LRT_fits[[7]]$table[!v15 & v23,1],col="firebrick3",pch=16)

points(LRT_fits[[5]]$table[v15 & v23,1],LRT_fits[[7]]$table[v15 & v23,1],col="plum",pch=16)

legend(.5,4,c("DE at 15°C","DE at both T°","DE at 23°C"),col=c("cornflowerblue","plum","firebrick3"),pch=16,bty="n")

#chi-square test for the correlation
chisq.test(sign(LRT_fits[[5]]$table[v15 | v23,1]),sign(LRT_fits[[7]]$table[v15 | v23,1]))
chisq.test(sign(LRT_fits[[5]]$table[v15l | v23l,1]),sign(LRT_fits[[7]]$table[v15l | v23l,1]))

chisq.test(sign(LRT_fits[[5]]$table[inc_plast,1]),sign(LRT_fits[[7]]$table[inc_plast,1]))
chisq.test(sign(LRT_fits[[5]]$table[dec_plast,1]),sign(LRT_fits[[7]]$table[dec_plast,1]))

sum(v15 | v23)
sum(inc_plast)


