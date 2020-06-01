rm(list=ls())
library(data.table)
library(edgeR)

## This file load the count table containing all the samples, produce the MDS plot
## and merge the libraries, then export the final table.


# Read the count table
counts <- read.table(file="Counts_all.txt",sep='\t',h=TRUE)


# Extract the sample ID from the count header
extract_names <- tstrsplit(names(counts),'_')
Pop <- extract_names[[1]]
Temp <- extract_names[[2]]
Repl <- extract_names[[3]]
Tech <- extract_names[[4]]


######## The MDS plot in supplement ########
Colors_temp=c("chartreuse3","cornflowerblue","firebrick3")
Repl2 <- rep(1:5,3)[as.numeric(as.factor(Repl))]

keep=log2(rowSums(cpm(counts,normalized.lib.sizes=TRUE))/dim(counts)[2])>0
summary(keep)

Pop_MDS <- Pop
Pop_MDS[Pop_MDS=="B"] <- "A"

quartz()
par(mfrow=c(1,2))
for(i in c(15,23)){
DGEList=DGEList(counts[keep,Temp==i])
DGEList=calcNormFactors(DGEList, method=c("TMM"))

plotMDS.DGEList(DGEList,method="logFC", xlab="LogFC dim1", ylab="LogFC dim2", col=Colors_temp[as.numeric(as.factor(substring(names(counts)[Temp==i],1,3)))],labels=paste(Pop_MDS, Repl2,sep="")[Temp==i],bty="n",las=1,xlim=c(-1.1,1.1),asp=1,main=paste0(i,"°C Common Garden"))
}
legend(-1.25,1.1,c("Ancestral","Hot evolved","Cold evolved"),pch=16,col=Colors_temp[c(1,3,2)],bty="n")

counts_merged <- NULL
kkk=0
vect_names <- NULL
for(i in unique(Pop)){
	for(j in unique(Temp)){
		for(k in unique(Repl[Pop==i & Temp==j])){
			kkk=kkk+1
			vect_names<- c(vect_names,paste(i,j,k,sep="_"))
			if(sum(Pop==i & Temp==j & Repl==k)>1){
			counts_merged <- cbind(counts_merged ,rowSums(counts[,Pop==i & Temp==j & Repl==k]))
			}else{
			counts_merged <- cbind(counts_merged ,counts[,Pop==i & Temp==j & Repl==k])
			}
		}
	}
}
counts_merged<- as.data.frame(counts_merged)
names(counts_merged) <- vect_names


write.table(counts_merged,file="Counts_merged.txt",sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
