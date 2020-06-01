
# This code produces Figure S2- It follows "02_Analyse_gene_expression.R"

genes=read.table("genesID2.txt",sep="\t")

ox_pho_genes=unique(read.table("GO_Oxidative_phosphorylation.txt",sep="\t")$V1)
length(ox_pho_genes) # 188

ox_pho= counts_cpm[expressed_genes%in% ox_pho_genes,]
for(i in 1:dim(ox_pho)[1]) ox_pho[i,]=ox_pho[i,]/max(ox_pho[i,])
dim(ox_pho) # 157 genes

ox_pho_DE=ox_pho[LRT_fits_DE[[7]][expressed_genes%in%row.names(ox_pho)]!=0,]
dim(ox_pho_DE) # 29

length(ox_pho_genes) #188
gg_sav=subset(genes,V1%in%row.names(ox_pho))
dim(gg_sav) # 146

############

# Too many genes, the graph is split in three

gg=gg_sav[1:48,]
gg=gg_sav[49:97,]
gg=gg_sav[98:146,]


### Produce the barplot on the "gg" subset

DE_t1= LRT_fits[[7]]$table[expressed_genes%in%gg$V1,]
DE_t1=DE_t1[order(row.names(DE_t1)),]
DE_t1=DE_t1[order(order(gg$V1)),]

DE_t2= LRT_fits[[5]]$table[expressed_genes%in%gg$V1,]
DE_t2=DE_t2[order(row.names(DE_t2)),]
DE_t2=DE_t2[order(order(gg$V1)),]

LRT_fits_DE_low[[7]]= LRT_fits_DE_low[[7]]
LRT_fits_DE_low[[5]]=LRT_fits_DE_low[[5]]


######################################################################################################

#quartz(height=5.3,width=8)
barplot(rbind(DE_t1,DE_t2)[rep(c(1:length(gg$V1)),each=2)+rep(c(0,38),length(gg$V1)),]$logFC,col=c("firebrick3","cornflowerblue"),ylab="logFC",main="",density=c(20,100),ylim=c(-1.05,.5),width=I(2/2.5),space=c(0,rep(c(0,0.5),(nrow(gg)-1)),0),xlim=c(1,(2*nrow(gg))))

lines(rep(0,2),c(0.1,-.85),lty=2)
for(i in 1:nrow(gg)) lines(rep(I(2*i-0.2),2),c(0.1,-.85),lty=2)
text(I(c(1:nrow(gg))*2-1),-0.8,labels=gg$V4,srt=90)

for(i in c(1:nrow(gg))){
	temp=LRT_fits_DE[[7]][expressed_genes==gg[i,]$V1]
	temp2=LRT_fits_DE_low[[7]][expressed_genes==gg[i,]$V1]
	if(temp==-1){
		points(I(2*i+1/2.5)-2,0.1,pch=8,col="firebrick3")
		points(I(2*i+1/2.5)-2,0.15,pch=8,col="firebrick3")
		}else if(temp2==-1){
		points(I(2*i+1/2.5)-2,0.1,pch=8,col="firebrick3")			
		}
	if(temp==1){
		 points(I(2*i+1/2.5)-2,-0.1,pch=8,col="firebrick3")
		 points(I(2*i+1/2.5)-2,-0.15,pch=8,col="firebrick3")
		 }else if(temp2==1){
		 points(I(2*i+1/2.5)-2,-0.1,pch=8,col="firebrick3")		 	
		 }

	temp=LRT_fits_DE[[5]][expressed_genes==gg[i,]$V1]
	temp2=LRT_fits_DE_low[[5]][expressed_genes==gg[i,]$V1]
	if(temp==-1){
		points(I(2*i+1/2.5)-2+2/2.5,0.15,pch=8,col="cornflowerblue")
		points(I(2*i+1/2.5)-2+2/2.5,0.1,pch=8,col="cornflowerblue")
		}else if(temp2==-1){
			points(I(2*i+1/2.5)-2+2/2.5,0.1,pch=8,col="cornflowerblue")
			}
	if(temp==1){
		points(I(2*i+1/2.5)-2+2/2.5,-0.1,pch=8,col="cornflowerblue")
		points(I(2*i+1/2.5)-2+2/2.5,-0.15,pch=8,col="cornflowerblue")
		}else if(temp2==1){
		points(I(2*i+1/2.5)-2+2/2.5,-0.1,pch=8,col="cornflowerblue")	
		}

}
