
# This code produces Figure 3A- It follows "02_Analyse_gene_expression.R"

genes=read.table("genesID2.txt",sep="\t")

glyc=c("Hex-A","Hex-C","Pgi","Pfk","fbp","Tpi","Ald","Gapdh1","Gapdh2","Pgk","Pglym78","Eno","PyK","ImpL3","CG11876","CG5261","CG7430","l(1)G0334")

gg=subset(genes,V4%in%glyc)

gg=gg[order(gg$V4),]
gg=gg[order(order(glyc)),]

summary(gg$V4==glyc)

DE_t1= LRT_fits[[7]]$table[expressed_genes%in%gg$V1,]
DE_t1=DE_t1[order(row.names(DE_t1)),]
DE_t1=DE_t1[order(order(gg$V1)),]

DE_t2= LRT_fits[[5]]$table[expressed_genes%in%gg$V1,]
DE_t2=DE_t2[order(row.names(DE_t2)),]
DE_t2=DE_t2[order(order(gg$V1)),]

BCH_LRT_BH_23_DE_lowFDR= LRT_fits_DE_low[[7]]
BCH_LRT_BH_15_DE_lowFDR=LRT_fits_DE_low[[5]]

quartz(height=5.3,width=8)
barplot(rbind(DE_t1,DE_t2)[rep(c(1:length(gg$V1)),each=2)+rep(c(0,18),length(gg$V1)),]$logFC,col=c("firebrick3","cornflowerblue"),ylab="logFC",main="",density=c(20,100),ylim=c(-.95,.5),width=I(2/2.5),space=c(0,rep(c(0,0.5),17)),xlim=c(1,36))

lines(rep(0,2),c(0.1,-.85),lty=2)
for(i in 1:18) lines(rep(I(2*i-0.2),2),c(0.1,-.85),lty=2)
text(I(c(1:18)*2-1),-0.8,labels=gg$V4,srt=90)

for(i in c(1:18)){
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
	temp2= LRT_fits_DE_low[[5]][expressed_genes==gg[i,]$V1]
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

