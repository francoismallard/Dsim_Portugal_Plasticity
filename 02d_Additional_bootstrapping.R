

# This code produces the analysis shown in Figure S3 - It follows "02_Analyse_gene_expression.R"


plast_low <- LRT_fits_DE_low[[8]]!=0
sum(inc_plast) # 241
sum(dec_plast) # 84
sum(plast_low) # 417

BH15_plast <- LRT_fits[[5]]$table[plast_low,1]
BH23_plast <- LRT_fits[[7]]$table[plast_low,1]

temp_stat <- NULL
temp_cor <- NULL
temp_stat2 <- NULL
for(i in 1:10000){
vect_sample <- sample(c(1:length(BH15_plast)),sum(inc_plast),replace=FALSE)

temp_stat <-c(temp_stat, chisq.test(
sign(BH15_plast[vect_sample]),
sign(BH23_plast[vect_sample])
)$statistic)


temp_stat2 <-c(temp_stat2, chisq.test(
(BH15_plast[vect_sample]),
(BH23_plast[vect_sample])
)$statistic)
temp_cor <- c(temp_cor ,cor(BH15_plast[vect_sample],BH23_plast[vect_sample]))
}

true_stat <- chisq.test(
sign(LRT_fits[[5]]$table[inc_plast,1]),
sign(LRT_fits[[7]]$table[inc_plast,1]),
)$statistic

hist(temp_stat,n=100,main="",xlab=expression(paste(chi^2," statistics")))
abline(v=(true_stat),col="red")
abline(v=sort(temp_stat)[9500],col="blue")
sum(temp_stat> true_stat) # 75







