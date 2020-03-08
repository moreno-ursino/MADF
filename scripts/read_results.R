##############################################################
#
# Script to build result table for the selected scenario
#
# ssss      the interger number associated to the scenario (called scen in other script)
# nnnn      NTcrm + NT3p3
##############################################################

ssss=100
nnnn=10


MTD = function(x, theta){
  
  order(abs(x-theta))[1]
}

theta_MTD = 0.33

eval(parse(text = paste("mMADF=read.table(\"res_median_MADFscen",ssss,"_",nnnn,".txt\")", sep="")))
eval(parse(text = paste("mMADF1=read.table(\"res_median_MADF1scen",ssss,"_",nnnn,".txt\")", sep="")))
eval(parse(text = paste("mMADF2=read.table(\"res_median_MADF2scen",ssss,"_",nnnn,".txt\")", sep="")))
eval(parse(text = paste("mMADF3=read.table(\"res_median_MADF3scen",ssss,"_",nnnn,".txt\")", sep="")))
eval(parse(text = paste("mMADF4=read.table(\"res_median_MADF4scen",ssss,"_",nnnn,".txt\")", sep="")))
eval(parse(text = paste("mZKO=read.table(\"res_ZKO_ZKO",ssss,"_",nnnn,".txt\")", sep="")))


# MTD selection
MTD_mMADF <- apply(mMADF, 1, MTD, theta=theta_MTD)
MTD_mMADF1 <- apply(mMADF1, 1, MTD, theta=theta_MTD)
MTD_mMADF2 <- apply(mMADF2, 1, MTD, theta=theta_MTD)
MTD_mMADF3 <- apply(mMADF3, 1, MTD, theta=theta_MTD)
MTD_mMADF4 <- apply(mMADF4, 1, MTD, theta=theta_MTD)
MTD_mZKO <- apply(mZKO, 1, MTD, theta=theta_MTD)

# proportion
xMADF = prop.table(table(MTD_mMADF))
xMADF1 = prop.table(table(MTD_mMADF1))
xMADF2 = prop.table(table(MTD_mMADF2))
xMADF3 = prop.table(table(MTD_mMADF3))
xMADF4 = prop.table(table(MTD_mMADF4))
xZKO = prop.table(table(MTD_mZKO))


res_tot_s <- matrix(0, ncol=7, nrow=6)
res_tot_s[1, as.numeric(as.character(dimnames(xMADF)$MTD_mMADF))] = xMADF
res_tot_s[2, as.numeric(as.character(dimnames(xMADF1)$MTD_mMADF1))] = xMADF1
res_tot_s[3, as.numeric(as.character(dimnames(xMADF2)$MTD_mMADF2))] = xMADF2
res_tot_s[4, as.numeric(as.character(dimnames(xMADF3)$MTD_mMADF3))] = xMADF3
res_tot_s[5, as.numeric(as.character(dimnames(xMADF4)$MTD_mMADF4))] = xMADF4
res_tot_s[6, as.numeric(as.character(dimnames(xZKO)$MTD_mZKO))] = xZKO


# to add the number of patients at each dose
eval(parse(text = paste("pat_res=read.table(\"patient_scen",ssss,"_",nnnn,".txt\")", sep="")))

mat_toplot=apply(pat_res,2, quantile, probs=c(0.5, 0.25, 0.75))
toplot = NULL
for(i in 1:dim(mat_toplot)[2]) toplot[i] = paste(mat_toplot[1,i], " (", mat_toplot[2,i], ", ", mat_toplot[3,i], ")", sep="")

library(xtable)
final_table = xtable(rbind(res_tot_s,toplot), digits = 3, caption = "Results", align="rccccccc")

rownames( final_table ) <- c( "MADF", "MADF1", "MADF2", "MADF3", "MADF4", "ZKO", "#patients")
final_table


