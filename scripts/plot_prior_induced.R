##############################################################
#
# Script to check the induced priors (with ESS) on fixed effects
#
#
##############################################################

library(MASS)

invlogit <- function(x) 1/(1+exp(-x))


I = 7  #number of doses

Ntot = 10000  # number of samples

delta_doses = c(2,3)  # to be used for jump less than 2: plot the first, 2 jumps and 5 jumps

##################################################################

mu_norm = -2
sigma_norm = 5
var_coeff= 0.5
slope_gamma = 0.667 

target_MTD=0.33

set.seed(125896)


shape_gamma = 1/var_coeff^2
rate_gamma = (var_coeff^2*slope_gamma)^-1


######## first dose
res = NULL
for (n in 1:Ntot){
  res_sample1 = c(rnorm(1,mu_norm, sigma_norm), rgamma(length(delta_doses), 
                                                       shape=shape_gamma*delta_doses, rate=rate_gamma))
  res_sample1 = cumsum(res_sample1)
  res= rbind(res, invlogit(res_sample1))
}

round(apply(res, 2, median),3)

res=data.frame(res)
dimnames(res)[[2]]<- c("first_dose", "two_jumps", "five_jumps")


x=seq(0.001,0.999,0.001)
library(ggplot2)
ab=fitdistr(res[,1],"beta",list(shape1=1,shape2=1)) 
plotd = data.frame(cbind(x=seq(0.001,0.999,0.001), y=dbeta(x, ab$estimate[1], ab$estimate[2])))
pprior1 =  ggplot(res, aes(x=first_dose)) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.title = element_blank(), legend.text = element_text(size = 12, face = "bold"),
                     plot.title = element_text(size=20, face="bold", hjust = 0.5),
                     axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = element_text(size=16, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = element_text(face = "bold", size = 14),
                     axis.text.y = element_text(face = "bold", size = 14),
                     legend.position =  "none",
                     legend.key.size = unit(1.5, 'lines'),
                     strip.text = element_text(size=15)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  ggtitle("First dose") +  
  xlab("probability of toxicity") +
  geom_line(data=plotd, aes(x=x, y=y), size=0.5, linetype = "dashed") +
  ylim(0,10) +
  geom_text(x=0.5, y=7, label=paste("ESS = ", round(sum(ab$estimate),3),sep=""), aes(fontface=2), size=6) +
  geom_text(x=0.5, y=6, label=paste("prior median = ", round(median(res[,1]),3),sep=""), aes(fontface=2), size=6)

#pprior1 

ab=fitdistr(res[,2],"beta",list(shape1=1,shape2=1)) 
plotd = data.frame(cbind(x=seq(0.001,0.999,0.001), y=dbeta(x, ab$estimate[1], ab$estimate[2])))
pprior2 =  ggplot(res, aes(x=two_jumps)) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.title = element_blank(), legend.text = element_text(size = 12, face = "bold"),
                     plot.title = element_text(size=20, face="bold", hjust = 0.5),
                     axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = element_blank(),
                     axis.text.x = element_text(face = "bold", size = 14),
                     axis.text.y = element_text(face = "bold", size = 14),
                     legend.position =  "none",
                     legend.key.size = unit(1.5, 'lines'),
                     strip.text = element_text(size=15)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  ggtitle("Two units of increment") +  
  xlab("probability of toxicity") +
  geom_line(data=plotd, aes(x=x, y=y), size=0.5, linetype = "dashed") +
  ylim(0,10) +
  geom_text(x=0.5, y=7, label=paste("ESS = ", round(sum(ab$estimate),3),sep=""), aes(fontface=2), size=6) +
  geom_text(x=0.5, y=6, label=paste("prior median = ", round(median(res[,2]),3),sep=""), aes(fontface=2), size=6)

#pprior2 


ab=fitdistr(res[,3],"beta",list(shape1=1,shape2=1)) 
plotd = data.frame(cbind(x=seq(0.001,0.999,0.001), y=dbeta(x, ab$estimate[1], ab$estimate[2])))
pprior3 =  ggplot(res, aes(x=five_jumps)) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.title = element_blank(), legend.text = element_text(size = 12, face = "bold"),
                     plot.title = element_text(size=20, face="bold", hjust = 0.5),
                     axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = element_blank(),
                     axis.text.x = element_text(face = "bold", size = 14),
                     axis.text.y = element_text(face = "bold", size = 14),
                     legend.position =  "none",
                     legend.key.size = unit(1.5, 'lines'),
                     strip.text = element_text(size=15)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  ggtitle("Five units of increment") +  
  xlab("probability of toxicity") +
  geom_line(data=plotd, aes(x=x, y=y), size=0.5, linetype = "dashed") +
  ylim(0,10) +
  geom_text(x=0.5, y=7, label=paste("ESS = ", round(sum(ab$estimate),3),sep=""), aes(fontface=2), size=6) +
  geom_text(x=0.5, y=6, label=paste("prior median = ", round(median(res[,3]),3),sep=""), aes(fontface=2), size=6)

#pprior3 

library(gridExtra)


jpeg("p_less_2_only3.jpeg", width = 20, height = 10, unit="in", res=300)
grid.arrange(pprior1, pprior2, pprior3, nrow = 1)
dev.off()



############################################# jump > 2
Ntot = 10000  # number of samples

delta_doses = c(5,4)

##################################################################

mu_norm = -4
sigma_norm = 3.5
var_coeff= 0.5
slope_gamma = 0.642 

target_MTD=0.33

set.seed(125896)


shape_gamma = 1/var_coeff^2
rate_gamma = (var_coeff^2*slope_gamma)^-1


######## first dose
res = NULL
for (n in 1:Ntot){
  res_sample1 = c(rnorm(1,mu_norm, sigma_norm), rgamma(length(delta_doses), 
                                                       shape=shape_gamma*delta_doses, rate=rate_gamma))
  res_sample1 = cumsum(res_sample1)
  res= rbind(res, invlogit(res_sample1))
}

round(apply(res, 2, median),3)

res=data.frame(res)
dimnames(res)[[2]]<- c("first_dose", "five_jumps", "nine_jumps")



library(ggplot2)
ab=fitdistr(res[,1],"beta",list(shape1=1,shape2=1)) 
plotd = data.frame(cbind(x=seq(0.001,0.999,0.001), y=dbeta(x, ab$estimate[1], ab$estimate[2])))
pprior1 =  ggplot(res, aes(x=first_dose)) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.title = element_blank(), legend.text = element_text(size = 12, face = "bold"),
                     plot.title = element_text(size=20, face="bold", hjust = 0.5),
                     axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = element_text(size=16, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = element_text(face = "bold", size = 14),
                     axis.text.y = element_text(face = "bold", size = 14),
                     legend.position =  "none",
                     legend.key.size = unit(1.5, 'lines'),
                     strip.text = element_text(size=15)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  ggtitle("First dose") +  
  xlab("probability of toxicity") +
  geom_line(data=plotd, aes(x=x, y=y), size=0.5, linetype = "dashed") +
  ylim(0,15) +
  geom_text(x=0.5, y=10, label=paste("ESS = ", round(sum(ab$estimate),3),sep=""), aes(fontface=2), size=6) +
  geom_text(x=0.5, y=9, label=paste("prior median = ", round(median(res[,1]),3),sep=""), aes(fontface=2), size=6)

#pprior1 

ab=fitdistr(res[,2],"beta",list(shape1=1,shape2=1)) 
plotd = data.frame(cbind(x=seq(0.001,0.999,0.001), y=dbeta(x, ab$estimate[1], ab$estimate[2])))
pprior2 =  ggplot(res, aes(x=five_jumps)) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.title = element_blank(), legend.text = element_text(size = 12, face = "bold"),
                     plot.title = element_text(size=20, face="bold", hjust = 0.5),
                     axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = element_blank(),
                     axis.text.x = element_text(face = "bold", size = 14),
                     axis.text.y = element_text(face = "bold", size = 14),
                     legend.position =  "none",
                     legend.key.size = unit(1.5, 'lines'),
                     strip.text = element_text(size=15)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  ggtitle("Five units of increment") +  
  xlab("probability of toxicity") +
  geom_line(data=plotd, aes(x=x, y=y), size=0.5, linetype = "dashed") +
  ylim(0,15) +
  geom_text(x=0.5, y=10, label=paste("ESS = ", round(sum(ab$estimate),3),sep=""), aes(fontface=2), size=6) +
  geom_text(x=0.5, y=9, label=paste("prior median = ", round(median(res[,2]),3),sep=""), aes(fontface=2), size=6)

#pprior2 


ab=fitdistr(res[,3],"beta",list(shape1=1,shape2=1)) 
plotd = data.frame(cbind(x=seq(0.001,0.999,0.001), y=dbeta(x, ab$estimate[1], ab$estimate[2])))
pprior3 =  ggplot(res, aes(x=nine_jumps)) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.title = element_blank(), legend.text = element_text(size = 12, face = "bold"),
                     plot.title = element_text(size=20, face="bold", hjust = 0.5),
                     axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = element_blank(),
                     axis.text.x = element_text(face = "bold", size = 14),
                     axis.text.y = element_text(face = "bold", size = 14),
                     legend.position =  "none",
                     legend.key.size = unit(1.5, 'lines'),
                     strip.text = element_text(size=15)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  ggtitle("Nine units of increment") +  
  xlab("probability of toxicity") +
  geom_line(data=plotd, aes(x=x, y=y), size=0.5, linetype = "dashed") +
  ylim(0,15) +
  geom_text(x=0.5, y=10, label=paste("ESS = ", round(sum(ab$estimate),3),sep=""), aes(fontface=2), size=6) +
  geom_text(x=0.5, y=9, label=paste("prior median = ", round(median(res[,3]),3),sep=""), aes(fontface=2), size=6)

#pprior3 

library(gridExtra)

jpeg("p_greater_2_only3.jpeg", width = 20, height = 10, unit="in", res=300)
grid.arrange(pprior1, pprior2, pprior3, nrow = 1)
dev.off()

