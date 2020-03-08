##############################################################
#
# Script to run MADF on Sorafenib and Irinotecan + S-1 examples
#
#
##############################################################

require("rstan")
require(Iso)
install.packages("ggrepel")
library(ggrepel)

# non centered parametrization, it is good
sm_meta <- stan_model(file="MADF.stan", model_name='sm-meta1',verbose=FALSE)

priors_matrix = data.frame(rbind(
  c(-2.0,5.0,0.667,0.5),
  c(-4.0,3.5,0.642,0.5)
))

colnames(priors_matrix) <-c("mu_norm", "sigma_norm", "slope_gamma", "var_gamma")

SED4=1478963

# for stan
nchains=4


#### data reading
data <- read.csv("sorafenib2.csv")

# derive data characteristics:
N         <- nrow(data)              # number of data points (study/dose combinations)
studyIdx  <- as.numeric(data$study)  # study indices
studyName <- levels(data$study)      # study names
K         <- max(studyIdx)           # number of studies

# show some study-level summary statistics:
sumstats <- cbind("doses"=table(data$study),
                  "patients"=tapply(data$total, data$study, sum),
                  "events"=tapply(data$events, data$study, sum))
print(sumstats)


doses_true <<- c(100, 200, 300, 400, 600, 800, 1000)    #in mg for example
doses_modif <<- doses_true/mean(doses_true)
unit_dose= doses_modif[2]-doses_modif[1] 
# specify data provided to Stan

delta_dose <- (doses_modif[2:length(doses_modif)] - doses_modif[1:(length(doses_modif)-1)]) / unit_dose


nDLT_meta <- function(doses, data) {
  sum(data[data$dose == doses, "events"])
}
nPT_meta <- function(doses, data) {
  sum(data[data$dose == doses, "total"])
}
ti <- round(sapply(doses_true, nDLT_meta, data = data), 0)
ni <- round(sapply(doses_true, nPT_meta, data = data), 0)

EB_prior= pava(ti/ni)
pos_EB = order(abs(EB_prior-0.33))[1]
pos_prior = ( ((doses_modif[pos_EB] - doses_modif[1])/unit_dose) >=2 )*1


  doselevel_tot = doses_true
  I = length(doselevel_tot)
  doselevel=match(data$dose, doselevel_tot)
  events =  data$events
  total = data$total
  
  var_coeff= priors_matrix[pos_prior+1, "var_gamma"]
  slope_gamma = priors_matrix[pos_prior+1, "slope_gamma"]
  mu_first = priors_matrix[pos_prior+1, "mu_norm"]
  sd_first = priors_matrix[pos_prior+1, "sigma_norm"]
  
  shape_gamma = 1/var_coeff^2
  rate_gamma = (var_coeff^2*slope_gamma)^-1
  
  stanData <- list(K        = K,
                   N       = N,
                   I       = I,
                   studyIdx = studyIdx,
                   doselevel     = doselevel,
                   events  = events,
                   total   = total,
                   doses_modif = doses_modif,
                   shape_gamma  = shape_gamma,
                   rate_gamma  = rate_gamma,
                   delta_dose = delta_dose,
                   mu_first = mu_first,
                   sd_first = sd_first)
  
  set.seed(SED4)
  
  fit_mod <- sampling(sm_meta, data=stanData ,iter=4000, chains=nchains, control = list(adapt_delta = 0.95))
  results <- extract(fit_mod, pars="ptox")$ptox
  

  
  
  ptox_res <- apply(results,2,median)
  ptox_res2 <- apply(results,2,mean)
  q1_res<- apply(results,2,quantile, probs=0.25)
  q3_res <- apply(results,2,quantile, probs=0.75)
  


results_Sorafenib = results

library(ggplot2)
str(results_Sorafenib)
doses_true <<- c(100, 200, 300, 400, 600, 800, 1000)  
forplot <- NULL
for (i in 1:dim(results_Sorafenib)[2]) forplot <- rbind(forplot, cbind(results_Sorafenib[,i], 
                                                rep(doses_true[i],dim(results_Sorafenib)[1] )))
dimnames(forplot)[[2]] <- c("Probability", "Dose")
forplot <- data.frame(forplot)
#forplot$Dose <- as.factor(forplot$Dose)

p <- ggplot(forplot, aes(x=Dose, y=Probability, group=factor(Dose))) +  
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     plot.title = element_text(size=20, face="bold", hjust = 0.5),
                     axis.title.x = element_text(size=14, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = element_text(size=14, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = element_text(face = "bold", size = 12),
                     axis.text.y = element_text(face = "bold", size = 12),
                     legend.position = "none") +
  ggtitle("Sorafenib") +
  xlab("dose (mg)") + ylab("Probability of toxicity") +
  geom_violin(trim=FALSE, size=1.1, fill="#D3D3D3") +
  geom_hline(yintercept=0.33, color = "black", size=0.8, linetype="dashed") +
  geom_hline(yintercept=0.25, color = "black", size=0.8, linetype="dashed") +
  geom_hline(yintercept=0.20, color = "black", size=0.8, linetype="dashed") +
  coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks=c(0,0.20,0.25,0.33, 0.5,0.75,1)) +
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  scale_x_continuous(breaks = doses_true )
p

jpeg("plot_Sorafenib.jpeg", width = 8, height = 5, unit="in", res=300)
p 
dev.off()


#######  EWOC part

results_Ewok = results_Sorafenib

checking <- function(x, tau) mean(x>=tau)

round(apply(results_Ewok, 2, checking, tau=0.33),3)
# 0.000 0.000 0.000 0.000 0.369 0.991 1.000
round(apply(results_Ewok, 2, checking, tau=0.25),3)
# 0.000 0.000 0.000 0.002 0.832 1.000 1.000
round(apply(results_Ewok, 2, checking, tau=0.20),3)
# 0.000 0.000 0.001 0.016 0.964 1.000 1.000

##################### Irinotecan

priors_matrix = data.frame(rbind(
  c(-2.0,5.0,0.667,0.5),
  c(-4.0,3.5,0.642,0.5)
))

colnames(priors_matrix) <-c("mu_norm", "sigma_norm", "slope_gamma", "var_gamma")

SED4=1478963

# for stan
nchains=4


#### data reading
data <- read.csv("irinotecan.csv")

# derive data characteristics:
N         <- nrow(data)              # number of data points (study/dose combinations)
studyIdx  <- as.numeric(data$study)  # study indices
studyName <- levels(data$study)      # study names
K         <- max(studyIdx)           # number of studies

# show some study-level summary statistics:
sumstats <- cbind("doses"=table(data$study),
                  "patients"=tapply(data$total, data$study, sum),
                  "events"=tapply(data$events, data$study, sum))
print(sumstats)


doses_true <<- sort(unique(data$dose))    #in mg for example
doses_modif <<- doses_true/mean(doses_true)
unit_dose= doses_modif[2]-doses_modif[1] 
# specify data provided to Stan

delta_dose <- (doses_modif[2:length(doses_modif)] - doses_modif[1:(length(doses_modif)-1)]) / unit_dose


nDLT_meta <- function(doses, data) {
  sum(data[data$dose == doses, "events"])
}
nPT_meta <- function(doses, data) {
  sum(data[data$dose == doses, "total"])
}
ti <- round(sapply(doses_true, nDLT_meta, data = data), 0)
ni <- round(sapply(doses_true, nPT_meta, data = data), 0)

EB_prior= pava(ti/ni)
pos_EB = order(abs(EB_prior-0.33))[1]
pos_prior = ( ((doses_modif[pos_EB] - doses_modif[1])/unit_dose) >=2 )*1


doselevel_tot = doses_true
I = length(doselevel_tot)
doselevel=match(data$dose, doselevel_tot)
events =  data$events
total = data$total

var_coeff= priors_matrix[pos_prior+1, "var_gamma"]
slope_gamma = priors_matrix[pos_prior+1, "slope_gamma"]
mu_first = priors_matrix[pos_prior+1, "mu_norm"]
sd_first = priors_matrix[pos_prior+1, "sigma_norm"]

shape_gamma = 1/var_coeff^2
rate_gamma = (var_coeff^2*slope_gamma)^-1

stanData <- list(K        = K,
                 N       = N,
                 I       = I,
                 studyIdx = studyIdx,
                 doselevel     = doselevel,
                 events  = events,
                 total   = total,
                 doses_modif = doses_modif,
                 shape_gamma  = shape_gamma,
                 rate_gamma  = rate_gamma,
                 delta_dose = delta_dose,
                 mu_first = mu_first,
                 sd_first = sd_first)

set.seed(SED4)
# only to see that with the wrong parametrization, you have some issues in sampling
#fit <- sampling(sm_meta, data=stanData ,iter=4000, chains=4, control = list(adapt_delta = 0.9))

fit_mod <- sampling(sm_meta, data=stanData ,iter=4000, chains=nchains, control = list(adapt_delta = 0.95))
results <- extract(fit_mod, pars="ptox")$ptox




ptox_res <- apply(results,2,median)
ptox_res2 <- apply(results,2,mean)
q1_res<- apply(results,2,quantile, probs=0.25)
q3_res <- apply(results,2,quantile, probs=0.75)



results_Irinotecan = results

doses_true <- c(40,50 ,60,70,80, 90, 100, 120, 125, 150)

forplot <- NULL
for (i in 1:dim(results_Irinotecan)[2]) forplot <- rbind(forplot, cbind(results_Irinotecan[,i], 
                                                                       rep(doses_true[i],dim(results_Irinotecan)[1] )))
dimnames(forplot)[[2]] <- c("Probability", "Dose")
forplot <- data.frame(forplot)
#forplot$Dose <- as.factor(forplot$Dose)

p2 <- ggplot(forplot, aes(x=Dose, y=Probability, group=factor(Dose))) +  
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     plot.title = element_text(size=20, face="bold", hjust = 0.5),
                     axis.title.x = element_text(size=14, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = element_text(size=14, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = element_text(face = "bold", size = 12),
                     axis.text.y = element_text(face = "bold", size = 12),
                     legend.position = "none") +
  ggtitle("Irinotecan + S1") +
  xlab(expression(bold(dose~(mg/m^2)))) + ylab("Probability of toxicity") +
  geom_violin(trim=FALSE, size=1.1, fill="#D3D3D3") +
  geom_hline(yintercept=0.33, color = "black", size=0.8, linetype="dashed") +
  geom_hline(yintercept=0.25, color = "black", size=0.8, linetype="dashed") +
  geom_hline(yintercept=0.20, color = "black", size=0.8, linetype="dashed") +
  coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks=c(0,0.20,0.25,0.33,0.5,0.75,1)) +
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  scale_x_continuous(breaks = doses_true )

p2

jpeg("plot_Irinotecan.jpeg", width = 15, height = 6, unit="in", res=300)
p2
dev.off()


results_Ewok = results_Irinotecan

checking <- function(x, tau) mean(x>=tau)

round(apply(results_Ewok, 2, checking, tau=0.33),3)
# 0.000 0.000 0.000 0.004 0.061 0.349 0.773 0.990 0.996 1.000
round(apply(results_Ewok, 2, checking, tau=0.25),3)
# 0.000 0.000 0.002 0.027 0.238 0.677 0.944 0.998 1.000 1.000
round(apply(results_Ewok, 2, checking, tau=0.20),3)
# 0.000 0.001 0.008 0.082 0.466 0.866 0.984 1.000 1.000 1.000
