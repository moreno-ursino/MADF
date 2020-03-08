##############################################################
#
# Script to extract the number of patients
# at each dose for each run in the selected scenario
#
#
##############################################################


# run "scenario_generation.R" and then the following code


n_pat_scenario <- NULL
for (tr in 1:TR){
  data_analysis <- datagen[[tr]]
  n_pat_trial = NULL
  for (i in p_true){
    index = which(data_analysis$res[,1]==i)
    if (length(index)==0) {
      n_pat_trial = c(n_pat_trial, 0)
    } else {
      n_pat_trial = c(n_pat_trial, sum(data_analysis$res[index,4]))
    }
  }
 n_pat_scenario = rbind(n_pat_scenario, n_pat_trial)
}

#table_res=apply(n_pat_scenario, 2, quantile, probs=c(0.5, 0.25, 0.75))

write.table(n_pat_scenario, file=paste("patient_scen",scen,"_",NTcrm+NT3p3,".txt", sep=""), row.names=FALSE,
            col.names = FALSE)
