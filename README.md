
Here, you can find the code to run the models and the simulations presented in the article entitled “Random-effects meta-analysis of Phase I dose-finding studies using stochastic process priors”. 
Link: https://clicktime.symantec.com/3MnsciFmpCx2GTRRtjmLJoD6H2?u=https%3A%2F%2Farxiv.org%2Fabs%2F1908.06733

The scripts can be run in all operating system. R must be installed with the following libraries: 
    rstan (we used rstan_2.17.3), UBCRM, MASS, dfcrm, Iso, ggplot2, ggrepel, gridExtra, parallel, dfmeta, xtable.  

We worked on the following environment:
R version 3.5.0 
Platform: macOS High Sierra 10.13.6

Folders:

*	**stan_models**: the Stan models needed to compute all proposed methods are grouped in this folder. The user must to add them in the same folder of the R analysis session. 
    * MADF.stan;
    * MADF1.stan;
    * MADF2.stan;
    * MADF3.stan (to note, this is used also for MADF4). 

*	**data**: data used in the work are collected in this folder. 
    * sorafenib2.csv;
    * irinotecan.csv.

*   **scripts**: R scripts used to run methods and simulations.
    * scenario_generation.r: to generate datasets in a pre-specified scenario (see details in the script);
    * scenario_npatients.r: to count the number of patients allocated to each dose at each run;
    * MADF_sim.R: to run simulation using MADF method; 
    * MADF1_sim.R: to run simulation using MADF1 method; 
    * MADF2_sim.R: to run simulation using MADF2 method; 
    * MADF3-4_sim.R: to run simulation using MADF3 and MADF4 method; 
    * ZKO_sim.R: to run simulation using ZKO method; 
    * read_results.R: to create the results table with percentage of final MTD selection; 
    * plot_prior_induced.R: to run plot the induced priors and ESS on fixed effects; 
    * Estimation_Examples.R: to run MADF on the datasets in data folder. 


For any issue, please contact moreno.ursino@gmail.com.
