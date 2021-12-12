# COVID-19-Nonparametric-Inference (https://doi.org/10.1101/2021.02.01.21250936)

The code for the multi-compartment stochastic dynamic model,
described in tha article "Inference on the dynamics of the COVID pandemic from observational data" and 
applied to the evolution of the COVID-19 pandemic in various states of USA, is available here. The model is decribed in the Figure "diagram.png" and the File "model.pdf".

The proposed  model is expressed through a system of difference equations that describe their temporal dynamics.
Information on the social distancing measures and diagnostic testing rates are incorporated to characterize the dynamics of the various compartments of the model.
The model involves interpretable temporally static and dynamic epidemiological rate parameters. 
A model fitting strategy built upon nonparametric smoothing is employed for estimating the time-varying parameters, while profiling over the time-independent parameters.
Confidence bands of the parameters are obtained through a residual bootstrap procedure. 
A key feature of the methodology is its ability to estimate latent unobservable compartments such as the number of asymptomatic but infected individuals who are known to be the key vectors of COVID-19 spread.
The nature of the disease dynamics is further quantified by relevant epidemiological markers that make use of the estimates of latent compartments.
The methodology is applied to understand the true extent and dynamics of the pandemic in various states within the USA. 

However, the same model can be applied to describe the population dynamics for any state/country/city as long as the data on the number of cases, deaths, current hospitalizations, recovery and some quantification of the social distancing measure taken in the particular context is available. 


# Description of the files

(A) The folder "USA data processing" contains the data processing steps: 

1. The number of cases, deaths are obtained from https://github.com/CSSEGISandData/COVID-19.git. These are stored in the folder ```"data" -> "COVID-19"```.
2. The data on current hospitalization, and recovery are obtained from https://covidtracking.com/data/download. These are stored in the folder ```"data" -> "covid_tracking"``` (according to the USA states).
3. The social mobility data in terms of reduced (/increased) activity in different sectors are obtained from https://www.google.com/covid19/mobility. 
These are stored in the folder ```"data" -> "EconomicTracker"``` (cloned from the Github Repository https://github.com/OpportunityInsights/EconomicTracker.git)
4. The above data are merged into a data-frame in R using the code in ```"state_data_processing.R"``` and the output is stored in the RData file ```"covid_df.Rda"```.

As long as the relevant data on the number of cases, deaths, hospitalization, recovery, and social mobility are available, the same steps in the R code can be retraced.

(B) The folder ```"code"``` contains the estimation procedure:

1. ```fitting_functions.R``` -- Code for optimizing the profile loss function to estimate γ and ρ<sub>A</sub>.
2. ```est_plot_func.R``` -- Code for estimation of the other parameters  for a given value of γ and ρ<sub>A</sub>. 
Also contains the ggplot codes for plotting the estimated values (displays the goodness of the fits).
3. ```est_module.R ```-- The function "estimation_module" takes in observed dataframe (consisting of the information on the number of cases, deaths, hospitalization, recovery, and social mobility); 
the tuning parameters ("numBlock" and "wd.length") and returns the estimated parameters. It also outputs the estimates for the compartments C,Q,R,D,H.
Based on the input argument "boot == FALSE" (i.e., based on the observed sample, not the bootstrapped sample), the ggplots of the profile loss vs gamma.grid and profile loss vs rhoA.grid are saved.
The function "plot_for_estimation" runs the function estimation_module for boot == FALSE (observed data only) and saves the estimated parameters and relevant plots.
4. ```individual_plots.R``` -- Code for running the estimation module for any given state in the USA and saving the respective plots.  
5. ```resi_boot.R``` -- Code for performing residual bootstrap method. Saves the bootstrap estimates for the time independent and time varying parameters 
as well as the estimated states C,Q,R,H,D (as an output of the function estimation_module) for each bootstrap sample into .Rda files.
It calls the estimation module function using the argument "boot == TRUE". Requires parallel computing to generate the bootstrap samples for any given state of the USA. 
6. ```ggplot_ci.R``` -- Based on the output of the "resi_boot.R", plots the CIs for different parameters across different states. 

7. ```opt_wd.R```-- Computes the optimal tuning parameters ("wd.length" and "numBlock") for carrying out the nonparametric estimation process in "fitting_functions.R" and "est_plot_func.R". It also stores the optimal parameters for a given state.

(C) Aggregate Simulation

A detailed numerical simulation under the Poisson process model is also carried out to validate the statistical performance of the proposed estimation procedure.
The proposed method and estimation procedure do not explicitly use the underlying assumption of a Poisson process. 
However, we use an ensemble of independent Poisson processes to simulate data from the proposed model. 
These aggregated data sets are then used to accurately estimate various parameters, which validate our estimation procedure. 
The aggregation has the effect of increasing the number of observations in the compartments and thereby improving estimation accuracy.
If the number of individuals in the symptomatic or quarantined compartments are low, e.g. at the onset of the pandemic, inherent biases are introduced 
in the estimated trajectories. A bigger sample size is required to correct such contaminants.

The results and code are found in the folder ```"aggregation_simulation"```.


# How to run the code :

* Step 1 - Use the data processing code (```"state_data_processing.R"``` in the folder "USA data processing" ) to generate the ```"covid_df.Rda"``` file.

* Step 2 - Find the optimal tuning paramters ("numBlock" and "wd.length") for one/multiple state(s) in the USA that is of interest, using ```"opt_wd.R" ```(in the folder "code").

* Step 3 - For the estimation of different parameters and compartments of the model for one/multiple state(s) in the USA that is of interest, run the ```"individual_plots.R"``` script (in the folder "code"). 
This R script requires the tuning paramters to be provided for the given US state (the optimal paramters are found in step 2) and sources the file ```"est_module.R"```, 
which in turn uses the code from "fitting_functions.R" and "est_plot_func.R". 
The relevant plots are generated and the estimation outputs are stored in the RData file ```"<state name>_est_st.Rda"```.

* Step 4 - For performing the residual bootstrap, use ```"resi_boot.R"``` (in the folder "code"). This requires the estimation outputs from Step 3.
Bootstrap estimate outputs are stored in the RData files ```"<state name>_est_bs.Rda"``` and ```"<state name>_est_boot_data.Rda"```.

* Step 5 - For generating the confidence interval and histogram plots, use ```"ggplot_ci.R"``` (in the folder "code"). This requires the output files from Step 3 and Step 4.

For example: for the state ```st = "Utah"```, the optimal tuning parameters found using Step 1 and Step 2 are stored in opt_wd_list.csv. For finding the estimators  for the relevant epi-markers, and save the relevant plots we run the wrapper function ```"individual_plots.R"``` in Step 3. The boostrap CIs and histogram plots can be generated using the wrapper function ```"ggplot_ci.R"``` in Step 5 (which uses the boostrapped data and estimates from Step 4, which are already stored in RData files ```"Utah_est_bs.Rda"``` and ```"Utah_est_boot_data.Rda"``` in the folder "Utah_updated_boot_final").


# Installing dependencies

The following packages are required additionally: 
```
library(magrittr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(plyr)
library(gridExtra)
library(parallel)
```
# Results

The estimation results and relevant plots are stored in the folder "Results for USA States" according to the state names. These consists of the following:

1. ``` <state name>_est_st.Rda ``` -- contains the estimated parameters and compartments for each state called "state name".

2. Contains plots for the estimated parameters and relevant epi-markers.


3. ```<state name>_est_bs.Rda ``` -- contains the estimated parameters using bootstrap samples for each state called "state name" (in the corresponding folder "<state name>_updated_boot_final").

  
4. ```<state name>_est_boot_data.Rda ``` -- contains the estimated compartments using bootstrap samples for each state called "state name" (in the corresponding folder "<state name>_updated_boot_final").

5. Contains plots for the confience intervals for the relevant paramters and relevant epi-markers, generated using the residual bootstrap method (in the corresponding folder "<state name>_updated_boot_final").



For example, for the US state, ```<state name> == Utah```, relevant estimation plots are stored in the sub-folder named "Utah" (in the folder "Results for USA States").
The corresponding Rdata file containing the estimation results using the observed data is named  as ```Utah_est_st.Rda```. The corresponding revised bootstrap code is stored in ```Utah_est_bs.Rda```,  and ```Utah_est_boot_data.Rda ```, respectively, in the folder "Utah__updated_boot_final".
The relevant plots for the bootstrap confidence intervals are stored in the folder "Results for USA States" -> "Utah" -> "Utah__updated_boot_final" -> "CI_Utah". 



