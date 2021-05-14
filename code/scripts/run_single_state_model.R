# R Script for running the single-region model
# at the **state** level

##############################
# Load packages

library(readr)
library(dplyr)
library(stringr)
library(splines)
library(tidyr)
library(rstan)
library(scriptmisc) # devtools::install_github("jpkeller/scriptmisc")

#############################################
# Preliminary Setup

rm(list=ls())
Sys.setenv(USE_CXX14 = 1) # Necessary for rstan

results_dir <- "results/results_states/"
data_dir <- "data/"
stan_dir <- "code/stan/"

#############################################
# Read in the command-line arguments

## Reading in the arguments
args <- commandArgs(TRUE)
print(args)

## Parse the arguments
if(length(args)==0){
    print("No arguments supplied.")
} else {
    for (i in 1:length(args)){
        eval(parse(text=args[i]))
    }
}

# Expected parameters passed as command-line arguments and their defaults:

# Parameters for Data Setup
data_setup_defaults <- list(model_state="colorado") # State for model fit)

# Parameters for Model Fitting
model_fitting_defaults <- list(stan_model="seiurd_fixEpsSD_betaepi", # STAN model name
                               model_name_string="",
                               start_date_cases=100, # First day of model is day on which this number of cases  (or more) reached
                               start_date="2020-03-17", # First day of data to use (only if start_date_cases is 0)
                               stop_date=NULL, # Last day of data to use
                               betaU_d=10, #df for betaU spline
                               betaU_mobility=FALSE, # If TRUE, include mobility data
                               mobility_smooth="none",
                               betaU_degree=1, # degree of polynomial spline for betaU
                               betaU_natspline=FALSE, # use natural cubic splines
                               betaU_interval=30, # interval of piecewise constant components
                               betaU_addinterceptbasis=FALSE, # If TRUE, betaU basis includes intercept
                               omega_d=5, # df for omega spline
                               omega_degree=1, # degree of polynomial spline for betaU
                               omega_natspline=FALSE, # use natural cubic splines
                               omega_interval=30, # interval of piecewise constnat components
                               omega_addinterceptbasis=FALSE, # If TRUE, omega basis includes intercept
                               psi_d=7, # df for eta spline
                               psi_degree=1, # degree of polynomial spline for betaU
                               psi_natspline=FALSE, # use natural cubic splines
                               psi_interval=30, # interval of piecewise constnat components
                               psi_addinterceptbasis=FALSE # If TRUE, omega basis includes intercept
)


# Prior Distribution Parameters
prior_defaults <- list(betaU_b0_mean=0,
                       betaU_b0_sd=1,
                       betaU_b_mu=0,
                       betaU_b_sigma_mean=0.5, #0.3,
                       betaU_b_sigma_sd=0.1, # 0.05,
                       tau_a=1,
                       tau_b=2, #6
                       psi_a_a=3,
                       psi_a_b=1,
                       psi_a_mean=NULL,
                       psi_a_sd=NULL,
                       psi_a_limit_low=0,
                       psi_a_limit_high=1,
                       psi_b0_a=1,
                       psi_b0_b=1,
                       psi_c_a=5,
                       psi_c_b=100,
                       psi_b0_mean=0,
                       psi_b0_sd=1,
                       psi_b_mu=0,
                       psi_b0_limit_high=4,
                       psi_b0_limit_low=-4,
                       psi_b_limit_high=4,
                       psi_b_limit_low=-4,
                       psi_b_sigma_limit=1,
                       psi_b_sigma_mean=0.5,
                       psi_b_sigma_sd=0.1,
                       omega_b0_mean=-1,
                       omega_b0_sd=1,
                       omega_b_mu=0,
                       omega_b_sigma_mean=0.5, 
                       omega_b_sigma_sd=0.1,
                       rescale_prior_sd=FALSE,
                       alpha_a=12,
                       alpha_b=22.5,
                       alpha_mean=NULL,
                       alpha_sd=NULL,
                       alpha_limit_low=0,
                       alpha_limit_high=1,
                       eta0_a=20,
                       eta0_b=26,
                       eta0_mean=NULL, #0.4,
                       eta0_sd=NULL, #0.1,
                       eta0_limit_low=0,
                       eta0_limit_high=1,
                       nu_a=120, #7,
                       nu_b=720, #42,
                       nu_mean=NULL, #0.14,
                       nu_sd=NULL, #0.1,
                       nu_limit_low=0,
                       nu_limit_high=1,
                       gamma0_a=300, #25,
                       gamma0_b=5010, #275,
                       gamma0_mean=NULL, #0.1,
                       gamma0_sd=NULL, #0.1,
                       gamma0_limit_low=0,
                       gamma0_limit_high=1,
                       delta0_a=300, # 13,
                       delta0_b=2940, # 156,
                       delta0_mean=NULL, #0.1,
                       delta0_sd=NULL, #0.1,
                       delta0_limit_low=0,
                       delta0_limit_high=1,
                       EEIU0_a=1,  #20
                       EEIU0_b=1,  # 80
                       odC_a=50, #1,
                       odC_b=100, #1,
                       odD_a=50,# 1,
                       odD_b=100, # 1
                       phiC_a=50, #1,
                       phiC_b=100, #1,
                       phiD_a=50,# 1,
                       phiD_b=100, # 1
                       kappa_a=25,#125,
                       kappa_b=5,#25,
                       epsilonI_sd=0.1,
                       epsilonD_sd=0.1,
                       epsilonI_sd_a=1, #5,
                       epsilonI_sd_b=10, #100
                       epsilonD_sd_a=1,
                       epsilonD_sd_b=10
) #1)




# Algorithm Parameters
algorithm_defaults <- list(iter=20,
                           warmup_iter=20,
                           nchains=1, 
                           ncores=1,
                           adapt_delta=0.9,
                           treedepth=12,
                           save_warmup=FALSE)

vars_defaults <- c(data_setup_defaults,
                   model_fitting_defaults,
                   prior_defaults,
                   algorithm_defaults)

# Set variables to their defaults, if they aren't provided
sapply(names(vars_defaults),
       check_and_set_Default,
       defaults= vars_defaults,
       envir=parent.frame())
#############################################
# Rescale the priors

# Only done if rescale_prior_sd is a  number not FALSE
# *and* the standard deviation parameters are provided.
if (rescale_prior_sd){
  if (is.null(alpha_sd) || is.null(eta0_sd) || is.null(nu_sd) || is.null(gamma0_sd) || is.null(delta0_sd)  || is.null(psi_a_sd) ) stop("rescaling but no sd provided.")
  alpha_sd <- alpha_sd*rescale_prior_sd
  eta0_sd <- eta0_sd*rescale_prior_sd
  nu_sd <- nu_sd*rescale_prior_sd
  gamma0_sd <-gamma0_sd*rescale_prior_sd
  delta0_sd <- delta0_sd*rescale_prior_sd
  psi_a_sd <- psi_a_sd*rescale_prior_sd
}

recalc_beta_a <- function(mean, sd){
  a=(1-mean)*mean^2/sd^2 - mean
  a
}
recalc_beta_b <- function(mean,a){
  b=a*(1/mean -1)
  b
}

if (!is.null(alpha_sd)){
  alpha_a <- recalc_beta_a(mean=alpha_mean,
                sd=alpha_sd)
  alpha_b <- recalc_beta_b(mean=alpha_mean,
                           a=alpha_a)
}
if (!is.null(eta0_sd)){
  eta0_a <- recalc_beta_a(mean=eta0_mean,
                           sd=eta0_sd)
  eta0_b <- recalc_beta_b(mean=eta0_mean,
                           a=eta0_a)
}
if (!is.null(nu_sd)){
  nu_a <- recalc_beta_a(mean=nu_mean,
                           sd=nu_sd)
  nu_b <- recalc_beta_b(mean=nu_mean,
                           a=nu_a)
}
if (!is.null(gamma0_sd)){
  gamma0_a <- recalc_beta_a(mean=gamma0_mean,
                           sd=gamma0_sd)
  gamma0_b <- recalc_beta_b(mean=gamma0_mean,
                           a=gamma0_a)
}
if (!is.null(delta0_sd)){
  delta0_a <- recalc_beta_a(mean=delta0_mean,
                           sd=delta0_sd)
  delta0_b <- recalc_beta_b(mean=delta0_mean,
                           a=delta0_a)
}
if (!is.null(psi_a_sd)){
  psi_a_a <- recalc_beta_a(mean=psi_a_mean,
                           sd=psi_a_sd)
  psi_a_b <- recalc_beta_b(mean=psi_a_mean,
                           a=psi_a_a)
}

#############################################
# Read in the data

all_cases <- read_csv(paste0(data_dir, "cases_deaths_allstates_preprocessed.csv"))

###########################################
# Setup the Data

# Case/Death Data
state_cases <- subset(all_cases,
                      state==model_state)
state_population <- state_cases$state_population[1]

# Start date is first day after cases over 100
if (start_date_cases>0){
  start_date <- state_cases %>%
    filter(cases>start_date_cases) %>%
    filter(date==min(date)) %>%
    .$date
}

# Take early case numbers for initial conditions
early_case_count <- state_cases %>%
  filter(date < start_date & date >= start_date - 5) %>%
  .$new_cases %>%
  sum()

# sum(state_cases$new_cases[state_cases$date<start_date])
state_cases <- subset(state_cases,
                      date>=start_date) %>%
  arrange(date)

if (!is.null(stop_date)){
  state_cases <- subset(state_cases,
                        date<=stop_date) %>%
    arrange(date)
}
current_dates <- state_cases$date

# Splines
# Set up splines for time-varying parameters
#
if (!betaU_natspline){
  if (betaU_degree>0) {
    betaU_Mx <- bs(current_dates,
                   degree=betaU_degree,
                   df=betaU_d,
                   intercept=betaU_addinterceptbasis)
    # betaU_Mx <- scale(betaU_Mx, center=FALSE)
  } else if (betaU_degree==0) {
    n_betaU_pieces <- as.numeric(max(current_dates) - min(current_dates)) %/% betaU_interval
    n_betaU_pieces <- n_betaU_pieces + 1 # add the final partial period
    betaU_Mx <- model.matrix(~cut(current_dates,
                                  breaks=as.Date(start_date) + betaU_interval*(0:n_betaU_pieces),
                                  include.lowest=TRUE))
    betaU_Mx <- betaU_Mx[, -1, drop=FALSE]
    betaU_d <- ncol(betaU_Mx)
  }
} else {
  betaU_Mx  <- ns(current_dates,
                  df=betaU_d,
                  intercept=betaU_addinterceptbasis)
  betaU_Mx <- scale(betaU_Mx)
}
# Mobility data
if (betaU_mobility!=FALSE){
  # Mobility data
  betaU_mobility <- str_split(betaU_mobility, "\\.")[[1]]
  for (i in 1:length(betaU_mobility)){
    if (betaU_mobility[i] %in% c("sgFT","sgPT", "sgHD", "sgCH", "sgRES", "aggLEX")){
      if (substr(betaU_mobility[i], 1, 2)=="sg"){
        state_mobility <- read_csv(paste0(data_dir, "state_mobility_sg_smoothed.csv")) %>%
          mutate(state=str_replace(state, pattern = " ", replacement = "")) %>%
          mutate(state=str_to_lower(state))
  
        state_mobility <- state_mobility %>%
          filter(state==model_state)
        # Fill in missing dates via merge
        
        if (mobility_smooth %in% c("roll7", "ks", "sm")){
          state_mobility <- state_mobility %>%
            select(date,
                   mob=paste0(betaU_mobility[i], "_", mobility_smooth))
        } else if (mobility_smooth=="none") {
          mob_var <- case_when(betaU_mobility[i]=="sgCH"~"completely_home_prop",
                               betaU_mobility[i]=="sgFT" ~"full_time_work_prop",
                               betaU_mobility[i]=="sgPT" ~"part_time_work_prop",
                               betaU_mobility[i]=="sgHD" ~"median_home_dwell_time",
                               betaU_mobility[i]=="sgRES"~"restaurants_visit_prop")
          state_mobility <- state_mobility %>%
            select(date, mob=mob_var) %>%
            zoo::na.locf()
        } else {
          stop("Selected mobility smooth not currently implemented.")
        }
      } else if (betaU_mobility[i]=="aggLEX") {
        state_mobility <- read_csv(paste0(data_dir, "state_from_county_mobility_aggLEX_smoothed.csv")) %>%
          mutate(state=str_replace(state, pattern = " ", replacement = ""))
        state_mobility <- state_mobility %>%
          filter(state==model_state)
        if (mobility_smooth %in% c("roll7", "ks", "sm")){
          state_mobility <- state_mobility %>%
            select(date,
                   mob=paste0("aggLEX_", mobility_smooth))
        } else if (mobility_smooth=="none") {
          # Carry forward for missing values
          state_mobility <- state_mobility %>%
            select(date, mob=aggLEX) %>%
            zoo::na.locf()
          
        } else {
          stop("Selected mobility smooth not currently implemented.")
        }
        
      }
      state_mobility <- left_join(state_cases[, "date"], state_mobility)
      betaU_Mx<- cbind(betaU_Mx, scale(state_mobility$mob))
      betaU_d <- betaU_d + 1
    }
  }
}
if (!omega_natspline){
  if (omega_degree>0) {
    omega_Mx <- bs(current_dates,
                   degree=omega_degree,
                   df=omega_d,
                   intercept=omega_addinterceptbasis)
  } else if (omega_degree==0) {
    n_omega_pieces <- as.numeric(max(current_dates) - min(current_dates)) %/% omega_interval
    n_omega_pieces <- n_omega_pieces + 1 # add the final partial period
    omega_Mx <- model.matrix(~cut(current_dates,
                                  breaks=as.Date(start_date) + omega_interval*(0:n_omega_pieces),
                                  include.lowest=TRUE))
    omega_Mx <- omega_Mx[, -1, drop=FALSE]
    omega_d <- ncol(omega_Mx)
  }} else {
    omega_Mx  <- ns(current_dates,
                    df=omega_d,
                    intercept=omega_addinterceptbasis)
    omega_Mx <- scale(omega_Mx)
  }
if (!psi_natspline){
  if (psi_degree>0) {
    psi_Mx <- bs(current_dates,
                 degree=psi_degree,
                 df=psi_d,
                 intercept=psi_addinterceptbasis)
  } else if (psi_degree==0) {
    n_psi_pieces <- as.numeric(max(current_dates) - min(current_dates)) %/% psi_interval
    n_psi_pieces <- n_psi_pieces + 1 # add the final partial period
    psi_Mx <- model.matrix(~cut(current_dates,
                                breaks=as.Date(start_date) + psi_interval*(0:n_psi_pieces),
                                include.lowest=TRUE))
    psi_Mx <- psi_Mx[, -1, drop=FALSE]
    psi_d <- ncol(psi_Mx)
  }}  else {
    psi_Mx  <- ns(state_cases$date,
                  df=psi_d,
                  intercept=psi_addinterceptbasis)
    psi_Mx <- scale(psi_Mx)
  }

# State Data
state_stan_data <- list(yc=pmax(state_cases$new_cases,0),
                        yd=pmax(state_cases$new_deaths,0),
                        nT=length(current_dates),
                        date=current_dates,
                        betaU_d=ncol(betaU_Mx),
                        betaU_Mx=betaU_Mx,
                        psi_d=ncol(psi_Mx),
                        psi_Mx=psi_Mx,
                        omega_d=ncol(omega_Mx),
                        omega_Mx=omega_Mx)

# Initial State
E_initial <- 0 # set in STAN code
IU_initial <- 0 # set in STAN code
ID_initial <- early_case_count
S_initial <- state_population - (E_initial + IU_initial + ID_initial)
state_init <- c(S_initial, #S
                E_initial, #E
                IU_initial, #IU
                0, # RU
                ID_initial, # ID
                0, # UD
                0, #RD
                0) #DD


# Prior distributions
current_prior <- list(tau_a=tau_a,
                      tau_b=tau_b,
                      betaU_b0_mean=betaU_b0_mean,
                      betaU_b0_sd=betaU_b0_sd,
                      betaU_b_mu=betaU_b_mu,
                      betaU_b_sigma_mean=betaU_b_sigma_mean,
                      betaU_b_sigma_sd=betaU_b_sigma_sd,
                      psi_a_limit_low=psi_a_limit_low,
                      psi_a_limit_high=psi_a_limit_high,
                      psi_a_a=psi_a_a,
                      psi_a_b=psi_a_b,
                      psi_a_mean=psi_a_mean,
                      psi_a_sd=psi_a_sd,
                      psi_b0_a=psi_b0_a,
                      psi_b0_b=psi_b0_b,
                      psi_c_a=psi_c_a,
                      psi_c_b=psi_c_b,
                      psi_b0_mean=psi_b0_mean,
                      psi_b0_sd=psi_b0_sd,
                      psi_b_mu=psi_b_mu,
                      psi_b_sigma_mean=psi_b_sigma_mean,
                      psi_b_sigma_sd=psi_b_sigma_sd,
                      psi_b_sigma_limit=psi_b_sigma_limit,
                      psi_b_limit_low=psi_b_limit_low,
                      psi_b_limit_high=psi_b_limit_high,
                      psi_b0_limit_low=psi_b0_limit_low,
                      psi_b0_limit_high=psi_b0_limit_high,
                      omega_b0_mean=omega_b0_mean,
                      omega_b0_sd=omega_b0_sd,
                      omega_b_mu=omega_b_mu,
                      omega_b_sigma_mean=omega_b_sigma_mean,
                      omega_b_sigma_sd=omega_b_sigma_sd,
                      epsilonI_sd=epsilonI_sd,
                      epsilonD_sd=epsilonD_sd,
                      epsilonI_sd_a=epsilonI_sd_a,
                      epsilonI_sd_b=epsilonI_sd_b,
                      epsilonD_sd_a=epsilonD_sd_a,
                      epsilonD_sd_b=epsilonD_sd_b,
                      alpha_a=alpha_a,
                      alpha_b=alpha_b,
                      alpha_mean=alpha_mean,
                      alpha_sd=alpha_sd,
                      alpha_limit_low=alpha_limit_low,
                      alpha_limit_high=alpha_limit_high,
                      nu_a=nu_a,
                      nu_b=nu_b,
                      nu_mean=nu_mean,
                      nu_sd=nu_sd,
                      nu_limit_low=nu_limit_low,
                      nu_limit_high=nu_limit_high,
                      eta0_a=eta0_a,
                      eta0_b=eta0_b,
                      eta0_mean=eta0_mean,
                      eta0_sd=eta0_sd,
                      eta0_limit_low=eta0_limit_low,
                      eta0_limit_high=eta0_limit_high,
                      gamma0_a=gamma0_a,
                      gamma0_b=gamma0_b,
                      gamma0_mean=gamma0_mean,
                      gamma0_sd=gamma0_sd,
                      gamma0_limit_low=gamma0_limit_low,
                      gamma0_limit_high=gamma0_limit_high,
                      delta0_a=delta0_a,
                      delta0_b=delta0_b,
                      delta0_mean=delta0_mean,
                      delta0_sd=delta0_sd,
                      delta0_limit_low=delta0_limit_low,
                      delta0_limit_high=delta0_limit_high,
                      kappa_a=kappa_a,
                      kappa_b=kappa_b,
                      EEIU0_a=EEIU0_a,
                      EEIU0_b=EEIU0_b,
                      odC_a=odC_a,
                      odC_b=odC_b,
                      odD_a=odD_a,
                      odD_b=odD_b,
                      phiC_a=phiC_a,
                      phiC_b=phiC_b,
                      phiD_a=phiD_a,
                      phiD_b=phiD_b)



# Data for STAN model
stan_data <- c(state_stan_data,
               state_init=list(state_init),
               current_prior)

###########################
# Model Fitting

filename_stem <- paste0(model_state, "_", stan_model, ifelse(nchar(model_name_string)>0, "_",  ""), model_name_string)
filename_fit <- paste0("fit_", filename_stem, ".RData")
filename_data <- paste0("data_", filename_stem, ".RData")

# Compile STAN Model
stan_model_obj <- stan_model(paste0(stan_dir, stan_model, ".stan"),  isystem=c(getwd(),paste0(getwd(), "/code/stan/")))



# Initial Values for parameters
initf1 <- function() {
  list(tau=runif(1, 0.1, 0.2),
       alpha=runif(1, 0.2, 0.4),
       # nu=runif(1, 0.12, 0.15),
       gamma0=gamma0_limit_low + (gamma0_limit_high - gamma0_limit_low)*runif(1, 0.1, 0.5),
       delta0=delta0_limit_low + (delta0_limit_high - delta0_limit_low)*runif(1, 0.05, 0.1),
       betaU_b0_raw=-2, # -2,
       betaU_b_raw = rep(0, betaU_d),
       betaU_b_sigma=runif(1, 0.4, 0.6),
       eta0=eta0_mean,
       psi_c=0.1,
       psi_b0_raw=-1,
       psi_b_raw=rep(0, psi_d),
       psi_b_sigma=runif(1, 0.4, 0.5),
       omega_b0_raw=0,
       omega_b_raw =rep(0, omega_d),
       omega_b_sigma=runif(1, 0.4, 0.6),
       epsilonI_sd=0.1,
       epsilonD_sd=0.1,
       epsilonI=rep(1, state_stan_data$nT),
       epsilonD=rep(1, state_stan_data$nT-1),
       extra_count=0,
       kappa=1,
       EEIU0=0.2,
       odC=0.3,
       odD=0.3)
}

save(stan_data, file=paste0(results_dir,filename_data ))


fit <- sampling(stan_model_obj,
                data=stan_data,
                chains=nchains,
                iter=iter + warmup_iter,
                warmup=warmup_iter,
                save_warmup=save_warmup,
                cores=ncores,
                control=list(max_treedepth=treedepth,
                             adapt_delta=adapt_delta),
                init=initf1,
                verbose=TRUE,
                init_r=0.4)

save(fit, file=paste0(results_dir, filename_fit))

print(fit,
      pars=c("muC", "muD", "muCorig", "muDorig", "state",  "eta", "betaU", "theta", "epsilonI", "epsilonD", "gamma", "psi", "rho",  "gammaU", "gammaD", "delta", "deltaD" ,"phiC", "phiD", "R0","R0_noepsilon", "omega", "yc_pred", "yd_pred"),
      include=FALSE)

print(fit,
      pars=c("R0[1]","R0[10]"),
      include=T)
