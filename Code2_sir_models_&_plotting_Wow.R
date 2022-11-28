############## NOTES ##################
### SIR Ebola virus disease model - apes
# Based on models in Keeling & Rohani, 2011
# & on code from: https://cran.r-project.org/web/packages/EpiDynamics/EpiDynamics.pdf
#
# the code is roughly set up in chunks as follows:
# load packages, set up some defaults ...
# then the models are written, with small test runs & plots to check they run
# then the models are run with some 'default' parameters, first singly, then with multiple runs
# then there are some functions to pull out metrics and make plots
# then the models are run across the full range of parameters
# then there's a bit of re-coding for the meta-population models
#
############## LOAD PACKAGES ##########

library(reshape)
library(EpiDynamics)
library(dplyr)     
library(reshape2) 
library(stringr)
library(emdbook)  
library(ggplot2); theme_set(theme_bw())

## MODELS ################################

############################## set up base values/ranges and empty arrays ###
# based on the pathogens:
beta_a <- seq(from = 1,to = 20, by = 2)
rho_a <- seq(from = 0, to = 0.9, by = 0.1)
d <- as.data.frame(matrix(NA, length(beta_a),length(rho_a)))
res<-array(unlist(d), dim=c(length(beta_a), length(rho_a)))

## set initial values ## note these need to change for meta-population analyses below #####

initials <- c(S = 1000, I = 10, R = 0)
end.time <- 20 * 365 # 20 years

############## MODEL 1 SIR WITH MORTALITY IN I ######
# unit == day (per day)
model1 =
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 6)
        change <- matrix(0, nrow = 6, ncol = 3)
        N <- S + I + R
        tau <- 1
        rate[1] <- beta * S * I/N
        change[1, ] <- c(-1, 1, 0)
        rate[2] <- (gamma+mu)/(1-rho) * I
        change[2, ] <- c(0, -1, 0)
        rate[3] <- mu * N
        change[3, ] <- c(1, 0, 0)
        rate[4] <- mu * S
        change[4, ] <- c(-1, 0, 0)
        rate[5] <- mu * R
        change[5, ] <- c(0, 0, -1)
        rate[6] <- gamma * I
        change[6, ] <- c(0, 0, 1)
        init <- c(S = S, I = I, R = R)
        for (i in 1:6) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    S <- I <- R <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      S <- c(S, init["S"])
      I <- c(I, init["I"])
      R <- c(R, init["R"])
      init <- tmp
    }
    return(list(pars = pars, init = init2, time = time, results = data.frame(time, 
                                                                             S, I, R)))
  }

############## --> TEST MODEL 1 WITH BASIC PLOT ######

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, 
                N = sum(initials), tau = 1, rho = 0)

res2 <- model1(pars = parameters, init = initials,
               end.time = end.time)
PlotMods(res2)
min(subset(res2$results,I==0)$time)

############## MODEL 2 SIR WITH MORTALITY IN I & IMPORTATION #########
#unit == day (per day)
model2=
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 8)
        change <- matrix(0, nrow = 8, ncol = 3)
        N <- S + I + R 
        tau <- 1 
        rate[1] <- beta * S * I/N
        change[1, ] <- c(-1, 1, 0)
        rate[2] <- (gamma+mu)/(1-rho) * I 
        change[2, ] <- c(0, -1, 0)
        rate[3] <- mu * N
        change[3, ] <- c(1, 0, 0)
        rate[4] <- mu * S
        change[4, ] <- c(-1, 0, 0)
        rate[5] <- mu * R
        change[5, ] <- c(0, 0, -1)
        rate[6] <- epsilon * S
        change[6, ] <- c(-1, +1, 0)
        rate[7] <- delta 
        change[7, ] <- c(0, +1, 0)
        rate[8] <- gamma * I
        change[8, ] <- c(0, 0, 1)
        init <- c(S = S, I = I, R = R)
        for (i in 1:8) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    S <- I <- R <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      S <- c(S, init["S"])
      I <- c(I, init["I"])
      R <- c(R, init["R"])
      init <- tmp
    }
    return(list(pars = pars, init = init2, time = time, results = data.frame(time, 
                                                                             S, I, R)))
  }

############## --> TEST MODEL 2 WITH BASIC PLOT #######

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = sum(initials), 
                tau = 1, rho = 0.5,
                epsilon = 2e-5, delta = 0.01)
res <- model2(pars = parameters, init = initials,
              end.time = end.time)
PlotMods(res)

############## MODEL 3 SIR WITH MORTALITY IN R ################

model3 =
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 7)
        change <- matrix(0, nrow = 7, ncol = 3)
        N <- S + I + R
        tau <- 1
        rate[1] <- beta * S * I/N
        change[1, ] <- c(-1, 1, 0)
        rate[2] <- gamma * I
        change[2, ] <- c(0, -1, 0)
        rate[3] <- mu * N
        change[3, ] <- c(1, 0, 0)
        rate[4] <- mu * S
        change[4, ] <- c(-1, 0, 0)
        rate[5] <- mu * I
        change[5, ] <- c(0, -1, 0)
        rate[6] <- mu * R
        change[6, ] <- c(0, 0, -1)
        rate[7] <- (1-rho) * gamma * I
        change[7, ] <- c(0, 0, 1)
        init <- c(S = S, I = I, R = R)
        for (i in 1:7) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    S <- I <- R <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      S <- c(S, init["S"])
      I <- c(I, init["I"])
      R <- c(R, init["R"])
      init <- tmp
    }
    return(list(pars = pars, init = init2, time = time, results = data.frame(time, 
                                                                             S, I, R)))
  }

############## --> TEST MODEL 3 WITH BASIC PLOT #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = sum(initials), 
                tau = 1, rho = 0.5)
res <- model3(pars = parameters, init = initials,
              end.time = end.time)
PlotMods(res)

############## MODEL 4 SIR WITH IMPORTS AND SPLIT I ###############################

model4 =
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 9)
        change <- matrix(0, nrow = 9, ncol = 3)
        N <- S + I + R
        tau <- 1
        rate[1] <- beta * S * I/N
        change[1, ] <- c(-1, 1, 0)
        rate[2] <- (gamma+mu)/(1-rho/2) * I
        change[2, ] <- c(0, -1, 0)
        rate[3] <- mu * N
        change[3, ] <- c(1, 0, 0)
        rate[4] <- mu * S
        change[4, ] <- c(-1, 0, 0)
        rate[5] <- mu * I
        change[5, ] <- c(0, -1, 0)
        rate[6] <- mu * R
        change[6, ] <- c(0, 0, -1)
        rate[7] <- epsilon * S
        change[7, ] <- c(-1, +1, 0)
        rate[8] <- delta
        change[8, ] <- c(0, +1, 0)
        rate[9] <- (1-rho/2) * gamma * I
        change[9, ] <- c(0, 0, 1)
        init <- c(S = S, I = I, R = R)
        for (i in 1:9) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] <
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    S <- I <- R <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      S <- c(S, init["S"])
      I <- c(I, init["I"])
      R <- c(R, init["R"])
      init <- tmp
    }
    return(list(pars = pars, init = init2, time = time, results = data.frame(time,
                                                                             S, I, R)))
  }

############## --> TEST MODEL 4 WITH BASIC PLOT ################

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = sum(initials), tau = 1, rho = 0.5,
                epsilon = 2e-5, delta = 0.01)
res <- model4(pars = parameters, init = initials,
              end.time = end.time)
PlotMods(res)

############## MODEL 5 SIR WITH IMPORT & MORTALITY IN R ###############################

model5 =
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 9)
        change <- matrix(0, nrow = 9, ncol = 3)
        N <- S + I + R
        tau <- 1
        rate[1] <- beta * S * I/N
        change[1, ] <- c(-1, 1, 0)
        rate[2] <- rho * gamma * I
        change[2, ] <- c(0, -1, 0)
        rate[3] <- mu * N
        change[3, ] <- c(1, 0, 0)
        rate[4] <- mu * S
        change[4, ] <- c(-1, 0, 0)
        rate[5] <- mu * I
        change[5, ] <- c(0, -1, 0)
        rate[6] <- mu * R
        change[6, ] <- c(0, 0, -1)
        rate[7] <- epsilon * S
        change[7, ] <- c(-1, +1, 0)
        rate[8] <- delta 
        change[8, ] <- c(0, +1, 0)
        rate[9] <- (1-rho) * gamma * I
        change[9, ] <- c(0, -1, 1)
        init <- c(S = S, I = I, R = R)
        for (i in 1:9) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    S <- I <- R <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      S <- c(S, init["S"])
      I <- c(I, init["I"])
      R <- c(R, init["R"])
      init <- tmp
    }
    return(list(pars = pars, init = init2, time = time, results = data.frame(time, 
                                                                             S, I, R)))
  }

############## --> TEST MODEL 5 WITH BASIC PLOT ######

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = sum(initials), 
                tau = 1, rho = 0.5,
                epsilon = 2e-5, delta = 0.01)
res <- model5(pars = parameters, init = initials,
              end.time = end.time)
PlotMods(res)

############## MODEL 6 SIR WITH MORTALITY IN I & 2 METAPOPULATIONS ####

model6=
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 28)
        change <- matrix(0, nrow = 28, ncol = 6)
        N <- S + I + R
        Na <- Sa + Ia + Ra
        tau <- 1
        
        rate[1] <- beta * S * I/N
        change[1, ] <- c(-1, 1, 0, 0, 0, 0)
        
        rate[2] <- (mu)* I
        change[2, ] <- c(0, -1, 0, 0, 0, 0)
        
        rate[3] <- mu * N
        change[3, ] <- c(1, 0, 0, 0, 0, 0)
        
        rate[4] <- mu * S
        change[4, ] <- c(-1, 0, 0, 0, 0, 0)
        
        rate[5] <- rho * gamma * I
        change[5, ] <- c(0, -1, 0, 0, 0, 0)
        
        rate[6] <- mu * R
        change[6, ] <- c(0, 0, -1, 0, 0, 0)
        
        rate[7] <- epsilon * S
        change[7, ] <- c(-1, +1, 0, 0, 0, 0)
        
        rate[8] <- delta * Sa
        change[8, ] <- c(1, 0, 0, 0, 0, 0)
        
        rate[9] <- delta * Ia
        change[9, ] <- c(0, 1, 0, 0, 0, 0)
        
        rate[10] <- delta * Ra
        change[10, ] <- c(0, 0, 1, 0, 0, 0)
        
        rate[11] <- delta * S
        change[11, ] <- c(-1, 0, 0, 0, 0, 0)
        
        rate[12] <- delta * I
        change[12, ] <- c(0, -1, 0, 0, 0, 0)
        
        rate[13] <- delta * R
        change[13, ] <- c(0, 0, -1, 0, 0, 0)
        
        rate[14] <- beta * Sa * Ia/Na
        change[14, ] <- c(0, 0, 0, -1, 1, 0)
        
        rate[15] <- (mu)* Ia
        change[15, ] <- c(0, 0, 0, 0, -1, 0)
        
        rate[16] <- mu * Na
        change[16, ] <- c(0, 0, 0, 1, 0, 0)
        
        rate[17] <- mu * Sa
        change[17, ] <- c(0, 0, 0, -1, 0, 0)
        
        rate[18] <- rho * gamma * Ia
        change[18, ] <- c(0, 0, 0, 0, -1, 0)
        
        rate[19] <- mu * Ra
        change[19, ] <- c(0, 0, 0, 0, 0, -1)
        
        rate[20] <- epsilon * Sa
        change[20, ] <- c(0, 0, 0, -1, +1, 0)
        
        rate[21] <- delta * Sa
        change[21, ] <- c(0, 0, 0, -1, 0, 0)
        
        rate[22] <- delta * Ia
        change[22, ] <- c(0, 0, 0, 0, -1, 0)
        
        rate[23] <- delta * Ra
        change[23, ] <- c(0, 0, 0, 0, 0, -1)
        
        rate[24] <- delta * S
        change[24, ] <- c(0, 0, 0, +1, 0, 0)
        
        rate[25] <- delta * I
        change[25, ] <- c(0, 0, 0, 0, +1, 0)
        
        rate[26] <- delta * R
        change[26, ] <- c(0, 0, 0, 0, 0, +1)
        
        
        rate[27] <- (1-rho) * gamma * I
        change[27, ] <- c(0, -1, 1, 0, 0, 0)
        
        rate[28] <- (1-rho) * gamma * Ia
        change[28, ] <- c(0, 0, 0, 0, -1, 1)
        
        
        init <- c(S = S, I = I, R = R, Sa = Sa, Ia = Ia, Ra = Ra)
        for (i in 1:28) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    S <- I <- R <- Sa <- Ia <- Ra <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      S <- c(S, init["S"])
      I <- c(I, init["I"])
      R <- c(R, init["R"])
      Sa <- c(Sa, init["Sa"])
      Ia <- c(Ia, init["Ia"])
      Ra <- c(Ra, init["Ra"])
      init <- tmp
    }
    return(list(pars = pars, init = init2, time = time, results = data.frame(time, 
                                                                             S, I, R, Sa, Ia, Ra)))
  }

############## --> TEST MODEL 6 WITH BASIC PLOT ######
initials <- c(S = 10000, I = 1, R = 0, Sa = 10000, Ia = 0, Ra = 0)

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, 
                N = sum(initials)/2, 
                Na =sum(initials)/2, tau = 1, rho = 0.5,
                epsilon = 2e-5, delta = 0.001)
res <- model6(pars = parameters, init = initials,
              end.time = end.time)
PlotMods(res)
df_meta_sir<-res$results

############## MODEL 7 SEIR WITH MORTALITY IN I ######

model7 =
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 9)
        change <- matrix(0, nrow = 9, ncol = 4)
        N <- S + E + I + R
        tau <- 1
        rate[1] <- beta * S * I/N
        change[1, ] <- c(-1, 1, 0, 0)
        
        rate[2] <- phi * E
        change[2, ] <- c(0, -1, 1, 0)
        
        rate[3] <- mu * N
        change[3, ] <- c(1, 0, 0, 0)
        
        rate[4] <- mu * I
        change[4, ] <- c(0, 0, -1, 0)
        
        rate[5] <- (1-rho) * gamma * I
        change[5, ] <- c(0, 0, -1, 1)
        
        rate[6] <- mu * S
        change[6, ] <- c(-1, 0, 0, 0)
        
        rate[7] <- mu * E
        change[7, ] <- c(0, -1, 0, 0)
        
        rate[8] <- mu * R
        change[8, ] <- c(0, 0, 0, -1)
        
        rate[9] <- rho * gamma * I
        change[9, ] <- c(0, 0, -1, 0)
        
        init <- c(S = S, E = E, I = I, R = R)
        
        for (i in 1:9) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    S <- E <- I <- R <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      S <- c(S, init["S"])
      E <- c(E, init["E"])
      I <- c(I, init["I"])
      R <- c(R, init["R"])
      init <- tmp
    }
    return(list(pars = pars, init = init2, time = time, results = data.frame(time, 
                                                                             S, E, I, R)))
  }

############## --> TEST MODEL 7 WITH BASIC PLOT ######

initials <- c(S = 10000, E = 0, I = 1, R = 0)
parameters <- c(beta = 2 / 10, phi = 1/9, gamma = 1 / 10, mu = 5e-4, N = sum(initials), tau = 1, rho = 0.5)

res <- model7(pars = parameters, init = initials,
              end.time = end.time)
PlotMods(res)
min(subset(res$results,I==0)$time)


## SINGLE RUNS - ALL MODELS ################################

############## SINGLE RUN MODEL 1 - NO MORTALITY #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0)
initials <- c(S = 10000, I = 1, R = 0)
res <- model1(pars = parameters, init = initials,
              end.time = end.time)

min(subset(res$results,I==0)$time)
df_model1_res1=res$results

############## SINGLE RUN MODEL 1 - MORTALITY #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5)
res <- model1(pars = parameters, init = initials,
              end.time = end.time)
min(subset(res$results,I==0)$time)
df_model1_res2=res$results

############## SINGLE RUN MODEL 1 - MORTALITY & DEMO CHANGE #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-3, N = 10000, tau = 1, rho = 0.5,
                epsilon = 2e-5, delta = 0.01)
res <- model1(pars = parameters, init = initials,
              end.time = end.time)
min(subset(res$results,I==0)$time)
df_model1_res3=res$results

############## SINGLE RUN MODEL 2 - NO MORTALITY #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0,
                epsilon = 2e-5, delta = 0.01)
res <- model2(pars = parameters, init = initials,
              end.time = end.time)
df_model2_res1=res$results

############## SINGLE RUN MODEL 2 - MORTALITY #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5,
                epsilon = 2e-5, delta = 0.01)
res <- model2(pars = parameters, init = initials,
              end.time = end.time)
df_model2_res2=res$results

############## SINGLE RUN MODEL 3 - MORTALITY IN R #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5)
res <- model3(pars = parameters, init = initials,
              end.time = end.time)
df_model3_res1=res$results

############## SINGLE RUN MODEL 4 - WITH IMPORTS AND SPLIT I #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5,
                epsilon = 2e-5, delta = 0.01)
res <- model4(pars = parameters, init = initials,
              end.time = end.time)
df_model4_res1=res$results

############## SINGLE RUN MODEL 5 - IMPORTS AND MORT IN R #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5,
                epsilon = 2e-5, delta = 0.01)
res <- model5(pars = parameters, init = initials,
              end.time = end.time)
df_model5_res1=res$results

############## SINGLE RUN MODEL 7 - No MORTALITY #####

initials_seir <- c(S = 10000, E = 0, I = 1, R = 0)
parameters <- c(beta = 2 / 10, phi = 1/9, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0)
res <- model7(pars = parameters, init = initials_seir,
              end.time = end.time)
df_model7_res1=res$results

############## SINGLE RUN MODEL 7 - MORTALITY #####

initials_seir <- c(S = 10000, E = 0, I = 1, R = 0)
parameters <- c(beta = 2 / 10, phi = 1/9, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5)
res <- model7(pars = parameters, init = initials_seir,
              end.time = end.time)
df_model7_res2=res$results

## MULTIPLE RUNS - ALL MODELS ################################

#### set initial values - change for meta populaions
initials <- c(S = 10000, I = 1, R = 0)
end.time <- 20 * 365
n_rep <- 100
############## MX RUNS MODEL 1 - NO MORTALITY #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0)
sim_run_m1_1<-replicate(n_rep,(model1(pars = parameters, init = initials,
                                      end.time = end.time)))

############## MX RUNS MODEL 1 - MORTALITY #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5)
sim_run_m1_2<-replicate(n_rep,(model1(pars = parameters, init = initials,
                                      end.time = end.time)))

############## MX RUNS MODEL 1 - MORTALITY & DEMO CHANGE #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-3, N = 10000, tau = 1, rho = 0.5)
sim_run_m1_3<-replicate(n_rep,(model1(pars = parameters, init = initials,
                                      end.time = end.time)))

############## MX RUNS MODEL 2 - NO MORTALITY #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0,
                epsilon = 2e-5, delta = 0.01)
sim_run_m2_1<-replicate(n_rep,(model2(pars = parameters, init = initials,
                                      end.time = end.time)))

############## MX RUNS MODEL 2 - MORTALITY #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5,
                epsilon = 2e-5, delta = 0.01)
sim_run_m2_2<-replicate(n_rep,(model2(pars = parameters, init = initials,
                                      end.time = end.time)))

############## MX RUNS MODEL 3 - MORTALITY IN R #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5)
sim_run_m3<-replicate(n_rep,(model3(pars = parameters, init = initials,
                                    end.time = end.time)))

############## MX RUNS MODEL 4 - MORTALITY in I split #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5,
                epsilon = 2e-5, delta = 0.01)
sim_run_m4<-replicate(n_rep,(model4(pars = parameters, init = initials,
                                    end.time = end.time)))

############## MX RUNS MODEL 5 - IMPORTS AND MORT IN R #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5,
                epsilon = 2e-5, delta = 0.01)
sim_run_m5<-replicate(n_rep,(model5(pars = parameters, init = initials,
                                    end.time = end.time)))

############## MX RUNS MODEL 6 - META #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000,
                tau = 1, rho = 0.5, epsilon = 0,
                delta = 0.01)
initials <- c(S = 9999, I = 1, R = 0,Sa = 10000, Ia = 0, Ra = 0)
end.time <- 20 * 365
n_rep <- 100
sim_run_m6<-replicate(n_rep,(model6(pars = parameters, init = initials,
                                    end.time = end.time)))
############## MX RUNS MODEL 6 - META * 2 #####

parameters <- c(beta = 2 / 10, gamma = 1 / 10, mu = 5e-4, N = 10000,
                tau = 1, rho = 0.5, epsilon = 2e-5,
                delta = 0.01)
initials <- c(S = 9999, I = 1, R = 0,Sa = 10000, Ia = 0, Ra = 0)
end.time <- 20 * 365
n_rep <- 100
sim_run_m6v<-replicate(n_rep,(model6(pars = parameters, init = initials,
                                     end.time = end.time)))


############## MX RUNS MODEL 7 - NO MORTALITY #####
initials <- c(S = 10000, E = 0, I = 1, R = 0)
parameters <- c(beta = 2 / 10, phi = 1/9, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0)
sim_run_m7_1<-replicate(n_rep,(model7(pars = parameters, init = initials,
                                      end.time = end.time)))

############## MX RUNS MODEL 7 - MORTALITY #####

parameters <- c(beta = 2 / 10, phi = 1/9, gamma = 1 / 10, mu = 5e-4, N = 10000, tau = 1, rho = 0.5)
sim_run_m7_2<-replicate(n_rep,(model7(pars = parameters, init = initials,
                                      end.time = end.time)))

## SET UP FUNCTIONS FOR MX METRICS ################################

############## FUNCTION  1 COUNT EXTINCTIONS ########

# minimum time that I == 0 
my_min_ext<-function(x=get_time$results, y=get_time$results$I,...){
  res_no<-vector()
  if(!is.na(min(subset(x,y==0)$time))){ #if there is no missing value
    res_no = min(subset(x,y==0)$time)   # res_no will collect the minimum time (e.g. unit = day) that I go to 0
  } else{
    res_no = end.time                   # else, res_no = the end.time
  }
  res_no
}

#N
my_min_ext_N<-function(x=get_time$results, y=get_time$results$N,...){
  res_no<-vector()
  if(!is.na(min(subset(x,y==0)))){
    res_no = min(subset(x,y==0)$time)
  } else{
    res_no = end.time
  }
  res_no
}

############## FUNCTION  2 COUNT EXTINCTIONS - NO INF #########

#at time [i], count I == 0 as 1, count I >=1 as 0

#1 = no infectious animal in the herd (count the event of disease extinction in the herd)
#0 = more than one infectious animal in the herd

my_imp_ext<-function(x=time,y=I,...){
  res_no<-vector()
  
  for (i in 2:length(x)){
    if(is.na(y[i]))
    {res_no[i-1] = NA} #if there is missing value == NA
    
    else
      #else follow this line
      
      if(y[i-1]>0 & y[i]==0 & !is.na(y[i])){ 
        res_no[i-1] = 1
      } else{
        res_no[i-1] = 0
      }
  }
  sum(res_no, na.rm = T)
}

my_imp_ext_N<-function(x=time,y=N,...){
  res_no<-vector()
  for (i in 2:length(x)){
    if(is.na(y[i]))
    {res_no[i-1] = NA}
    else
      if(y[i-1]>0 & y[i]==0 & !is.na(y[i])){
        res_no[i-1] = 1
      } else{
        res_no[i-1] = 0
      }
  }
  sum(res_no, na.rm = T)
}

############## FUNCTION  3 COUNT PERSISTENCE TIMES ########

my_imp_ext_na<-function(x=time,y=I,...){
  
  res_no<-vector()
  
  for (i in 2:length(x)) {
    if(is.na(y[i]))
    {res_no[i-1] = NA}
    else
      if(y[i-1]>0 & y[i]==0 & !is.na(y[i])){
        res_no[i-1] = 1
      } else{
        res_no[i-1] = 0
      }
  }
  length(res_no[!is.na(res_no)])
}

my_imp_ext_na_N<-function(x=time,y=N,...){
  res_no<-vector()
  for (i in 2:length(x)){
    if(is.na(y[i]))
    {res_no[i-1] = NA}
    else
      if(y[i-1]>0 & y[i]==0 & !is.na(y[i])){
        res_no[i-1] = 1
      } else{
        res_no[i-1] = 0
      }
  }
  length(res_no[!is.na(res_no)])
}


############## FUNCTION  4 COUNTS EXTINCTIONS METAPOP MODEL ##############

my_min_meta_ext<-function(x=get_time$results, y=get_time$results$I, y1=get_time$results$Ia,...){
  res_no<-vector()
  if(!is.na(min(subset(x,y==0 & y1==0)$time))){
    res_no = min(subset(x,y==0 & y1==0)$time)
  } else{
    res_no = end.time
  }
  res_no
}


## SET UP FUNCTIONS FOR MODEL RUNS LOOP THROUGH PARS ################################
#### set initial values - change for populations

#gaur poplation
N = 300 
#ratio c:sa:a
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

end.time <- 10 * 365
n_rep <- 10

############## FUNCTION  6 MODEL PAR LOOPS & MIN TIME TO EXTINCTION (no epsilon) #######

### model 2 Anthrax SI ####
my_fun_model2_ant<-function() {
  for (k in 1:length(beta)){
    for (i in 1:length(rho)) {
      parameters <- c(beta_a = beta[k], beta_sa = beta[k] , beta_c = beta[k],
                      gamma_c = (1/(1/24))/365,
                      gamma_sa = (1/(1/24))/365,
                      gamma_a = (1/(1/24))/365,
                      rho_c= rho[i],
                      rho_sa= rho[i],
                      rho_a=rho[i],
                      N = sum(initials),
                      tau=1,
                      mu_b = 0.34/365, 
                      mu_c = 0.27/365, 
                      mu_sa = 0.15/365,
                      mu_a = 0.165/365,
                      delta_c = 1/365,
                      delta_sa = 1/(3*365))
      
      initials <- c(Sc = c, Ic = 0, Ssa = sa, Isa = 0, Sa = (a-1), Ia = 1 )
      
      get_time <- model2(pars = parameters, init = initials,
                         end.time = end.time)
      res_min[k,i]<- my_min_ext_N(x=get_time$results,y=get_time$results$N)
      res_num_ext[k,i]<- my_imp_ext_N(x=get_time$time,y=get_time$results$N)
      res_time_inf[k,i]<- my_imp_ext_na_N(x=get_time$time,y=get_time$results$N)
    }}
  d <- list(res_min,res_num_ext,res_time_inf)
  d
  }
  
############## FUNCTION  7 MODEL PAR LOOPS & COUNT EXTINCTIONS (with epsilon) #######

##### model 2 Anthrax SI #####
my_fun_model2_ant_ep<-function() {
  for (k in 1:length(beta)){
  for (i in 1:length(rho)) {
    parameters <- c(beta_a = beta[k], 
                    beta_sa = beta[k], 
                    beta_c = beta[k],
                    gamma_c = (1/(1/24))/365,
                    gamma_sa = (1/(1/24))/365,
                    gamma_a = (1/(1/24))/365,
                    rho_c= rho[i],
                    rho_sa= rho[i],
                    rho_a= rho[i],
                    epsilon = 2e-5,
                    N = sum(initials),
                    tau=1,
                    mu_b = 0.34/365, 
                    mu_c = 0.27/365, 
                    mu_sa = 0.15/365,
                    mu_a = 0.165/365,
                    delta_c = 1/365,
                    delta_sa = 1/(3*365))
    
    initials <- c(Sc = c, Ic = 0, Ssa = sa, Isa = 0, Sa = (a-1), Ia = 1 )
    get_time <- model2(pars = parameters, init = initials,
                       end.time = end.time)
    res_min[k,i]<- my_min_ext_N(x=get_time$results,y=get_time$results$N)
    res_num_ext[k,i]<- my_imp_ext_N(x=get_time$time,y=get_time$results$N)
    res_time_inf[k,i]<- my_imp_ext_na_N(x=get_time$time,y=get_time$results$N)
  }}
  d <- list(res_min,res_num_ext,res_time_inf)
  d
}

## PLOT DATA PREPARATION #################################

out_put_fun_min<-function(x, par1, par2, ...){ # x = out, par1, par2 = beta_a, rho_a, ...
  output <- array(unlist(x), dim = c(length(par1),length(par2),n_rep*3))
  output_min<-output[,,seq(from=1,to=n_rep*3,by=3)]
  output_min
}

out_put_fun_num<-function(x, par1, par2, ...){ # x = out, par1, par2 = beta_a, rho_a, ...
  output <- array(unlist(x), dim = c(length(par1),length(par2),n_rep*3))
  output_num<-output[,,seq(from=2,to=n_rep*3,by=3)]
  output_num
}

out_put_fun_time<-function(x, par1, par2, ...){ # x = out, par1, par2 = beta_a, rho_a, ...
  output <- array(unlist(x), dim = c(length(par1),length(par2),n_rep*3))
  output_time<-output[,,seq(from=3,to=n_rep*3,by=3)]/end.time
  output_time
}

############## FUNCTION 15 PREP OUTPUT FOR PLOTTING TIME TO EXTINCTION or NUMBER OF EXTINCTIONS #############

my_plot_min<-function(x, par1, par2, par1_n, par2_n){ 
  # x = data run, e.g. big_run, par1, par2 are parameters varied, e.g. beta_a, rho_a, par1_n, par2_n are names, e.g. 'beta', 'rho', op is output, either time or extinctions
  plot_res<-apply(x, c(1,2), mean, na.rm = T)
  row.names(plot_res)<-par1; colnames(plot_res)<-par2
  df <- melt(plot_res)
  colnames(df)<-c(par1_n,par2_n,'duration')
  df
}

############## FUNCTION 16 PREP OUTPUT FOR PLOTTING PROPORTION OF OUTBREAKS PERSISTING ########

my_plot_num<-function(x, par1, par2, par1_n, par2_n){ # x = data run, e.g. big_run, par1, par2 are parameters varied, e.g. beta_a, rho_a
  x[!is.finite(x)] <- 0
  plot_res<-apply(x, c(1,2), sum, na.rm = T)/n_rep
  row.names(plot_res)<-par1; colnames(plot_res)<-par2
  df <- melt(plot_res)
  colnames(df)<-c(par1_n,par2_n,'outbreaks')
  df
}

############## FUNCTION 17 PREP OUTPUT FOR PLOTTING TIME TO EXTINCTION WITH NO INFINITE FROM PERSISTENCE ######

my_plot_time<-function(x, par1, par2, par1_n, par2_n){ # x = data run, e.g. big_run, par1, par2 are parameters varied, e.g. beta_a, rho_a
  x[!is.finite(x)] <- end.time
  plot_res<-apply(x, c(1,2), mean, na.rm = T)
  row.names(plot_res)<-par1; colnames(plot_res)<-par2
  df <- melt(plot_res)
  colnames(df)<-c(par1_n,par2_n,'time')
  df
}

############## FUNCTION 18 PREP OUTPUT FOR PLOTTING MX SIMULATIONS - SINGLE POP ####

single_pop_sim_prep <- function(x, n_rep, end.time){ # x = simulation of model, e.g. sim_run_m1
  mat = matrix(NA, nrow=n_rep, ncol = end.time+1)
  for (i in 1:n_rep){
    mat[i,]<-x[,i]$results$N
  }
  colnames(mat) = paste("time", seq(from=1,to=end.time+1,by=1), sep="")
  rownames(mat) = paste("run", seq(n_rep), sep="")
  dat = as.data.frame(mat)
  dat$run = rownames(dat)
  mdat = melt(dat, id.vars="run")
  mdat$time = as.numeric(gsub("time", "", mdat$variable))
  mdat
}

############## PARAMETER RANGE #######################

beta <- seq(from = 5e-5,to = 0.01, length.out=10)
rho <- seq(from = 1, to = 1, length.out=10)
#delta <- seq(from = 0.001,to = 0.01, by = 0.001)

d <- as.data.frame(matrix(NA, length(beta),length(rho)))
res_min<-array(unlist(d), dim=c(length(beta), length(rho)))
res_num_ext<-array(unlist(d), dim=c(length(beta), length(rho)))
res_time_inf<-array(unlist(d), dim=c(length(beta), length(rho)))

n_rep = 10
end.time = 5 * 365

############## RUN MODEL Anthrax #####
#no epsilon
big_model2_ant<-replicate(n_rep,my_fun_model2_ant())
#with epsilon
big_model2_ant_ep<-replicate(n_rep,my_fun_model2_ant_ep())

############## RUN MODEL 2 #####

big_run_model2<-replicate(n_rep,my_fun_model2())

############## RUN MODEL 3 ##############

big_run_model3<-replicate(n_rep,my_fun_model3())

############## RUN MODEL 4 ###########

big_run_model4<-replicate(n_rep,my_fun_model4())

############## RUN MODEL 5 ###########

big_run_model5<-replicate(n_rep,my_fun_model5())

########## for METAPOPULATION RUN ONLY IF LENGTH DELTA ! == BETA ######

d <- as.data.frame(matrix(NA, length(delta_a),length(rho_a)))
res_min<-array(unlist(d), dim=c(length(delta_a), length(rho_a)))
res_num_ext<-array(unlist(d), dim=c(length(delta_a), length(rho_a)))
res_time_inf<-array(unlist(d), dim=c(length(delta_a), length(rho_a)))

############## RUN MODEL 6 ###########

big_run_meta<-replicate(n_rep,my_fun_model6())

############## RUN MODEL 6 MIG ###########

big_run_model6_mig<-replicate(n_rep,my_fun_model6_mig())

## MODEL OUTPUTS PREPARATION #######################

############## MODELS 1-6 OUTPUT PREPARATION #############################

mod_res<-list(big_run_model1,
              big_run_model2,
              big_run_model3,
              big_run_model4,
              big_run_model5,
              big_run_meta)

for (i in 1:length(mod_res)){
  assign(paste0("Res_min_", i), out_put_fun_min(x=mod_res[[i]], par1 = beta, par2 = rho))
}

for (i in 1:length(mod_res)){
  assign(paste0("Res_num_", i), out_put_fun_num(x=mod_res[[i]],par1 = beta ,par2 = rho))
}

for (i in 1:length(mod_res)){
  assign(paste0("Res_time_", i), out_put_fun_time(x=mod_res[[i]],par1 = beta,par2 = rho))
}

plot_res_min<-list(Res_min_1,
                   Res_min_2)

plot_res_min <- lapply(plot_res_min,function(x) replace(x,is.infinite(x),end.time))

plot_res_num<-list(Res_num_1,
                   Res_num_2)

plot_res_time<-list(Res_time_1,
                    Res_time_2)

for (i in 1:length(plot_res_min)){
  assign(paste0("df_min_", i), my_plot_min(x=plot_res_min[[i]],par1 = beta_a, par2 = rho_a, par1_n = 'beta', par2_n = 'rho'))
}

for (i in 1:length(plot_res_num)){
  assign(paste0("df_num_", i), my_plot_num(x=plot_res_num[[i]],par1 = beta_a,par2 = rho_a, par1_n = 'beta', par2_n = 'rho'))
}

for (i in 1:length(plot_res_time)){
  assign(paste0("df_time_", i), my_plot_time(x=plot_res_time[[i]],par1 = beta_a,par2 = rho_a, par1_n = 'beta', par2_n = 'rho'))
}

n = length(plot_res_time)
for(i in 1:n){
  t1 <- do.call(cbind, mget(paste0("df_min_", 1:n) ) )
}
for(i in 1:n){
  t2 <- do.call(cbind, mget(paste0("df_num_", 1:n) ) )
}
for(i in 1:n){
  t3 <- do.call(cbind, mget(paste0("df_time_", 1:n) ) )
}

res<-cbind(t1,t2,t3)

colnames(res) <- colnames(res) %>% str_replace(".*.beta", "beta")
colnames(res) <- colnames(res) %>% str_replace(".*.rho", "rho")

res_all = melt(res, id.vars=c("beta",'rho'))

## PLOT MODEL OUTPUTS TILES #############################

for (k in unique(res_all$variable)){
  subdata <- subset(res_all, variable == k)
  if (str_detect(subdata$variable, "duration")==T){
    pname <- paste0("duration365-",k)
    p<-(ggplot(subdata, aes(x = beta, y = rho, fill = value))+
          geom_tile()+
          scale_fill_gradientn(colors=colorRampPalette(c("whitesmoke","royalblue","seagreen","orange","red","brown"))(500),
                               name="Average\noutbreak\nduration", na.value = "grey", limits = c(0,365)) +
          labs(x = expression(beta),y=expression(rho)) +
          theme_bw())
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)
    pname <- paste0("duration-",k)
    p<-(ggplot(subdata, aes(x = beta, y = rho, fill = value))+
          geom_tile()+
          scale_fill_gradientn(colors=colorRampPalette(c("whitesmoke","royalblue","seagreen","orange","red","brown"))(500),
                               name="Average\noutbreak\nduration", na.value = "grey") +
          labs(x = expression(beta),y=expression(rho)) +
          theme_bw())
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)}
  else if
  (str_detect(subdata$variable, "outbreaks")==T){
    pname <- paste0("extinctions100-",k)
    p<-(ggplot(subdata, aes(x = beta, y = rho, fill = value))+
          geom_tile()+
          scale_fill_gradientn(colors=colorRampPalette(c("brown","red","orange","seagreen","royalblue","whitesmoke"))(500),
                               name="Average\nnumber\nextinctions", limits = c(0,100)) +
          labs(x = expression(beta),y=expression(rho)) +
          theme_bw())
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)
    pname <- paste0("extinctions-",k)
    p<-(ggplot(subdata, aes(x = beta, y = rho, fill = value))+
          geom_tile()+
          scale_fill_gradientn(colors=colorRampPalette(c("brown","red","orange","seagreen","royalblue","whitesmoke"))(500),
                               name="Average\nnumber\nextinctions") +
          labs(x = expression(beta),y=expression(rho)) +
          theme_bw())
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)}
  else
  {
    pname <- paste0("time-",k)
    p<-(ggplot(subdata, aes(x = beta, y = rho, fill = value))+
          geom_tile()+
          scale_fill_gradientn(colors=colorRampPalette(c("brown","red","orange","seagreen","royalblue","whitesmoke"))(500),
                               name="Average\ntime", limits = c(0,1)) +
          labs(x = expression(beta),y=expression(rho)) +
          theme_bw())
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)}
}


## PLOT MODEL OUTPUTS TIME SERIES #############################

for (k in unique(res_alld$variable)){
  subdata <- subset(res_alld, variable == k)
  if (str_detect(subdata$variable, "duration")==T){
    pname <- paste0("metap_duration365_",k)
    p<-(ggplot(subdata, aes(x = delta, y = rho, fill = value))+
          geom_tile()+
          scale_fill_gradientn(colors=colorRampPalette(c("whitesmoke","royalblue","seagreen","orange","red","brown"))(500),
                               name="Average\noutbreak\nduration", na.value = "grey", limits = c(0,365)) +
          labs(x = expression(delta),y=expression(rho)) +
          theme_bw())
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)
    pname <- paste0("metap_duration_",k)
    p<-(ggplot(subdata, aes(x = delta, y = rho, fill = value))+
          geom_tile()+
          scale_fill_gradientn(colors=colorRampPalette(c("whitesmoke","royalblue","seagreen","orange","red","brown"))(500),
                               name="Average\noutbreak\nduration", na.value = "grey") +
          labs(x = expression(delta),y=expression(rho)) +
          theme_bw())
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)}
  else if
  (str_detect(subdata$variable, "outbreaks")==T){
    pname <- paste0("metap_extinction100_",k)
    p<-(ggplot(subdata, aes(x = delta, y = rho, fill = value))+
          geom_tile()+
          scale_fill_gradientn(colors=colorRampPalette(c("brown","red","orange","seagreen","royalblue","whitesmoke"))(500),
                               name="Average\nnumber\nextinctions", limits = c(0,100)) +
          labs(x = expression(delta),y=expression(rho)) +
          theme_bw())
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)
    pname <- paste0("metap_extinction_",k)
    p<-(ggplot(subdata, aes(x = delta, y = rho, fill = value))+
          geom_tile()+
          scale_fill_gradientn(colors=colorRampPalette(c("brown","red","orange","seagreen","royalblue","whitesmoke"))(500),
                               name="Average\nnumber\nextinctions") +
          labs(x = expression(delta),y=expression(rho)) +
          theme_bw())
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)}
  else
  {
    pname <- paste0("metap_time_",k)
    p<-(ggplot(subdata, aes(x = delta, y = rho, fill = value))+
          geom_tile()+
          scale_fill_gradientn(colors=colorRampPalette(c("brown","red","orange","seagreen","royalblue","whitesmoke"))(500),
                               name="Average\ntime", limits = c(0,1)) +
          labs(x = expression(delta),y=expression(rho)) +
          theme_bw())
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)}
}

## SINGLE POPULATION PRINT PREPS #######

res_mx<-list(sim_run_m1_1,
             sim_run_m1_2,
             sim_run_m1_3,
             sim_run_m2_1,
             sim_run_m2_2,
             sim_run_m3,
             sim_run_m4,
             sim_run_m5,
             sim_run_m7_1,
             sim_run_m7_2)

for (i in 1:length(res_mx)){
  assign(paste0("df_mx", i), single_pop_sim_prep(x=res_mx[[i]],n_rep = n_rep, end.time = end.time))
}

res_mx_p<-rbind(df_mx1,
                df_mx2,
                df_mx3,
                df_mx4,
                df_mx5,
                df_mx6,
                df_mx7,
                df_mx8,
                df_mx9,
                df_mx10)

res_mx_p$model<-c(rep('1',dim(df_mx1)[1]),
                  rep('2',dim(df_mx2)[1]),
                  rep('3',dim(df_mx3)[1]),
                  rep('4',dim(df_mx4)[1]),
                  rep('5',dim(df_mx5)[1]),
                  rep('6',dim(df_mx6)[1]),
                  rep('7',dim(df_mx7)[1]),
                  rep('8',dim(df_mx8)[1]),
                  rep('9',dim(df_mx9)[1]),
                  rep('10',dim(df_mx10)[1]))

## SINGLE POPULATION PLOTS #######

for (i in unique(res_mx_p$model)){
  subdata <- subset(res_mx_p, model == i)
  pdf(paste("plot_ts_mx", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot(subdata, aes(x=time, y=value, group=run)) +
          theme_bw() +
          theme(panel.grid=element_blank()) +
          geom_line(size=0.2, alpha=0.15)+
          ylab('Numbers') + xlab('time')+
          stat_summary(aes(group = 1), fun.y=mean, geom="line", colour="black",size = 1.1))
  dev.off()
}


############### PLOT ALL TS SINGLE RUNS ##############

res_p_ts<-list(sim_run_m1_1,
               sim_run_m1_2,
               sim_run_m1_3,
               sim_run_m2_1,
               sim_run_m2_2,
               sim_run_m3,
               sim_run_m4,
               sim_run_m5,
               sim_run_m6,
               sim_run_m6v,
               sim_run_m7_1,
               sim_run_m7_2)

# Make plots. Single runs
for (i in 1:length(res_p_ts)) {
  pdf(paste("plotts", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot() + 
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$S), color = "black", size =1.2) +
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$I), color = "red", size =1.2) +
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$R), color = "seagreen", size =1.2, linetype = "dotted") +
          xlab('Time') +
          ylab('Infection State Numbers'))
  dev.off()
}

############### PLOT ALL TS SINGLE RUNS I ONLY ##############

# Make plots. Single runs
for (i in 1:length(res_p_ts)) {
  pdf(paste("plotts_i", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot() + 
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$I), color = "red", size =1.2) +
          xlab('Time') +
          ylab('Infection Numbers'))
  dev.off()
}

# Make plots. Single runs - SCALED
for (i in 1:length(res_p_ts)) {
  pdf(paste("plotts_i_scale", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot() + 
          geom_line(data = res_p_ts[[i]][[4]], 
                    aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$I), color = "red", size =1.2) +
          xlab('Time') +
          ylim(0, 500) + 
          ylab('Infection Numbers'))
  dev.off()
}

## plot extinction times

df.agg <- aggregate(time ~ run + value + model, res_mx_p, min)
df.ag<-(df.agg[df.agg$value==0,c('model','time')])

neworder <- c("1","2","3","4","5","9","10","6","7","8")
library(plyr)  ## or dplyr (transform -> mutate)
df.ag <- arrange(transform(df.ag,
                           model=factor(model,levels=neworder)),model)
labs <- c('1' = "SIR",
          '2' = "SI*R",
          '3' = "SI*R+",
          '4' = "S[I]R",
          '5' = "S[I]*R",
          '6' = "SIR*",
          '7' = 'S[I]*R',
          '8' = 'S[I]R*',
          '9' = 'SEIR',
          '10' = 'SEI*R')
p<-ggplot(df.ag, aes(x=time))+
  geom_histogram(color="black", fill="grey")+
  facet_wrap(model~., ncol = 5,labeller = labeller(model = labs))+
  scale_x_continuous(breaks = c(0, 800, 1600), labels = c("0", "800", "1600"))
pdf("extinctions.pdf", width = 5, height = 3)
p
dev.off()

p_yr<-ggplot(df.ag, aes(x=time))+
  geom_histogram(color="black", fill="grey")+
  facet_wrap(model~., ncol = 5,labeller = labeller(model = labs))+
  scale_x_continuous(limits = c(0,1000),breaks = c(0, 400, 800), labels = c("0", "400", "800")) +
  ylim(0,60)
pdf("extinctions_1yr.pdf", width = 5, height = 3)
p_yr
dev.off()