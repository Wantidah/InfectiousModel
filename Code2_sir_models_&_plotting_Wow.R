############## NOTES ##################
# Infectious disease modelling code for wild Bovidae population
# - 18 Nov 2022-
# # Reference R Code: https://github.com/dtsh2/ebola_model; 
# # Reference Article:https://doi.org/10.1098/rsif.2021.0638  (Transmission models indicate Ebola virus persistence in non-human primate populations is unlikely)

# All models have 3 age classes (c = calf, a = Subadult, a = adult)
# unit == day (per day)
# This code set up as:
# 1) Set up the infectious disease model's function (7 models, details below)
# 2) Single run the model's function with parameters 
# 3) Multiple run (default = 100 times) the models with parameters 
# 4) Function for run the model with parameter ranges
# 5) the outputs will be in matrix -> extract and plot the results

rm(list=ls())
set.seed(111)
############## LOAD PACKAGES ##########
library(EpiDynamics)
library(dplyr)   
library(tidyverse)
library(reshape2) 
library(stringr)
library(emdbook)  
library(ggplot2); theme_set(theme_bw())

############## Run the model function ############### 

############## MODEL 1 - population dynamic model, no infection ######
model1 =
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 6)
        change <- matrix(0, nrow = 6, ncol = 3)
        
        N <- c + sa + a
        tau <- 1
        
        rate[1] <- mu_b * a           
        change[1, ] <- c(1, 0, 0)
        rate[2] <- mu_c * c           
        change[2, ] <- c(-1, 0, 0)
        rate[3] <- delta_c * c    
        change[3, ] <- c(-1, 1, 0)
        rate[4] <- delta_sa * sa
        change[4, ] <- c(0, -1, 1)
        rate[5] <- mu_sa * sa 
        change[5, ] <- c(0, -1, 0)
        rate[6] <-  mu_a * a
        change[6, ] <- c(0, 0, -1)
        
        init <- c(c = c, sa = sa, a = a)
        for (i in 1:6) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    a <- sa <- c <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      c <- c(c, init["c"])
      sa <- c(sa, init["sa"])
      a <- c(a, init["a"])
      init <- tmp
    }
    
    #sum population based on column name
    results<-data.frame(time, 
                        c,  sa,  a)%>% 
      dplyr::mutate(N = rowSums(across(-c(time), na.rm=TRUE)))
    
    return(list(pars = pars, init = init2, time = time, results = results))
    
  }

############## MODEL 2 SI - Anthrax #####
model2 = 
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 15)
        change <- matrix(0, nrow = 15, ncol = 6)
        
        N <- Sc+Ic +Ssa+Isa +Sa+Ia 
        tau <- 1
        
        rate[1] <- mu_b * Sa
        change[1, ] <- c(1, 0, 0, 0, 0, 0)
        rate[2] <- beta_c * Sc * (Ic+Isa+Ia)/N
        change[2, ] <- c(-1, 1, 0, 0, 0, 0)
        rate[3] <- gamma_c * Ic * rho_c
        change[3, ] <- c(0, -1, 0, 0, 0, 0)
        rate[4] <- mu_c * Sc
        change[4, ] <- c(-1, 0, 0, 0, 0, 0)
        rate[5] <- delta_c * Sc
        change[5, ] <- c(-1, 0, 1, 0, 0, 0)  
        rate[6] <- epsilon * Sc
        change[6, ] <- c(-1, 1, 0, 0, 0, 0)
        
        #saubadult
        rate[7] <- beta_sa * Ssa * (Ic+Isa+Ia)/N
        change[7, ] <- c(0, 0, -1, 1, 0, 0)    
        rate[8] <-  gamma_sa * Isa * rho_sa
        change[8, ] <- c(0, 0, 0, -1, 0, 0)    
        rate[9] <- mu_sa * Ssa
        change[9, ] <- c(0, 0, -1, 0, 0, 0)
        rate[10] <- delta_sa * Ssa
        change[10, ] <- c(0, 0, -1, 0, 1, 0)  
        rate[11] <- epsilon * Ssa
        change[11, ] <- c(0, 0, -1, 1, 0, 0)
        
        #adult
        rate[12] <- beta_a * Sa * (Ic+Isa+Ia)/N
        change[12, ] <- c(0, 0, 0, 0, -1, 1)    
        rate[13] <-  gamma_a * Ia *rho_a
        change[13, ] <- c(0, 0, 0, 0, 0, -1)    
        rate[14] <- mu_a * Sa
        change[14, ] <- c(0, 0, 0, 0, -1, 0)
        rate[15] <- epsilon * Sa
        change[15, ] <- c(0, 0, 0, 0, -1, 1)
        
        init <- c(Sc = Sc, Ic = Ic, Ssa = Ssa, Isa = Isa, Sa = Sa, Ia = Ia)
        for (i in 1:15) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] <0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    Sa <- Ia <- Ssa <- Isa <- Sc <- Ic <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      Sa <- c(Sa, init["Sa"])
      Ia <- c(Ia, init["Ia"])
      
      Ssa <- c(Ssa, init["Ssa"])
      Isa <- c(Isa, init["Isa"])
      
      Sc <- c(Sc, init["Sc"])
      Ic <- c(Ic, init["Ic"])
      
      init <- tmp
    }
    
    #sum population based on column name
    results<-data.frame(time, 
                        Sc,  Ic,  Ssa, Isa, Sa, Ia)%>% 
      dplyr::mutate(N = rowSums(across(-c(time), na.rm=TRUE)))%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))
    
    return(list(pars = pars, init = init2, time = time, results = results))
    
  }
############## MODEL 3 SEI - Bovine tuberculosis  #####

model3=
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 29)
        change <- matrix(0, nrow = 29, ncol = 9)
        
        N <- Sc+Ec+Ic +Ssa+Esa+Isa +Sa+Ea+Ia 
        tau <- 1
        
        #calf
        rate[1] <- mu_b * (Sa + Ea) 
        change[1, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[2] <- beta_c * Sc * (Ic+Isa+Ia)/N
        change[2, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0)
        rate[3] <- phi_c * Ec 
        change[3, ] <- c(0, -1, 1, 0, 0, 0, 0, 0, 0)
        rate[4] <-  rho_c * gamma_c * Ic
        change[4, ] <- c(0, 0, -1, 0, 0, 0, 0, 0, 0)
        rate[5] <-  delta_c * Sc
        change[5, ] <- c(-1, 0, 0, 1, 0, 0, 0, 0, 0)  
        rate[6] <-  delta_c * Ec
        change[6, ] <- c(0, -1, 0, 0, 1, 0, 0, 0, 0) 
        rate[7] <-  delta_c * Ic
        change[7, ] <- c(0, 0, -1, 0, 0, 1, 0, 0, 0) 
        rate[8] <-  mu_c * Sc
        change[8, ] <- c(-1, 0, 0, 0, 0, 0, 0, 0, 0)  
        rate[9] <- mu_c * Ec
        change[9, ] <- c(0, -1, 0, 0, 0, 0, 0, 0, 0)
        rate[10] <- mu_c * Ic
        change[10, ] <- c(0, 0, -1, 0, 0, 0, 0, 0, 0)
        rate[11] <- epsilon * Sc
        change[11, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0)
        
        # saubadult
        rate[12] <- beta_sa * Ssa * (Ic+Isa+Ia)/N
        change[12, ] <- c(0, 0, 0, -1, 1, 0, 0, 0, 0)
        rate[13] <- phi_sa * Esa 
        change[13, ] <- c(0, 0, 0, 0,-1, 1, 0, 0, 0)
        rate[14] <-  rho_sa * gamma_sa * Isa
        change[14, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 0)
        rate[15] <-  delta_sa * Ssa
        change[15, ] <- c(0, 0, 0, -1, 0, 0, 1, 0, 0)  
        rate[16] <-  delta_sa * Esa
        change[16, ] <- c(0, 0, 0, 0, -1, 0, 0, 1, 0) 
        rate[17] <-  delta_sa *  Isa
        change[17, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 1) 
        rate[18] <-  mu_sa * Ssa
        change[18, ] <- c(0, 0, 0, -1, 0, 0, 0, 0, 0)  
        rate[19] <- mu_sa * Esa
        change[19, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 0)
        rate[20] <- mu_sa * Isa
        change[20, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0)
        rate[21] <- epsilon * Ssa
        change[21, ] <- c(0, 0, 0, -1, 1, 0, 0, 0, 0)
        
        # adult
        rate[22] <- beta_a * Sa * (Ic+Isa+Ia)/N
        change[22, ] <- c(0, 0, 0, 0, 0, 0, -1, 1, 0)
        rate[23] <- phi_a * Ea 
        change[23, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 1)
        rate[24] <-   rho_a * gamma_a * Ia
        change[24, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[25] <-  mu_a * Sa
        change[25, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0)  
        rate[26] <- mu_a * Ea
        change[26, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0)
        rate[27] <- mu_a * Ia
        change[27, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[28] <- epsilon * Sa
        change[28, ] <- c(0, 0, 0, 0, 0, 0, -1, 1, 0)
        # birth rate from infected mother 
        rate[29] <- mu_bI * Ia
        change[29, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
        
        init <- c(Sc = Sc, Ec = Ec, Ic = Ic, Ssa = Ssa, Esa = Esa, Isa = Isa, Sa = Sa, Ea = Ea, Ia = Ia)
        for (i in 1:29) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    Sa <- Ea <- Ia <- Ssa <- Esa <- Isa <- Sc <- Ec <- Ic <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      Sa <- c(Sa, init["Sa"])
      Ea <- c(Ea, init["Ea"])
      Ia <- c(Ia, init["Ia"])
      
      Ssa <- c(Ssa, init["Ssa"])
      Esa <- c(Esa, init["Esa"])
      Isa <- c(Isa, init["Isa"])
      
      Sc <- c(Sc, init["Sc"])
      Ec <- c(Ec, init["Ec"])
      Ic <- c(Ic, init["Ic"])
      
      init <- tmp
    }
    
    #sum population based on column name
    results = data.frame(time, 
                         Sc, Ec, Ic, Ssa, Esa, Isa, Sa, Ea, Ia)%>% 
      dplyr::mutate(N = rowSums(across(-c(time), na.rm=TRUE)))%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))
    
    return(list(pars = pars, init = init2, time = time, results = results))
  }

############## MODEL 4 SIRS - Hemorrhagic septicemia #####
model4 =
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 31)
        change <- matrix(0, nrow = 31, ncol = 9)
        N <- Sc + Ic + Rc + Ssa + Isa + Rsa + Sa + Ia + Ra
        
        tau <- 1
        #calf
        rate[1] <- mu_b * (Sa + Ia + Ra)
        change[1, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[2] <- beta_c * Sc * (Ic+Isa+Ia)/N
        change[2, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0)
        rate[3] <- (1-rho_c) * gamma_c * Ic
        change[3, ] <- c(0, -1, 1, 0, 0, 0, 0, 0, 0)
        rate[4] <- rho_c *  gamma_c * Ic
        change[4, ] <- c(0, -1, 0, 0, 0, 0, 0, 0, 0)
        rate[5] <-  omega_c *  Rc
        change[5, ] <- c(1, 0, -1, 0, 0, 0, 0, 0, 0)
        rate[6] <-  delta_c * Sc
        change[6, ] <- c(-1, 0, 0, 1, 0, 0, 0, 0, 0)  
        rate[7] <-  delta_c * Ic
        change[7, ] <- c(0, -1, 0, 0, 1, 0, 0, 0, 0) 
        rate[8] <-  delta_c * Rc
        change[8, ] <- c(0, 0, -1, 0, 0, 1, 0, 0, 0) 
        rate[9] <-  mu_c * Sc
        change[9, ] <- c(-1, 0, 0, 0, 0, 0, 0, 0, 0)  
        rate[10] <- mu_c * Ic
        change[10, ] <- c(0, -1, 0, 0, 0, 0, 0, 0, 0)
        rate[11] <- mu_c * Rc
        change[11, ] <- c(0, 0, -1, 0, 0, 0, 0, 0, 0)
        rate[12] <- epsilon * Sc
        change[12, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0)
        
        #subadult
        rate[13] <- beta_sa * Ssa * (Ic+Isa+Ia)/N
        change[13, ] <- c(0, 0, 0, -1, 1, 0, 0, 0, 0)
        rate[14]<- (1-rho_sa) * gamma_sa * Isa
        change[14, ] <- c(0, 0, 0, 0, -1, 1, 0, 0, 0)
        rate[15] <-  rho_sa * gamma_sa * Isa
        change[15, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 0)
        rate[16] <-  omega_sa *  Rsa
        change[16, ] <- c(0, 0, 0, 1, 0, -1, 0, 0, 0)
        rate[17] <-  delta_sa * Ssa
        change[17, ] <- c(0, 0, 0, -1, 0, 0, 1, 0, 0)  
        rate[18] <-  delta_sa * Isa
        change[18, ] <- c(0, 0, 0, 0, -1, 0, 0, 1, 0) 
        rate[19] <-  delta_sa *  Rsa
        change[19, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 1) 
        rate[20] <-  mu_sa * Ssa
        change[20, ] <- c(0, 0, 0, -1, 0, 0, 0, 0, 0)  
        rate[21] <- mu_sa * Isa
        change[21, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 0)
        rate[22] <- mu_sa * Rsa
        change[22, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0)
        rate[23] <- epsilon * Ssa
        change[23, ] <- c(0, 0, 0, -1, 1, 0, 0, 0, 0)
        
        #adult
        rate[24] <- beta_a * Sa * (Ic+Isa+Ia)/N
        change[24, ] <- c(0, 0, 0, 0, 0, 0, -1, 1, 0)
        rate[25] <- (1-rho_a) * gamma_a * Ia
        change[25, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 1)
        rate[26] <-  rho_a * gamma_a * Ia
        change[26, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[27] <-  omega_a *  Ra
        change[27, ] <- c(0, 0, 0, 0, 0, 0, 1, 0, -1)
        rate[28] <-  mu_a * Sa
        change[28, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0)  
        rate[29] <- mu_a * Ia
        change[29, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0)
        rate[30] <- mu_a * Ra
        change[30, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[31] <- epsilon * Sa
        change[31, ] <- c(0, 0, 0, 0, 0, 0, -1, 1, 0)
        
        init <- c(Sc = Sc, Ic = Ic, Rc = Rc, 
                  Ssa = Ssa, Isa = Isa,  Rsa = Rsa, 
                  Sa = Sa, Ia = Ia, Ra = Ra)
        for (i in 1:31) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    Sa <- Ia <- Ra<- Ssa <- Isa <- Rsa <- Sc <- Ic <- Rc <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      Sa <- c(Sa, init["Sa"])
      Ia <- c(Ia, init["Ia"])
      Ra <- c(Ra, init["Ra"])
      
      Ssa <- c(Ssa, init["Ssa"])
      Isa <- c(Isa, init["Isa"])
      Rsa <- c(Rsa, init["Rsa"])
      
      Sc <- c(Sc, init["Sc"])
      Ic <- c(Ic, init["Ic"])
      Rc <- c(Rc, init["Rc"])
      
      init <- tmp
    }
    
    #sum population based on column name
    results<-data.frame(time, 
                        Sc,  Ic, Rc, Ssa, Isa, Rsa, Sa, Ia, Ra)%>% 
      dplyr::mutate(N = rowSums(across(-c(time), na.rm=TRUE)))%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
  }

############## MODEL 5 SEIRS - Lumpy skin disease #####
model5=
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 40)
        change <- matrix(0, nrow = 40, ncol = 12)
        
        N <- Sc+Ec+Ic+Rc +Ssa+Esa+Isa+Rsa +Sa+Ea+Ia+Ra
        tau <- 1
        
        #calf
        rate[1] <- mu_b * (Sa+Ea+Ra) 
        change[1, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[2] <- mu_bI * Ia 
        change[2, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[3] <- beta_c * Sc * (Ic+Isa+Ia)/N
        change[3, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[4] <- phi_c * Ec 
        change[4, ] <- c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[5] <-  (1-rho_c) * gamma_c * Ic
        change[5, ] <- c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[6] <-  rho_c * gamma_c * Ic
        change[6, ] <- c(0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[7] <-  omega_c *  Rc
        change[7, ] <- c(1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[8] <-  delta_c * Sc
        change[8, ] <- c(-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)  
        rate[9] <-  delta_c * Ec
        change[9, ] <- c(0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0) 
        rate[10] <-  delta_c * Ic
        change[10, ] <- c(0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0) 
        rate[11] <-  delta_c * Rc
        change[11, ] <- c(0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0) 
        rate[12] <-  mu_c * Sc
        change[12, ] <- c(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)  
        rate[13] <- mu_c * Ec
        change[13, ] <- c(0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[14] <- mu_c * Ic
        change[14, ] <- c(0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[15] <- mu_c * Rc
        change[15, ] <- c(0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[16] <- epsilon * Sc
        change[16, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        
        # saubadult
        rate[17] <- beta_sa * Ssa * (Ic+Isa+Ia)/N
        change[17, ] <- c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0)
        rate[18] <- phi_sa * Esa 
        change[18, ] <- c(0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0)
        rate[19] <-  (1-rho_sa) * gamma_sa * Isa
        change[19, ] <- c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0)
        rate[20] <-  rho_sa * gamma_sa * Isa
        change[20, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0)
        rate[21] <-  omega_sa *  Rsa
        change[21, ] <- c(0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0)
        rate[22] <-  delta_sa * Ssa
        change[22, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0)  
        rate[23] <-  delta_sa * Esa
        change[23, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0) 
        rate[24] <-  delta_sa *  Isa
        change[24, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0) 
        rate[25] <-  delta_sa *  Rsa
        change[25, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1) 
        rate[26] <-  mu_sa * Ssa
        change[26, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0)  
        rate[27] <- mu_sa * Esa
        change[27, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0)
        rate[28] <- mu_sa * Isa
        change[28, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0)
        rate[29] <- mu_sa * Rsa
        change[29, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0)
        rate[30] <- epsilon * Ssa
        change[30, ] <- c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0)
        
        # adult
        rate[31] <- epsilon * Sa
        change[31, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0) 
        rate[32] <- beta_a * Sa * (Ic+Isa+Ia)/N
        change[32, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0)
        rate[33] <- phi_a * Ea 
        change[33, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0)
        rate[34] <-   (1- rho_a) * gamma_a * Ia
        change[34, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1)
        rate[35] <-   rho_a * gamma_a * Ia
        change[35, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0)
        rate[36] <-  omega_a *  Ra
        change[36, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1)
        rate[37] <-  mu_a * Sa
        change[37, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0) 
        rate[38] <- mu_a * Ea
        change[38, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0) 
        rate[39] <- mu_a * Ia
        change[39, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0) 
        rate[40] <- mu_a * Ra
        change[40, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1) 
        
        init <- c(Sc = Sc, Ec = Ec, Ic = Ic, Rc = Rc, 
                  Ssa = Ssa, Esa = Esa, Isa = Isa, Rsa = Rsa,
                  Sa = Sa, Ea = Ea, Ia = Ia, Ra = Ra)
        for (i in 1:40) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    
    Sa <- Ea <- Ia <- Ra <- Ssa <- Esa <- Isa <- Rsa <- Sc <- Ec <- Ic <- Rc <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      Sa <- c(Sa, init["Sa"])
      Ea <- c(Ea, init["Ea"])
      Ia <- c(Ia, init["Ia"])
      Ra <- c(Ra, init["Ra"])
      Ssa <- c(Ssa, init["Ssa"])
      Esa <- c(Esa, init["Esa"])
      Isa <- c(Isa, init["Isa"])
      Rsa <- c(Rsa, init["Rsa"])
      Sc <- c(Sc, init["Sc"])
      Ec <- c(Ec, init["Ec"])
      Ic <- c(Ic, init["Ic"])
      Rc <- c(Rc, init["Rc"])
      
      init <- tmp
    }
    
    #sum population based on column name
    results<-data.frame(time, 
                        Sc, Ec, Ic, Rc ,Ssa, Esa, Isa, Rsa, Sa, Ea, Ia, Ra )%>% 
      dplyr::mutate(N = rowSums(across(-c(time), na.rm=TRUE)))%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
  }

############## MODEL 6 SEIRMS/E - Foot and mouth disease  #####
model6=
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 47)
        change <- matrix(0, nrow = 47, ncol = 14)
        
        N <- Sc+Ec+Ic+Rc +Ssa+Esa+Isa+Rsa +Sa+Ea+Ia+Ra +M +Sm
        tau <- 1
        
        #calf
        rate[1] <- mu_b * (Sa+Ea) 
        change[1, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[2] <- (1-alpha) * mu_bI * Ia 
        change[2, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[3] <- beta_c * Sc * (Ic+Isa+Ia)
        change[3, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[4] <- phi_c * Ec 
        change[4, ] <- c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[5] <-  (1-rho_c) * gamma_c * Ic
        change[5, ] <- c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[6] <-  rho_c * gamma_c * Ic
        change[6, ] <- c(0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[7] <-  omega_c *  Rc
        change[7, ] <- c(1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[8] <-  delta_c * Sc
        change[8, ] <- c(-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)  
        rate[9] <-  delta_c * Ec
        change[9, ] <- c(0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0) 
        rate[10] <-  delta_c * Ic
        change[10, ] <- c(0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0) 
        rate[11] <-  delta_c * Rc
        change[11, ] <- c(0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0) 
        rate[12] <-  mu_c * Sc
        change[12, ] <- c(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)  
        rate[13] <- mu_c * Ec
        change[13, ] <- c(0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[14] <- mu_c * Ic
        change[14, ] <- c(0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[15] <- mu_c * Rc
        change[15, ] <- c(0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[16] <- alpha * mu_bI * Ia
        change[16, ] <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[17] <- mu_b * Ra
        change[17, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
        rate[18] <- omega_m * M
        change[18, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1)
        rate[19] <- mu_c * M
        change[19, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0)
        #Sm go back to Ssa
        rate[20] <- (delta_c * Sc) + (omega_m * M)  
        change[20, ] <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        #Sm go back to Ec
        rate[21] <-  beta_c * Sm * (Ic+Isa+Ia)
        change[21, ] <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[22] <- mu_c * Sm
        change[22, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[23] <- epsilon * Sc  
        change[23, ] <- c(0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0)
        
        # saubadult
        rate[24] <- beta_sa * Ssa * (Ic+Isa+Ia)
        change[24, ] <- c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[25] <- phi_sa * Esa 
        change[25, ] <- c(0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0)
        rate[26] <-  (1-rho_sa) * gamma_sa * Isa
        change[26, ] <- c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0)
        rate[27] <-  rho_sa * gamma_sa * Isa
        change[27, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0)
        rate[28] <-  omega_sa *  Rsa
        change[28, ] <- c(0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0)
        rate[29] <-  delta_sa * Ssa
        change[29, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0)  
        rate[30] <-  delta_sa * Esa
        change[30, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0) 
        rate[31] <-  delta_sa *  Isa
        change[31, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0) 
        rate[32] <-  delta_sa *  Rsa
        change[32, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0) 
        rate[33] <-  mu_sa * Ssa
        change[33, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0)  
        rate[34] <- mu_sa * Esa
        change[34, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[35] <- mu_sa * Isa
        change[35, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0)
        rate[36] <- mu_sa * Rsa
        change[36, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0)
        rate[37] <- epsilon * Ssa
        change[37, ] <- c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
        
        # adult
        rate[38] <- epsilon * Sa
        change[38, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0) 
        rate[39] <- beta_a * Sa * (Ic+Isa+Ia)
        change[39, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0)
        rate[40] <- phi_a * Ea 
        change[40, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0)
        rate[41] <-   (1- rho_a) * gamma_a * Ia
        change[41, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0)
        rate[42] <-   rho_a * gamma_a * Ia
        change[42, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0)
        rate[43] <-  omega_a *  Ra
        change[43, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0)
        rate[44] <-  mu_a * Sa
        change[44, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0) 
        rate[45] <- mu_a * Ea
        change[45, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0) 
        rate[46] <- mu_a * Ia
        change[46, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0) 
        rate[47] <- mu_a * Ra
        change[47, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0) 
        
        init <- c(Sc = Sc, Ec = Ec, Ic = Ic, Rc = Rc, 
                  Ssa = Ssa, Esa = Esa, Isa = Isa, Rsa = Rsa,
                  Sa = Sa, Ea = Ea, Ia = Ia, Ra = Ra, M = M, Sm = Sm)
        for (i in 1:47) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    
    Sa <- Ea <- Ia <- Ra <- Ssa <- Esa <- Isa <- Rsa <- Sc <- Ec <- Ic <- Rc <- M <- Sm <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      Sa <- c(Sa, init["Sa"])
      Ea <- c(Ea, init["Ea"])
      Ia <- c(Ia, init["Ia"])
      Ra <- c(Ra, init["Ra"])
      Ssa <- c(Ssa, init["Ssa"])
      Esa <- c(Esa, init["Esa"])
      Isa <- c(Isa, init["Isa"])
      Rsa <- c(Rsa, init["Rsa"])
      Sc <- c(Sc, init["Sc"])
      Ec <- c(Ec, init["Ec"])
      Ic <- c(Ic, init["Ic"])
      Rc <- c(Rc, init["Rc"])
      M <- c(M, init["M"])
      Sm <- c(Sm, init["Sm"])
      init <- tmp
    }
    #sum population based on column name
    results<-data.frame(time, 
                        Sc, Ec, Ic, Rc, Ssa, Esa, Isa, Rsa, Sa, Ea, Ia, Ra, M, Sm)%>% 
      dplyr::mutate(N = rowSums(across(-c(time), na.rm=TRUE)))%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc,Sm)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
    
  }

############## MODEL 7 SEIRMS/E - Bovine brucellosis #####
model7=
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 47)
        change <- matrix(0, nrow = 47, ncol = 14)
        
        N <- Sc+Ec+Ic+Rc +Ssa+Esa+Isa+Rsa +Sa+Ea+Ia+Ra +M +Sm
        tau <- 1
        
        #calf
        rate[1] <- mu_b * (Sa+Ea) 
        change[1, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[2] <- (1-alpha) * mu_bI * Ia 
        change[2, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[3] <- beta_c * Sc * (Ic+Isa+Ia)/N
        change[3, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[4] <- phi_c * Ec 
        change[4, ] <- c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[5] <-  (1-rho_c) * gamma_c * Ic
        change[5, ] <- c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[6] <-  rho_c * gamma_c * Ic
        change[6, ] <- c(0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[7] <-  omega_c *  Rc
        change[7, ] <- c(1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[8] <-  delta_c * Sc
        change[8, ] <- c(-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)  
        rate[9] <-  delta_c * Ec
        change[9, ] <- c(0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0) 
        rate[10] <-  delta_c * Ic
        change[10, ] <- c(0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0) 
        rate[11] <-  delta_c * Rc
        change[11, ] <- c(0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0) 
        rate[12] <-  mu_c * Sc
        change[12, ] <- c(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)  
        rate[13] <- mu_c * Ec
        change[13, ] <- c(0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[14] <- mu_c * Ic
        change[14, ] <- c(0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[15] <- mu_c * Rc
        change[15, ] <- c(0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[16] <- alpha * mu_bI * Ia
        change[16, ] <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[17] <- mu_b * Ra
        change[17, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
        rate[18] <- omega_m * M
        change[18, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1)
        rate[19] <- mu_c * M
        change[19, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0)
        #Sm go back to Ssa
        rate[20] <- (delta_c * Sc) + (omega_m * M)  
        change[20, ] <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        #Sm go back to Ec
        rate[21] <-  beta_c * Sm * (Ic+Isa+Ia)/N
        change[21, ] <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[22] <- mu_c * Sm
        change[22, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[23] <- epsilon * Sc  
        change[23, ] <- c(0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0)
        
        # saubadult
        rate[24] <- beta_sa * Ssa * (Ic+Isa+Ia)/N
        change[24, ] <- c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[25] <- phi_sa * Esa 
        change[25, ] <- c(0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0)
        rate[26] <-  (1-rho_sa) * gamma_sa * Isa
        change[26, ] <- c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0)
        rate[27] <-  rho_sa * gamma_sa * Isa
        change[27, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0)
        rate[28] <-  omega_sa *  Rsa
        change[28, ] <- c(0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0)
        rate[29] <-  delta_sa * Ssa
        change[29, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0)  
        rate[30] <-  delta_sa * Esa
        change[30, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0) 
        rate[31] <-  delta_sa *  Isa
        change[31, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0) 
        rate[32] <-  delta_sa *  Rsa
        change[32, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0) 
        rate[33] <-  mu_sa * Ssa
        change[33, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0)  
        rate[34] <- mu_sa * Esa
        change[34, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[35] <- mu_sa * Isa
        change[35, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0)
        rate[36] <- mu_sa * Rsa
        change[36, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0)
        rate[37] <- epsilon * Ssa
        change[37, ] <- c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
        
        # adult
        rate[38] <- epsilon * Sa
        change[38, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0) 
        rate[39] <- beta_a * Sa * (Ic+Isa+Ia)/N
        change[39, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0)
        rate[40] <- phi_a * Ea 
        change[40, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0)
        rate[41] <-   (1- rho_a) * gamma_a * Ia
        change[41, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0)
        rate[42] <-   rho_a * gamma_a * Ia
        change[42, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0)
        rate[43] <-  omega_a *  Ra
        change[43, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0)
        rate[44] <-  mu_a * Sa
        change[44, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0) 
        rate[45] <- mu_a * Ea
        change[45, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0) 
        rate[46] <- mu_a * Ia
        change[46, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0) 
        rate[47] <- mu_a * Ra
        change[47, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0) 
        
        init <- c(Sc = Sc, Ec = Ec, Ic = Ic, Rc = Rc, 
                  Ssa = Ssa, Esa = Esa, Isa = Isa, Rsa = Rsa,
                  Sa = Sa, Ea = Ea, Ia = Ia, Ra = Ra, M = M, Sm = Sm)
        for (i in 1:47) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    
    Sa <- Ea <- Ia <- Ra <- Ssa <- Esa <- Isa <- Rsa <- Sc <- Ec <- Ic <- Rc <- M <- Sm <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      Sa <- c(Sa, init["Sa"])
      Ea <- c(Ea, init["Ea"])
      Ia <- c(Ia, init["Ia"])
      Ra <- c(Ra, init["Ra"])
      Ssa <- c(Ssa, init["Ssa"])
      Esa <- c(Esa, init["Esa"])
      Isa <- c(Isa, init["Isa"])
      Rsa <- c(Rsa, init["Rsa"])
      Sc <- c(Sc, init["Sc"])
      Ec <- c(Ec, init["Ec"])
      Ic <- c(Ic, init["Ic"])
      Rc <- c(Rc, init["Rc"])
      M <- c(M, init["M"])
      Sm <- c(Sm, init["Sm"])
      init <- tmp
    }
    
    #sum population based on column name
    results<-data.frame(time, 
                        Sc, Ec, Ic, Rc, Ssa, Esa, Isa, Rsa, Sa, Ea, Ia, Ra,  M, Sm )%>% 
      dplyr::mutate(N = rowSums(across(-c(time), na.rm=TRUE)))%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc,Sm)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
    
  }


############# SET up population parameters before running #####
#gaur population 
N = 300 

#estimate the age structure proportion
#calf:subadult:adult ratio
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

end.time = 100 * 365 #duration of prediction
n_rep = 100

############## SINGLE RUNS with plotting output - ALL MODELS  ###############################
# > MODEL 1 - NO infection #####
initials_m1 <- c(c = c, sa = sa, a = a )
parameters_m1 <- c(mu_b = 0.34/365, 
                mu_c = 0.27/365, 
                mu_sa = 0.15/365,
                mu_a = 0.165/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials_m1), 
                tau = 1)

res_model1 <- model1(pars = parameters_m1, init = initials_m1,
                end.time = end.time)

PlotMods(res_model1)

# > MODEL 2 SI - Anthrax #####
initials_m2 <- c(Sc = c, Ic = 0, Ssa = sa, Isa = 0, Sa = (a-1), Ia = 1 )
parameters_m2 <- c(
  beta_c = 0.0001, beta_sa = 0.0001, beta_a = 0.0001,
  gamma_c = 1, gamma_sa = 1, gamma_a = 1, #assume that animal die in 1 day
  rho_c = 1, rho_sa = 1,  rho_a = 1,
  epsilon = 2e-5,
  mu_b = 0.34/365, 
  mu_c = 0.27/365, mu_sa = 0.15/365, mu_a = 0.165/365,
  delta_c = 1/365, delta_sa = 1/(3*365),
  N = sum(initials_m2),
  tau=1)

res_model2<-model2(pars = parameters_m2, init = initials_m2,
       end.time = end.time)
PlotMods(res_model2)

# >  MODEL 3 SEI - Bovine tuberculosis #####
initials_m3 <- c(Sc = c, Ec = 0, Ic = 0, Ssa = sa, Esa = 0, Isa = 0, Sa = (a-1), Ea = 0, Ia = 1)
parameters_m3 <- c( 
  beta_c = 0.043/30, beta_sa = 0.043/30, beta_a = 0.043/30,
  phi_c = 0.21/30, phi_sa = 0.21/30, phi_a = 0.21/30,
  gamma_c = 0, gamma_sa = 0, gamma_a = 0,
  rho_c = 0, rho_sa = 0, rho_a = 0.1, 
  epsilon = 2e-5,
  mu_b = 0.34/365,  mu_bI = (0.34/365)*(1-0.27), #Ia birth rate reduce by = 27%   
  mu_c = 0.27/365,  mu_sa = 0.15/365, mu_a = 0.165/365,
  delta_c = 1/365, delta_sa = 1/(3*365),
  N = sum(initials_m3), tau=1)

res_model3<-model3(pars = parameters_m3, init = initials_m3,
                   end.time = end.time)
PlotMods(res_model3)

# >  MODEL 4 SIRS - Hemorrhagic septicemia #####
initials_m4 <- c(Sc = c,  Ic = 0, Rc = 0, Ssa = sa, Isa = 0, Rsa = 0, Sa = (a-1), Ia = 1, Ra = 0)
parameters_m4 <- c( 
  beta_c = 0.33/365,beta_sa = 0.33/365,beta_a = 0.33/365,
  gamma_c  =1/3, gamma_sa =1/3, gamma_a =1/3,
  rho_c = 0.9,rho_sa = 0.43,rho_a = 0.43,
  omega_c = 1/180,omega_sa = 1/180, omega_a = 1/180,
  epsilon = 2e-5,
  mu_b = 0.34/365, 
  mu_c = 0.27/365, 
  mu_sa = 0.15/365,
  mu_a = 0.165/365,
  delta_c = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1
)

res_model4<-model4(pars = parameters_m4, init = initials_m4,
                   end.time = end.time)
PlotMods(res_model4)

# >  MODEL 5 SEIRS - Lumpy skin disease #####
# FIX PARAMATERS before running#
initials_m5 <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = (a-1), Ea = 0, Ia = 1, Ra = 0)
parameters_m5 <- c( 
  beta_c = 0.038, beta_sa = 0.038, beta_a = 0.038, 
  phi_c = 1/7,phi_sa = 1/7,phi_a = 1/7,
  gamma_c = 1/35, gamma_sa = 1/35,gamma_a = 1/35,
  rho_c = 0.05,rho_sa = 0.03, rho_a = 0.01, 
  omega_c = 1/180, omega_sa =  1/180, omega_a = 1/180,
  epsilon = 2e-5,
  mu_b = 0.34/365,  mu_bI = (0.34/365)*(0.9), #Ia birth rate reduce by = 27%   
  mu_c = 0.27/365,  mu_sa = 0.15/365, mu_a = 0.165/365,
  delta_c = 1/365, delta_sa = 1/(3*365),
  N = sum(initials_m5),
  tau=1)

res_model5<-model5(pars = parameters_m5, init = initials_m5,
                   end.time = end.time)
PlotMods(res_model5)

# > MODEL 6 SEIRMS/E - Foot and mouth disease  #####
initials_m6 <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = a, Ea = 0, Ia = 1, Ra = 0)
parameters_m6 <- c( 
  beta_c = 0.52/365, beta_sa = 0.52/365, beta_a = 0.52/365,
  phi_c = 1/8, phi_sa = 1/6, phi_a = 1/6,
  gamma_c = 1/5, gamma_sa = 1/5, gamma_a = 1/5,
  rho_c = 0.1, rho_sa = 0.05, rho_a = 0.03, 
  alpha = 0.5,
  omega_c = (1/120), omega_sa =  (1/120), omega_a = (1/565), omega_m = (1/144),
  epsilon = 2e-5,
  mu_b = 0.34/365, mu_bI = (0.34/365)*(0.9), #Ia birth rate reduce by = 10%  (assume)
  mu_c = 0.27/365, mu_sa = 0.15/365, mu_a = 0.165/365,
  delta_c = 1/365, delta_sa = 1/(3*365),
  N = sum(initials_m6),
  tau=1)

res_model6<-model6(pars = parameters_m6, init = initials_m6,
                   end.time = end.time)
PlotMods(res_model6)

# > MODEL 7 SEIRMS/E - Brucellosis  #####
initials_m7 <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, 
                 Sa = (a-1), Ea = 0, Ia = 1, Ra = 0)
parameters_m7 <- c( 
  beta_c = 2/365, beta_sa = 2/365, beta_a = 2/365,
  phi_c = 1/14, phi_sa = 1/14, phi_a = 1/14,
  gamma_c = 1/(2*365), gamma_sa = 1/(2*365), gamma_a = 1/(2*365),
  rho_c = 0.1, rho_sa = 0.05, rho_a = 0.03, 
  alpha = 0.9,   
  omega_c = 1/180, omega_sa =  1/180, omega_a = 1/180, omega_m = 1/180,
  epsilon = 2e-5,
  mu_b = 0.34/365, 
  mu_bI = (0.34/365)*(0.5), #Ia birth rate reduce by = 10%  (assume)
  mu_c = 0.27/365, 
  mu_sa = 0.15/365,
  mu_a = 0.165/365,
  delta_c = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m7),
  tau=1
)
res_model7<-model7(pars = parameters_m7, init = initials_m7,
                   end.time = end.time)
PlotMods(res_model7)

## MULTIPLE RUNS - ALL MODELS ################################
#### set end.time and replication times 
end.time <- 100 * 365
n_rep <- 100

# use the parameters same as a single run

# > MX RUNS MODEL 1 - NO infection #####
sim_rep_m1<-replicate(n_rep,(model1(pars = parameters_m1, init = initials_m1,
                                    end.time = end.time)))

# > MX RUNS MODEL 2 SI - Anthrax #####
sim_rep_m2<-replicate(n_rep,(model2(pars = parameters_m2, init = initials_m2,
                                    end.time = end.time)))

# > MX RUNS MODEL 3 SEI - Bovine tuberculosis #####
sim_rep_m3<-replicate(n_rep,(model3(pars = parameters_m3, init = initials_m3,
                                    end.time = end.time)))

# > MX RUNS MODEL 4 SIRS - Hemorrhagic septicemia #####
sim_rep_m4<-replicate(n_rep,(model4(pars = parameters_m4, init = initials_m4,
                                    end.time = end.time)))

# > MX RUNS MODEL 5 SEIRS - Lumpy skin disease #####
sim_rep_m5<-replicate(n_rep,(model5(pars = parameters_m5, init = initials_m5,
                                    end.time = end.time)))

# > MX RUNS MODEL 6 SEIRMS/E - Foot and mouth disease #####
sim_rep_m6<-replicate(n_rep,(model6(pars = parameters_m6, init = initials_m6,
                                    end.time = end.time)))

# > MX RUNS MODEL 7 SEIRMS/E - Brucellosis #####
sim_rep_m7<-replicate(n_rep,(model7(pars = parameters_m7, init = initials_m7,
                                    end.time = end.time)))

# FUNCTION for plotting model's simulations #####
single_pop_sim_prep <- function(x, n_rep, end.time, melt){ # x = simulation of model, e.g. sim_run_m1
  df<-list()
  #loop for storing new df
  for (i in 1:n_rep){
    run <- paste("run", seq(n_rep), sep="")
    names(df)[i]
    df[[i]]<- x[,i]$results[,-c(1)]
    df[[i]]$time_d <- seq(from=1,to=end.time+1,by=1)
    
  }
  
  df<-map2(df,run, ~cbind(.x, run = .y))   # adding n_rep to the column
  df2<- data.table::rbindlist(df)          # binding row
  
  if (melt == T) {  #option for melting the data in case we need...
    
    df3 <- gather(df2, key = class, value = value, -c(time,run))
    
    return(df3 = data.frame(df3))
  } 
  
  else  {
    return( df2 = data.frame(df2)) 
  }
  
}

#creating the model simulation list 
sim_rep_m<-list(sim_rep_m1,
                sim_rep_m2,
                sim_rep_m3,
                sim_rep_m4,
                sim_rep_m5,
                sim_rep_m6,
                sim_rep_m7)
#creating name list for each model (in an order)
# from sim_rep_m1 to sim_rep_m7
nam<-c('pop_dynamic',
       'Anthrax',
       'bTB',
       'HS',
       'LSD',
       'FMD',
       'Brucellosis')
nam<-c('LSD')
m<-list()

#group and calculate total population change (%) loop--------
for (i in 1:length(sim_rep_m)) {
  m[[i]]<- single_pop_sim_prep(x = sim_rep_m[[i]], n_rep=n_rep, end.time= end.time, melt = F)
  m[[i]]<- m[[i]]%>%
    group_by(run) %>%
    mutate(Ndiff = ((last(N)-first(N))/first(N))*100)%>% #calculate change in the total population at year100, and year0
    mutate(time_y = time_d/365) %>% #convert day to year for plotting
    as.data.frame()
  
  m[[i]]$model <- paste0(nam[[i]])
  
}

#loop for saving the data.frame as .rds
for (i in 1:length(m)) {
  saveRDS(m[[i]], file = paste0("df_",nam[[i]],".rds")) 
}

## PLOT MODEL OUTPUTS: population line graphs #############################


## SET UP FUNCTIONS FOR MX METRICS ################################

############## FUNCTION  1 COUNT EXTINCTIONS ########
# minimum time that I == 0 
my_min_ext_I<-function(x=get_time$results, y=get_time$results$I,...){
  res_no<-vector()
  if(!is.na(min(subset(x,y==0)$time))){ #if there is no missing value
    res_no = min(subset(x,y==0)$time)   # res_no will collect the minimum time (e.g. unit = day) that I go to 0
  } else{
    res_no = end.time                   # else, res_no = the end.time
  }
  res_no
}




# minimum time that N == 0 
my_min_ext_N<-function(x=get_time$results, y=get_time$results$N,...){
  res_no<-vector()
  if(!is.na(min(subset(x,y==0)))){
    res_no = min(subset(x,y==0)$time)
  } else{
    res_no = end.time
  }
  res_no
}

############## FUNCTION  2 COUNT EXTINCTIONS  #########
#count I == 0 as 1, count I >=1 as 0 at time [i]

#1 = no infectious animal in the herd (count the event of disease extinction in the herd)
#0 = more than one infectious animal in the herd

my_imp_ext_I<-function(x=time,y=I,...){
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

#count N == 0 as 1, count I >=1 as 0 at time [i]
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
# at time (i-1), and (i) I or N  == 0 , count == 1 
# 1 = disease/population extinct from the population at time i
# 0 = no disease/population extinction 
# same as function #2 but the result is  == length(res_no[!is.na(res_no)])

# disease persistence
my_imp_ext_na_I<-function(x=time,y=I,...){
  
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

# population persistence 
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

## SET UP FUNCTIONS FOR MODEL RUNS LOOP THROUGH PARS ################################
#### set initial values - change for populations

############## FUNCTION  6  MODEL PAR LOOPS MIN TIME & COUNT I EXTINCTIONS #######

model2_anth_extinct <-function() {
  for (k in 1:length(beta)){
  for (i in 1:length(rho)) {
    parameters <- c(beta_a = beta[k], beta_sa = beta[k],  beta_c = beta[k],
                    gamma_c = (1/(1/24))/365,gamma_sa = (1/(1/24))/365, gamma_a = (1/(1/24))/365,
                    rho_c= rho[i], rho_sa= rho[i],rho_a= rho[i],
                    mu_b = 0.34/365,  mu_c = 0.27/365, mu_sa = 0.15/365, mu_a = 0.165/365,
                    delta_c = 1/365,delta_sa = 1/(3*365), 
                    epsilon = 2e-5,
                    N = sum(initials), 
                    tau=1)
    
    initials <- c(Sc = c, Ic = 0, Ssa = sa, Isa = 0, Sa = (a-1), Ia = 1 )
    
    get_time <- model2(pars = parameters, init = initials,
                                   end.time = end.time)
    
    res_min_N[k,i]<- my_min_ext_N(x=get_time$results,y=get_time$results$N)
    res_num_ext_N[k,i]<- my_imp_ext_N(x=get_time$time,y=get_time$results$N)
    res_time_inf_N[k,i]<- my_imp_ext_na_N(x=get_time$time,y=get_time$results$N)
    
    res_min_I[k,i]<- my_min_ext_I(x=get_time$results,y=get_time$results$I)
    res_num_ext_I[k,i]<- my_imp_ext_I(x=get_time$time,y=get_time$results$I)
    res_time_inf_I[k,i]<- my_imp_ext_na_I(x=get_time$time,y=get_time$results$I)
  }}
  d <- list(res_min_N,res_num_ext_N,res_time_inf_N,res_min_I,res_num_ext_I,res_time_inf_I)
  d
}

model3_tb_extinct<-function() {
  for (k in 1:length(beta)){
    for (i in 1:length(rho)) {
     
       parameters <- c( beta_a = beta[k], beta_sa = beta[k], beta_c = beta[k],
                        phi_c = 0.21/30, phi_sa = 0.21/30, phi_a = 0.21/30,
                        gamma_c = 0, gamma_sa = 0, gamma_a = 0,
                        rho_c= rho[i], rho_sa= rho[i],rho_a= rho[i],
                        mu_b = 0.34/365,  mu_bI = (0.34/365)*(1-0.27), #Ia birth rate reduce by = 27%   
                        mu_c = 0.27/365,  mu_sa = 0.15/365, mu_a = 0.165/365,
                        delta_c = 1/365, delta_sa = 1/(3*365),
                        epsilon = 2e-5,
                        N = sum(initials), tau=1)
      
      initials <- c(Sc = c, Ec = 0, Ic = 0, Ssa = sa, Esa = 0, Isa = 0, Sa = a, Ea = 0, Ia = 1)
      get_time <- model3(pars = parameters, init = initials,
                         end.time = end.time)
      
      res_min_N[k,i]<- my_min_ext_N(x=get_time$results,y=get_time$results$N)
      res_num_ext_N[k,i]<- my_imp_ext_N(x=get_time$time,y=get_time$results$N)
      res_time_inf_N[k,i]<- my_imp_ext_na_N(x=get_time$time,y=get_time$results$N)
      
      res_min_I[k,i]<- my_min_ext_I(x=get_time$results,y=get_time$results$I)
      res_num_ext_I[k,i]<- my_imp_ext_I(x=get_time$time,y=get_time$results$I)
      res_time_inf_I[k,i]<- my_imp_ext_na_I(x=get_time$time,y=get_time$results$I)
    }}
  d <- list(res_min_N,res_num_ext_N,res_time_inf_N,res_min_I,res_num_ext_I,res_time_inf_I)
  d
}

model4_hs_extinct<-function() {
  for (k in 1:length(beta)){
    for (i in 1:length(rho)) {
      
      parameters <- c( beta_a = beta[k], beta_sa = beta[k], beta_c = beta[k],
                       gamma_c  =1/3/365, gamma_sa =1/3/365, gamma_a =1/3/365,
                       rho_c= rho[i], rho_sa= rho[i], rho_a= rho[i],
                       omega_c = 1/180/365, omega_sa = 1/180/365, omega_a = 1/180/365,
                       epsilon = 2e-5,
                       mu_b = 0.34/365, mu_c = 0.27/365, mu_sa = 0.15/365, mu_a = 0.165/365,
                       delta_c = 1/365, delta_sa = 1/(3*365),
                       N = sum(initials),
                       tau=1)
      
      initials <- c(Sc = c,  Ic = 0, Rc = 0, Ssa = sa, Isa = 0, Rsa = 0, Sa = a, Ia = 1, Ra = 0)
      
      get_time <- model4_ep(pars = parameters, init = initials,
                            end.time = end.time)
      res_min_N[k,i]<- my_min_ext_N(x=get_time$results,y=get_time$results$N)
      res_num_ext_N[k,i]<- my_imp_ext_N(x=get_time$time,y=get_time$results$N)
      res_time_inf_N[k,i]<- my_imp_ext_na_N(x=get_time$time,y=get_time$results$N)
      
      res_min_I[k,i]<- my_min_ext_I(x=get_time$results,y=get_time$results$I)
      res_num_ext_I[k,i]<- my_imp_ext_I(x=get_time$time,y=get_time$results$I)
      res_time_inf_I[k,i]<- my_imp_ext_na_I(x=get_time$time,y=get_time$results$I)
    }}
  d <- list(res_min_N,res_num_ext_N,res_time_inf_N,res_min_I,res_num_ext_I,res_time_inf_I)
  d
}

model5_lsd_extinct<-function() {
  for (k in 1:length(beta)){
    for (i in 1:length(rho)) {
      
      parameters <- c( beta_a = beta[k], beta_sa = beta[k], beta_c = beta[k],
                       phi_c = 1/7/365, phi_sa = 1/7/365,phi_a = 1/7/365,
                       gamma_c = 0, gamma_sa = 0,gamma_a = 0,
                       rho_c= rho[i], rho_sa= rho[i], rho_a= rho[i],
                       omega_c = (1/365)/365, omega_sa =  (1/365)/365, omega_a = (1/365)/365,
                       mu_b = 0.34/365,  mu_bI = (0.34/365)*(1-0.27), 
                       mu_c = 0.27/365,  mu_sa = 0.15/365, mu_a = 0.165/365,
                       delta_c = 1/365, delta_sa = 1/(3*365),
                       epsilon = 2e-5,
                       N = sum(initials),
                       tau=1)
      
      c(Sc = c, Ec = 0, Ic = 0, Rc = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = a, Ea = 0, Ia = 1, Ra = 0)
      
      get_time <- model5(pars = parameters, init = initials,
                            end.time = end.time)
      res_min_N[k,i]<- my_min_ext_N(x=get_time$results,y=get_time$results$N)
      res_num_ext_N[k,i]<- my_imp_ext_N(x=get_time$time,y=get_time$results$N)
      res_time_inf_N[k,i]<- my_imp_ext_na_N(x=get_time$time,y=get_time$results$N)
      
      res_min_I[k,i]<- my_min_ext_I(x=get_time$results,y=get_time$results$I)
      res_num_ext_I[k,i]<- my_imp_ext_I(x=get_time$time,y=get_time$results$I)
      res_time_inf_I[k,i]<- my_imp_ext_na_I(x=get_time$time,y=get_time$results$I)
    }}
  d <- list(res_min_N,res_num_ext_N,res_time_inf_N,res_min_I,res_num_ext_I,res_time_inf_I)
  d
}

model6_fmd_extinct<-function() {
  for (k in 1:length(beta)){
    for (i in 1:length(rho)) {
      
      parameters <- c( beta_a = beta[k], beta_sa = beta[k], beta_c = beta[k],
                       phi_c = 1/8/365, phi_sa = 1/6/365, phi_a = 1/6/365,
                       gamma_c = 1/5/365, gamma_sa = 1/5/365, gamma_a = 1/5/365,
                       rho_c= rho[i], rho_sa= rho[i], rho_a= rho[i],
                       omega_c = (1/120)/365, omega_sa =  (1/120)/365, omega_a = (1/565)/365, omega_m = (1/144)/365,
                       epsilon = 2e-5,
                       mu_b = 0.34/365, mu_bI = (0.34/365)*(0.9), #Ia birth rate reduce by = 10%  (assume)
                       mu_c = 0.27/365, mu_sa = 0.15/365, mu_a = 0.165/365,
                       delta_c = 1/365, delta_sa = 1/(3*365),
                       N = sum(initials),
                       tau=1)
      
      initials <- c c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = a, Ea = 0, Ia = 1, Ra = 0)
      
      get_time <- model6(pars = parameters, init = initials,
                            end.time = end.time)
      
      res_min_N[k,i]<- my_min_ext_N(x=get_time$results,y=get_time$results$N)
      res_num_ext_N[k,i]<- my_imp_ext_N(x=get_time$time,y=get_time$results$N)
      res_time_inf_N[k,i]<- my_imp_ext_na_N(x=get_time$time,y=get_time$results$N)
      
      res_min_I[k,i]<- my_min_ext_I(x=get_time$results,y=get_time$results$I)
      res_num_ext_I[k,i]<- my_imp_ext_I(x=get_time$time,y=get_time$results$I)
      res_time_inf_I[k,i]<- my_imp_ext_na_I(x=get_time$time,y=get_time$results$I)
    }}
  d <- list(res_min_N,res_num_ext_N,res_time_inf_N,res_min_I,res_num_ext_I,res_time_inf_I)
  d
}

model7_bru_extinct<-function() {
  for (k in 1:length(beta)){
    for (i in 1:length(rho)) {
      
      parameters <- c(beta_c = 2/365, beta_sa = 2/365,beta_a = 2/365,
                      phi_c = 1/14/365, phi_sa = 1/14/365, phi_a = 1/14/365,
                      gamma_c = 1/2*365, gamma_sa = 1/2*365, gamma_a = 1/2*365,
                      rho_c= rho[i], rho_sa= rho[i], rho_a= rho[i],
                      alpha = 0.9,
                      omega_c = (1/180)/365, omega_sa =  (1/180)/365, omega_a = (1/180)/365, omega_m = (1/180)/365,
                      epsilon = 2e-5,
                      mu_b = 0.34/365, mu_bI = (0.34/365)*(0.5), #Ia birth rate reduce by = 50%  (assume)
                      mu_c = 0.27/365, mu_sa = 0.15/365, mu_a = 0.165/365,
                      delta_c = 1/365, delta_sa = 1/(3*365),
                      N = sum(initials),
                      tau=1)
      
      initials <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = a, Ea = 0, Ia = 1, Ra = 0)
      get_time <- model6(pars = parameters, init = initials,
                         end.time = end.time)
      
      res_min_N[k,i]<- my_min_ext_N(x=get_time$results,y=get_time$results$N)
      res_num_ext_N[k,i]<- my_imp_ext_N(x=get_time$time,y=get_time$results$N)
      res_time_inf_N[k,i]<- my_imp_ext_na_N(x=get_time$time,y=get_time$results$N)
      
      res_min_I[k,i]<- my_min_ext_I(x=get_time$results,y=get_time$results$I)
      res_num_ext_I[k,i]<- my_imp_ext_I(x=get_time$time,y=get_time$results$I)
      res_time_inf_I[k,i]<- my_imp_ext_na_I(x=get_time$time,y=get_time$results$I)
    }}
  
  d <- list(res_min_N,res_num_ext_N,res_time_inf_N,res_min_I,res_num_ext_I,res_time_inf_I)
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
single_pop_sim_prep_I <- function(x, n_rep, end.time){ # x = simulation of model, e.g. sim_run_m1
  mat = matrix(NA, nrow=n_rep, ncol = end.time+1)
  for (i in 1:n_rep){
    mat[i,]<-x[,i]$results$I
  }
  colnames(mat) = paste("time", seq(from=1,to=end.time+1,by=1), sep="")
  rownames(mat) = paste("run", seq(n_rep), sep="")
  dat = as.data.frame(mat)
  dat$run = rownames(dat)
  mdat = melt(dat, id.vars="run")
  mdat$time = as.numeric(gsub("time", "", mdat$variable))
  mdat
}

single_pop_sim_prep_N <- function(x, n_rep, end.time){ # x = simulation of model, e.g. sim_run_m1
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

#Anthrax
beta <- seq(from = 5e-5,to = 0.01, length.out=10)
rho <- seq(from = 0.1, to = 1, length.out=10)

#TB
beta <- seq(from = 0.0003,to = 0.0276, length.out=10)
rho <- seq(from = 0, to = 0.11, length.out=10)

d <- as.data.frame(matrix(NA, length(beta),length(rho)))
res_min<-array(unlist(d), dim=c(length(beta), length(rho)))
res_num_ext<-array(unlist(d), dim=c(length(beta), length(rho)))
res_time_inf<-array(unlist(d), dim=c(length(beta), length(rho)))

n_rep = 100
end.time = 100 * 365

############## RUN Big MODEL #####
big_model2<-replicate(n_rep,my_fun_model2_ant_ep())

############## RUN MODEL 2 #####

big_run_model2<-replicate(n_rep,my_fun_model2())

############## RUN MODEL 3 ##############

big_run_model3<-replicate(n_rep,my_fun_model3())

############## RUN MODEL 4 ###########

big_run_model4<-replicate(n_rep,my_fun_model4())

############## RUN MODEL 5 ###########

big_run_model5<-replicate(n_rep,my_fun_model5())


## MODEL OUTPUTS PREPARATION #######################

############## MODELS 1-6 OUTPUT PREPARATION #############################
mod_res<-list(big_model2_ant_no_ep_N,
              big_model2_ant_ep)
              

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
  assign(paste0("df_min_", i), my_plot_min(x=plot_res_min[[i]],par1 = beta, par2 = rho, par1_n = 'beta', par2_n = 'rho'))
}

for (i in 1:length(plot_res_num)){
  assign(paste0("df_num_", i), my_plot_num(x=plot_res_num[[i]],par1 = beta,par2 = rho, par1_n = 'beta', par2_n = 'rho'))
}

for (i in 1:length(plot_res_time)){
  assign(paste0("df_time_", i), my_plot_time(x=plot_res_time[[i]],par1 = beta,par2 = rho, par1_n = 'beta', par2_n = 'rho'))
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
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)
    }
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
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)
    }
  else
  {
    pname <- paste0("time-",k)
    p<-(ggplot(subdata, aes(x = beta, y = rho, fill = value))+
          geom_tile()+
          scale_fill_gradientn(colors=colorRampPalette(c("brown","red","orange","seagreen","royalblue","whitesmoke"))(500),
                               name="Average\ntime", limits = c(0,1)) +
          labs(x = expression(beta),y=expression(rho)) +
          theme_bw())
    ggsave(paste0(pname,".pdf"),p, width = 4.5, height = 3.9)
    }
}

## SINGLE POPULATION PRINT PREPS #######

res_mx<-list(sim_rep_m2,
              sim_rep_m3,
              sim_rep_m4,
              sim_rep_m5,
              sim_rep_m6,
              sim_rep_m7)

res_mx<-list(sim_rep_m2,
             sim_rep_m3)
#change single_pop_sim_prep from I to N

for (i in 1:length(res_mx)){
  assign(paste0("df_mx_I", i), single_pop_sim_prep(x=res_mx[[i]],n_rep = n_rep, end.time = end.time))
}

for (i in 1:length(res_mx)){
  assign(paste0("df_mx_I", i), single_pop_sim_prep(x=res_mx[[i]],n_rep = n_rep, end.time = end.time))
}
df_mx_I_m3<-single_pop_sim_prep_I(x=sim_rep_m3,n_rep = n_rep, end.time = end.time)
df_mx_N_m3<-single_pop_sim_prep_N(x=sim_rep_m3,n_rep = n_rep, end.time = end.time)

res_mx_p<-rbind( df_mx1,
                df_mx2,
                df_mx3,
                df_mx4,
                df_mx5,
                df_mx6)

res_mx_p<-rbind(df_mx1,
                df_mx2)

res_mx_p$model<-c(rep('1',dim(df_mx1)[1]),
                  rep('2',dim(df_mx2)[1]),
                  rep('3',dim(df_mx3)[1]),
                  rep('4',dim(df_mx4)[1]),
                  rep('5',dim(df_mx5)[1]),
                  rep('6',dim(df_mx6)[1]))

res_mx_p$model<-c(rep('1',dim(df_mx1)[1]),
                  rep('2',dim(df_mx2)[1]))
## SINGLE POPULATION PLOTS #######

for (i in unique(res_mx_p$model)){
  subdata <- subset(res_mx_p, model == 1)
  pdf(paste("plot_ts_tb_I.pdf", sep = ""), width = 4, height = 3)
  print(ggplot(df_mx_I_m3, aes(x=time, y=value, group=run)) +
          theme_bw() +
          theme(panel.grid=element_blank()) +
          geom_line(size=0.2, alpha=0.15)+
          ylab('I (Numbers)') + xlab('time')+
          ggtitle("bTB infectious of gaur population, 100 runs")+
          theme(plot.title = element_text(size=8))+
          
          stat_summary(aes(group = 1), fun=mean, geom="line", colour="black",size = 1.1))
  dev.off()
}

#for (i in unique(res_mx_p$model)){
#  subdata <- subset(res_mx_p, model == 2)
  pdf(paste("plot_ts_btb_N.pdf", sep = ""), width = 4, height = 3)
  print(ggplot(df_mx_N_m3, aes(x=time, y=value, group=run)) +
          theme_bw() +
          theme(panel.grid=element_blank()) +
          geom_line(size=0.2, alpha=0.15)+
          ylab('N (Numbers)') + xlab('time')+
          ggtitle("gaur total population with bTB infections, 100 runs")+
          theme(plot.title = element_text(size=8))+
          
          stat_summary(aes(group = 1), fun=mean, geom="line", colour="black",size = 1.1))
  dev.off()
#}


############### PLOT ALL TS SINGLE RUNS ##############

res_p_ts<-list(sim_rep_m2,
               sim_rep_m3,
               sim_rep_m4,
               sim_rep_m5,
               sim_rep_m6,
               sim_rep_m7)
res_p_ts<-list(sim_rep_m2,
                 sim_rep_m3)
# Make plots. Single runs
for (i in 1:length(res_p_ts)) {
  pdf(paste("plotts_all", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot() + 
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$N), color = "blue", size =1.2) +
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$S), color = "black", size =1.2) +
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$E), color = "yellow", size =1.2) +
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$I), color = "red", size =1.2) +
          #geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$R), color = "seagreen", size =1.2, linetype = "dotted") +
          xlab('Time') +
          ylab('Infection State Numbers'))
  dev.off()
}

#pdf(paste("plotts_all_anth.pdf", sep = ""), width = 4, height = 3)
png("gaur_btb_100y_all.png",width = 25, height = 15, units = 'cm', res = 600)
print(ggplot() +
        geom_line(data = res_p_ts[[2]][[4]], aes(x = res_p_ts[[2]][[4]]$time, y = res_p_ts[[2]][[4]]$S), color = "seagreen4", size =1.2) +
        geom_line(data = res_p_ts[[2]][[4]], aes(x = res_p_ts[[2]][[4]]$time, y = res_p_ts[[2]][[4]]$E), color = "darkorange2", size =1.2) +
        geom_line(data = res_p_ts[[2]][[4]], aes(x = res_p_ts[[2]][[4]]$time, y = res_p_ts[[2]][[4]]$I), color = "firebrick", size =1.2) +
        #geom_line(data = res_p_ts[[1]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$R), color = "seagreen", size =1.2, linetype = "dotted") +
        geom_line(data = res_p_ts[[2]][[4]], aes(x = res_p_ts[[2]][[4]]$time, y = res_p_ts[[2]][[4]]$N), color = "#153030", size =1.2) +
        xlab('Time') +
        ylab('Infection State Numbers'))+
        ggtitle("Gaur population with bTB, 100 years")+
  scale_color_manual( name = "class",
                      #labels = c("I"),
                      labels = c('S','E','I', 'total'),
                      values = c('S'='seagreen4',
                                  'E'='darkorange2',
                                  'I'='firebrick',
                                   #"R"='dodgerlblue3',
                                  "total"='#153030'))+ #blackgreen
                       
        theme_bw() +
        theme( plot.title = element_text(size = 18),
         axis.title.x = element_text(size = 15),
         axis.title.y = element_text(size = 15),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=13))
          guides(color = guide_legend(override.aes = list(alpha = 1,size=0.5)))

dev.off()

############### PLOT ALL TS SINGLE RUNS I ONLY ##############

# Make plots. Single runs
for (i in 1:length(res_p_ts)) {
  pdf(paste("plotts_i_I", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot() + 
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$I), color = "red", size =1.2) +
          xlab('Time') +
          ylab('Infection Numbers'))
  dev.off()
}

# Make plots. Single runs - SCALED
for (i in 1:length(res_p_ts)) {
  pdf(paste("plotts_i_scale_I", i, ".pdf", sep = ""), width = 4, height = 3)
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

df.agg <- aggregate(time ~ run + value + model, sim_rep_m3, min)

df.ag<-(df.agg[df.agg$value==0,c('model','time')])

neworder <- c("1","2","3","4","5","6")
library(plyr)  ## or dplyr (transform -> mutate)
df.ag <- arrange(transform(df.ag,
                           model=factor(model,levels=neworder)),model)
labs <- c('1' = "SI",
          '2' = "SEI",
          '3' = "SIRS+",
          '4' = "SEIRS",
          '5' = "SEIRM-FMD",
          '6' = "SEIRM-Bru")

p<-ggplot(df.ag, aes(x=time))+
  geom_histogram(color="black", fill="grey")+
  facet_wrap(model~., ncol = 6, labeller = labeller(model = labs))+
  scale_x_continuous(breaks = c(0, 800, 1600), labels = c("0", "800", "1600"))
pdf("extinctions.pdf", width = 8, height = 3)
p
dev.off()

p_yr<-ggplot(df.ag, aes(x=time))+
  geom_histogram(color="black", fill="grey")+
  facet_wrap(model~., ncol = 6,labeller = labeller(model = labs))+
  scale_x_continuous(limits = c(0,1000),breaks = c(0, 400, 800), labels = c("0", "400", "800")) +
  ylim(0,60)
pdf("extinctions_1yr.pdf", width = 5, height = 3)
p_yr
dev.off()
