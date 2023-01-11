############## NOTES ##################
# Infectious disease modelling code for wild Bovidae population
# - 18 Nov 2022-

# # Reference R Code: https://github.com/dtsh2/ebola_model; 
# # Reference Article:https://doi.org/10.1098/rsif.2021.0638  (Transmission models indicate Ebola virus persistence in non-human primate populations is unlikely)

# Working on:  R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
# RStudio 2022.07.1+554 "Spotted Wakerobin" Release (7872775ebddc40635780ca1ed238934c3345c5de, 2022-07-22) for macOS
# Mozilla/5.0 (Macintosh; Intel Mac OS X 13_0_1)

# All models have 3 age classes (c = calf, a = Subadult, a = adult)
# unit == day (per day)
# This code set up as:
# 1) Set up the infectious disease model's function (7 models, details below)
# 2) Set up parameters
# 3) Single run for all the models
# 4) Multiple run for all the models 
# 5) Plotting single and multiple runs

rm(list=ls())
set.seed(111)
# set library path
# path <- '...'

#install.packages(c('EpiDynamics','dplyr','tidyverse','reshape2',
#                   'stringr','hrbrthemes','viridis','ggstatsplot'), 
#                 lib = path)

############## LOAD PACKAGES ##########
library(EpiDynamics)
library(dplyr)   
library(tidyverse)
library(reshape2) 
library(stringr)
library(ggplot2); theme_set(theme_bw())
library(hrbrthemes)
library(viridis)
library(ggstatsplot)


# Frequency dependent transmission

############### 1) Set up the infectious disease model's function ############### 

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
    results<-data.frame(time, a,  sa,  c)%>% 
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
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(N = rowSums(across(-c(time,S,I), na.rm=TRUE)))
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
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(N = rowSums(across(c(S,E,I), na.rm=TRUE))) 
    
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
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))%>%
      dplyr::mutate(N = rowSums(across(c (S,I,R), na.rm=TRUE)))
    
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
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))%>% 
      dplyr::mutate(N = rowSums(across(c(S,E,I,R), na.rm=TRUE)))
    
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
                        Sc, Ec, Ic, Rc, Ssa, Esa, Isa, Rsa, Sa, Ea, Ia, Ra, M, Sm)%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc,Sm)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))%>%
      dplyr::mutate(N = rowSums(across(c(S,E,I,R,M), na.rm=TRUE)))
    
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
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc,Sm)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))%>%
      dplyr::mutate(N = rowSums(across(c(S,E,I,R,M), na.rm=TRUE)))
    
    return (list(pars = pars, init = init2, time = time, results = results))
    
  }

#############  2) Set up parameters before running #############  
#> population parameters #############
# gaur population 
N = 300 

#estimate the age structure proportion
#calf:subadult:adult ratio
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

#> disease parameters #############
#> MODEL 1 - NO infection #####
initials_m1 <- c(c = c, sa = sa, a = a )
parameters_m1 <- c(mu_b = 0.34/365, 
                mu_c = 0.27/365, 
                mu_sa = 0.15/365,
                mu_a = 0.165/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials_m1), 
                tau = 1)

#> MODEL 2 SI - Anthrax ####
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

#> MODEL 3 SEI - Bovine tuberculosis #####
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

#>  MODEL 4 SIRS - Hemorrhagic septicemia #####
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
  tau=1)

#>  MODEL 5 SEIRS - Lumpy skin disease #####
initials_m5 <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = (a-1), Ea = 0, Ia = 1, Ra = 0)
parameters_m5 <- c( 
  beta_c = 0.038, beta_sa = 0.038, beta_a = 0.038, 
  phi_c = 1/7,phi_sa = 1/7,phi_a = 1/7,
  gamma_c = 1/35, gamma_sa = 1/35,gamma_a = 1/35,
  rho_c = 0.05,rho_sa = 0.03, rho_a = 0.01, 
  omega_c = 1/180, omega_sa =  1/180, omega_a = 1/180,
  epsilon = 2e-5,
  mu_b = 0.34/365,  mu_bI = (0.34/365)*(0.9), #Ia birth rate reduce by = 10% (assumed)  
  mu_c = 0.27/365,  mu_sa = 0.15/365, mu_a = 0.165/365,
  delta_c = 1/365, delta_sa = 1/(3*365),
  N = sum(initials_m5),
  tau=1)

#> MODEL 6 SEIRMS/E - Foot and mouth disease  #####
initials_m6 <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, 
                 Sa = (a-1), Ea = 0, Ia = 1, Ra = 0)
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

#> MODEL 7 SEIRMS/E - Brucellosis  #####
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
  tau=1)

############## 3) SINGLE RUNS - ALL MODELS  ###############################
#duration of prediction
end.time = 100 * 365 

# run model 1 - 7
res_model1 <- model1(pars = parameters_m1, init = initials_m1,
                end.time = end.time)
PlotMods(res_model1)

res_model2<-model2(pars = parameters_m2, init = initials_m2,
       end.time = end.time)
PlotMods(res_model2)

res_model3<-model3(pars = parameters_m3, init = initials_m3,
                   end.time = end.time)
PlotMods(res_model3)

res_model4<-model4(pars = parameters_m4, init = initials_m4,
                   end.time = end.time)
PlotMods(res_model4)

res_model5<-model5(pars = parameters_m5, init = initials_m5,
                   end.time = end.time)
PlotMods(res_model5)

res_model6<-model6(pars = parameters_m6, init = initials_m6,
                   end.time = end.time)
res_model6$results<-res_model6$results %>% 
                    relocate(M,.after = R) 
PlotMods(res_model6)

res_model7<-model7(pars = parameters_m7, init = initials_m7,
                   end.time = end.time)
res_model7$results<-res_model7$results %>% 
                    relocate(M,.after = R) 
PlotMods(res_model7)

############## 4) MULTIPLE RUNS - ALL MODELS ################################
#### set end.time and replication times 
end.time <- 100 * 365
n_rep <- 100

# use the parameters same as a single run

#> MX RUNS MODEL 1 - NO infection #####
sim_rep_m1<-replicate(n_rep,(model1(pars = parameters_m1, init = initials_m1,
                                    end.time = end.time)))

#> MX RUNS MODEL 2 SI - Anthrax #####
sim_rep_m2<-replicate(n_rep,(model2(pars = parameters_m2, init = initials_m2,
                                    end.time = end.time)))

#> MX RUNS MODEL 3 SEI - Bovine tuberculosis #####
sim_rep_m3<-replicate(n_rep,(model3(pars = parameters_m3, init = initials_m3,
                                    end.time = end.time)))

#> MX RUNS MODEL 4 SIRS - Hemorrhagic septicemia #####
sim_rep_m4<-replicate(n_rep,(model4(pars = parameters_m4, init = initials_m4,
                                    end.time = end.time)))

#> MX RUNS MODEL 5 SEIRS - Lumpy skin disease #####
sim_rep_m5<-replicate(n_rep,(model5(pars = parameters_m5, init = initials_m5,
                                    end.time = end.time)))

#> MX RUNS MODEL 6 SEIRMS/E - Foot and mouth disease #####
sim_rep_m6<-replicate(n_rep,(model6(pars = parameters_m6, init = initials_m6,
                                    end.time = end.time)))

#> MX RUNS MODEL 7 SEIRMS/E - Brucellosis #####
sim_rep_m7<-replicate(n_rep,(model7(pars = parameters_m7, init = initials_m7,
                                    end.time = end.time)))

######## 5) PLOTTING Single and Multiple runs ######## 

# SET UP single run df&plots ########
# convert res_model to data.frame, change days -> years
s<-list(res_model1,
                res_model2,
                res_model3,
                res_model4,
                res_model5,
                res_model6,
                res_model7)
nam<-c('no_infection',
       'Anthrax',
       'bTB',
       'HS',
       'LSD',
       'FMD',
       'Brucellosis')

# arrange dataframe using melt function 
for(i in 1:length(s)){
  s[[i]]<-s[[i]]$results %>%
  mutate(time_y = time/365) %>% #convert day to year for plotting
  melt(id.vars = c('time','time_y'),
         value.name = 'value', variable.name = 'class')
  
  s[[i]]$model <- paste0(nam[[i]]) # adding name
}

  ms6<-res_model6$results%>%
    mutate(time_y = time/365) %>% #convert day to year for plotting
    melt(id.vars = c('time','time_y'),
         value.name = 'value', variable.name = 'class')
  
  ms6$model <- paste0('FMD') # adding name

str(ms6)

# save the data frame  (.rds) for working next time
for (i in 1:length(s)) {
  saveRDS(s[[i]], file = paste0("df_m",i,"_",nam[[i]],"_1run_fd.rds")) }
saveRDS(ms6,file = "df_m6_fmd_1run_fd.rds")

# SET UP mutiple run df&plots ########
# > FUNCTION for plotting model's simulations #####
pop_sim_prep <- function(x, n_rep, end.time, melt){ # x = simulation of model, e.g. sim_run_m1
  df<-list()
  for (i in 1:n_rep){
    run <- paste("run", seq(n_rep), sep="")
    names(df)[i]
    df[[i]]<- x[,i]$results[,-c(1)]
    df[[i]]$time_d <- seq(from=1,to=end.time+1,by=1)
  }
  
  df<-map2(df,run, ~cbind(.x, run = .y))   # adding number of replications to the column
  df2<- data.table::rbindlist(df)          # binding row
  
  if (melt == T) {  #option for melting the data in case we need...
    
    df3 <- melt(df2, id.vars = c('time_d','run'))
    return(df3 = data.frame(df3))
  } 
  
  else  {
    return( df2 = data.frame(df2)) 
  }
  
}

# creating the model simulation list 
sim_rep_m<-list(sim_rep_m1,
                sim_rep_m2,
                sim_rep_m3,
                sim_rep_m4,
                sim_rep_m5,
                sim_rep_m6,
                sim_rep_m7)

# creating name list for each model: sim_rep_m1 - sim_rep_m7
nam<-c('no_infection',
       'Anthrax',
       'bTB',
       'HS',
       'LSD',
       'FMD',
       'Brucellosis')

# > rearrange df and calculate total population change (%) loop--------

m<-list()

for (i in 1:length(sim_rep_m)) {
  m[[i]]<- pop_sim_prep(x = sim_rep_m[[i]], n_rep=n_rep, end.time= end.time, melt = F)
  m[[i]]<- m[[i]]%>%
    group_by(run) %>%
    mutate(Ndiff = ((last(N)-first(N))/first(N))*100)%>% #calculate change in the total population at year100, and year0
    mutate(time_y = time_d/365) %>% #convert day to year for plotting
    as.data.frame()
  m[[i]]$model <- paste0(nam[[i]])
}

m6<-pop_sim_prep(x = sim_rep_m6, n_rep=n_rep, end.time= end.time, melt = F)
m62<- m6%>%
  group_by(run) %>%
  mutate(Ndiff = ((last(N)-first(N))/first(N))*100)%>% #calculate change in the total population at year100, and year0
  mutate(time_y = time_d/365) %>% #convert day to year for plotting
  as.data.frame()
m62$model <- paste0("FMD")
str(m62)
m62<-m62%>% 
  relocate(M,.after = R) 
# relocate columns of FMD and Brucellosis models, this should be in order SEIRMN which will easier for plotting-labeling
m[[6]]<-m[[6]]%>% 
  relocate(M,.after = R) 
m[[7]]<-m[[7]]%>% 
  relocate(M,.after = R) 
# check columns before plotting: the population class should be in an order like SIRN, SIN, SIERN, if not back to the relocate
for (i in 1:length(m)){
  names(m)[i]
  print(head(m[[i]]))
}

# save the data frame  (.rds) for working next time
for (i in 1:length(m)) {
  saveRDS(m[[i]], file = paste0("df_m",i,"_",nam[[i]],"_100runs.rds")) }

# PLOTTING models #######

######## plot signle run ######
####### this can skip #######
# Load single runs .rds files and create a list
sl <- union(list.files(path = getwd(), pattern = "df_m1_no_infection_1run.rds"),
            list.files(path = getwd(), pattern = "_1run_fd.rds"))
sl
s = lapply(sl, readRDS)
str(s)
# # # #
#> p1 - plot 1 run no infection ########
p1<-ggplot(s[[1]]) + 
  geom_line(aes(x = time_y ,y = value,  color = class))  +
  labs(x="years", y= "population",
       title= 'A) No infection') +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  #  ylim(0, 1000) +
  scale_color_manual( name = "population",
                      labels = c( 'adult','subadult','calf','total' ),
                      values = c('a'='seagreen4',
                                 'sa'='firebrick',
                                 'c'='dodgerblue3',
                                 'N'='#153030'))+ 
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))

#print(p1)
#ggsave("gaur_m1_noinf_1run.png",p1,width = 22, height = 15, units = 'cm', dpi = 600)

#> p2 - plot 1 run SI Anthrax ######
s2 <- s[[2]] %>%  filter(class %in% c("S","I","N"))
p2<-ggplot(s2) + 
  geom_line(aes(x = time_y ,y = value,  color = class))  +
  labs(x="years", y= "population",
       title= 'B) Anthrax infection') +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  #  ylim(0, 1000) +
  scale_color_manual( name = "population",
                      labels = c('S','I', 'total') , #this one label is manual
                      values = c('S'='seagreen4',
                                 'I'='firebrick',
                                 'N'='#153030'))+
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))

#print(p2)
#ggsave("gaur_m2_anthrax_1run.png",p2, width = 22, height = 15, units = 'cm', dpi = 600)

#> p3 - plot 1 run SEI bTB ######
s3 <- s[[3]] %>% filter(class %in% c("S","E","I", "N"))
p3 <-ggplot(s3) + 
  geom_line(aes(x = time_y, y = value, color = class))+
  labs(x="years", y= "population",
       title= 'C) bTB infection') +
  # ylim(0,1000)+
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                      labels = c('S','E','I','total' ),
                      values = c('S'='seagreen4',
                                 'E'='darkorange2',
                                 'I'='firebrick',
                                 "N"='#153030'))+ 
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))

#print(p3)
#ggsave("gaur_m3_btb_1run.png",p3, width = 22, height = 15, units = 'cm', dpi = 600)

#> p4 plot 1 run SEI HS ######
s4 <- s[[4]] %>% filter(class %in% c("S","I","R", "N"))
p4 <-ggplot(s4) + 
  geom_line(aes(x = time_y, y = value, color = class))+
  labs(x="years", y= "population",
       title= 'D) HS infection') +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  #  ylim(0, 1000) +
  scale_color_manual( name = "population",
                      labels = c('S','I','R','total' ),
                      values = c('S'='seagreen4',
                                 'I'='firebrick',
                                 "R"='dodgerblue3',
                                 "N"='#153030'))+ 
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))

#print(p4)
#ggsave("gaur_m4_HS_1run.png",p4, width = 22, height = 15, units = 'cm', dpi = 600)

#> p5 plot 1 run SEIR LSD ######
s5<-s[[5]] |>filter(class %in% c("S","E","I","R","N"))
p5<-ggplot(s5) + 
  geom_line(aes(x = time_y, y = value, color = class))+
  labs(x="years", y= "population",
       title= 'E) LSD infection') +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  #  ylim(0, 1000) +
  scale_color_manual( name = "population",
                      labels = c('S','E','I',"R",'total' ),
                      values = c('S'='seagreen4',
                                 'E'='darkorange2',
                                 'I'='firebrick',
                                 "R"='dodgerblue3',
                                 "N"='#153030'))+ 
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))

#print(p5)
#ggsave("gaur_m5_LSD_1run.png",p5, width = 22, height = 15, units = 'cm', dpi = 600)

#> p6 plot 1 run SEIRMS/E FMD ######
s6<-s[[6]] |>filter(class %in% c("S","E","I","R","M","N"))
p6 <-ggplot(s6) + 
  geom_line(aes(x = time_y, y = value, color = class))+
  labs(x="years", y= "population",
       title= 'F) FMD infection') +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  #ylim(0, 1000) +
  scale_color_manual( name = "population",
                      labels = c('S','E','I',"R",'M','total' ),#'total change (%)'),
                      values = c('S'='seagreen4',
                                 'E'='darkorange2',
                                 'I'='firebrick',
                                 "R"='dodgerblue3',
                                 "M"='mediumorchid4',
                                 "N"='#153030'))+ 
  
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))
#print(p6)
ggsave("gaur_m6_fmd_1run_fd.png",p6, width = 22, height = 15, units = 'cm', dpi = 600)

#> p7 plot 1 run SEIRMS/E Brucellosis ######
s7<-s[[7]] |>filter(class %in% c("S","E","I","R","M","N"))
p7<-ggplot(s7) + 
  geom_line(aes(x = time_y, y = value, color = class))+
  labs(x="years", y= "population",
       title= 'G) Brucellosis infection') +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  #ylim(0, 1000) +
  scale_color_manual( name = "population",
                      labels = c('S','E','I',"R",'M','total' ),#'total change (%)'),
                      values = c('S'='seagreen4',
                                 'E'='darkorange2',
                                 'I'='firebrick',
                                 "R"='dodgerblue3',
                                 "M"='mediumorchid4',
                                 "N"='#153030'))+ 
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7))) 

#print(p7)
#ggsave("gaur_m7_bru_1run_fd.png",p7, width = 22, height = 15, units = 'cm', dpi = 600)

# multiple run plots ###### 

####### this can skip #######
# Load multiple runs .rds files and create a list
l<- list.files(path = getwd(), pattern = "_100runs.rds")
l
m = lapply(l, readRDS)
str(m)
# # # # #
#  > pl1 - plot 100 runs no infection ########
pl1<-ggplot(m[[1]]) + 
  geom_line(aes(x = time_y ,y = a,  group = run, color = 'adult'), linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = sa, group = run, color = 'subadult'),linewidth = 0.1, alpha = 0.12) +
  geom_line(aes(x = time_y, y = c,  group = run, color = 'calf' ), linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = N,  group = run, color = 'total'), linewidth = 0.1, alpha = 0.12) +
  labs(x="years", y= "population", title = 'H) No infection') +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                      labels = c( 'adult','subadult','calf','total'),
                      values = c('seagreen4',
                                 'firebrick',
                                 'dodgerblue3',
                                 '#153030'),
                      breaks = c('adult','subadult','calf','total'))+ 
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth =0.7 )))+
  
  stat_summary(m[[1]], mapping =aes( x = time_y, y = a, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m[[1]], mapping =aes( x = time_y, y = sa, group = 1), fun=mean, geom="line", colour= 'firebrick',linewidth = 0.5)+
  stat_summary(m[[1]], mapping =aes( x = time_y, y = c, group = 1), fun=mean, geom="line", colour= 'dodgerblue3',linewidth = 0.5)+
  stat_summary(m[[1]], mapping =aes( x = time_y, y = N, group = 1), fun=mean, geom="line", colour= '#153030',linewidth = 0.5)

#print(pl1)
#ggsave("gaur_m1_noinf_100runs.png",pl1,width = 22, height = 15, units = 'cm', dpi = 600)

#> pl2 - plot 100 runs SI Anthrax ######
pl2<-ggplot(m[[2]])+ 
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'), linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'), linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ), linewidth = 0.1, alpha = 0.12)+
  labs(x="years", y= "population", title= 'I) Anthrax infection') +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                      labels = c('S','I', 'total') , #this one label is manual
                      values = c('seagreen4',
                                 'firebrick',
                                 '#153030'),
                      breaks = c('S','I','total'))+
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title =element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))+
  stat_summary(m[[2]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
  stat_summary(m[[2]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m[[2]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom ="line", colour="firebrick",linewidth = 0.5)
 
#print(pl2)
#ggsave("gaur_m2_anthrax_100runs.png",pl2,width = 22, height = 15, units = 'cm', dpi = 600)

#> pl3 - plot 100 runs SEI bTB ######
pl3<-ggplot(m[[3]]) + 
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="years", y= "population", title = ('J) bTB infection')) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                      labels = c('S','E','I','total'),
                      values = c('seagreen4',
                                 'darkorange2',
                                 'firebrick',
                                 '#153030'),
                      breaks = c('S','E','I','total'))+ #blackgreen
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title =element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
  stat_summary(m[[3]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m[[3]], mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
  stat_summary(m[[3]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m[[3]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5) #blackgreen

#print(pl3)
#ggsave("gaur_m3_bTB_100runs.png",pl3,width = 22, height = 15, units = 'cm', dpi = 600)

#> pl4 - plot 100 runs SIRS HS ######
pl4 <-ggplot(m[[4]]) + 
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="years", y= "population", title = ('K) HS infection')) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                      labels = c('S','I','R','total'),
                      values = c('seagreen4',
                                 'firebrick',
                                 'dodgerblue3',
                                 '#153030'),
                      breaks = c('S','I','R','total'))+ #blackgreen
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title =element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
  stat_summary(m[[4]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
  stat_summary(m[[4]], mapping = aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m[[4]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m[[4]], mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)
 
#print(pl4)
#ggsave("gaur_m4_HS_100runs.png",pl4,width = 22, height = 15, units = 'cm', dpi = 600)

#> pl5 - plot 100 runs SEIRS LSD ######
pl5 <-ggplot(m[[5]]) + 
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="years", y= "population", title = ('L) LSD infection')) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                      labels = c('S','E','I','R','total'),
                      values = c('seagreen4',
                                 'darkorange2',
                                 'firebrick',
                                 'dodgerblue3',
                                 '#153030'),
                      breaks = c('S','E','I','R','total'))+ #blackgreen
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title =element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
  stat_summary(m[[5]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+
  stat_summary(m[[5]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m[[5]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
  stat_summary(m[[5]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m[[5]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)#blackgreen

#print(pl5)
#ggsave("gaur_m5_LSD_100runs.png",pl5,width = 22, height = 15, units = 'cm', dpi = 600)

#> pl6 - plot 100 runs SEIRMS/E FMD ######
pl6 <-ggplot(m[[6]]) + 
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="years", y= "population", title = ('M) FMD infection')) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                      labels = c('S','E','I','R','M','total'),
                      values = c('seagreen4',
                                 'darkorange2',
                                 'firebrick',
                                 'dodgerblue3',
                                 'mediumorchid4',
                                 '#153030'),
                      breaks = c('S','E','I','R','M','total'))+ #blackgreen
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title =element_text(size = 11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size = 11))+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
  stat_summary(m[[6]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
  stat_summary(m[[6]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m[[6]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
  stat_summary(m[[6]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m[[6]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
  stat_summary(m[[6]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)
  
#print(pl6)
#ggsave("gaur_m6_FMD_100runs_fd.png",pl6,width = 22, height = 15, units = 'cm', dpi = 600)

#> pl7 - plot 100 runs SEIRMS/E Brucellosis ######
pl7 <-ggplot(m[[7]]) + 
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="years", y= "population", title = ('N) Brucellosis infection')) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                      labels = c('S','E','I','R','M','total'),
                      values = c('seagreen4',
                                 'darkorange2',
                                 'firebrick',
                                 'dodgerblue3',
                                 'mediumorchid4',
                                 '#153030'),
                      breaks = c('S','E','I','R','M','total'))+ #blackgreen
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title =element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
  stat_summary(m[[7]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
  stat_summary(m[[7]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m[[7]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
  stat_summary(m[[7]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m[[7]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
  stat_summary(m[[7]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)
  
#print(pl7)
#ggsave("gaur_m7_Brucellosis_100runs.png",pl7,width = 22, height = 15, units = 'cm', dpi = 600)

#print(p7)
#ggsave("gaur_m7_bru_1run.png",p7, width = 22, height = 15, units = 'cm', dpi = 600)

# Combine plots using Patchwork #######
library(patchwork)
library(grid)
pp1 <- p1/p2/p3/p4/p5/p6/p7 & plot_annotation(title = 'Frequeny dependent, 1 simulation') & theme(plot.title = element_text(hjust = 0.1))
pp2 <- pl1/pl2/pl3/pl4/pl5/pl6/pl7 & plot_annotation (title = "100 simulations")  & theme(plot.title = element_text(hjust = 0.1))
pp3 <- pp1|pp2 #& plot_annotation(tag_levels = 'A', tag_suffix = ')')

png("patchwork_fd_1run.png",width = 20, height = 40, units = 'cm', res = 600)
print(pp1)
dev.off()

# > zoom in FMD in 10 years
z<-s6 |>
  filter(between(time_y,10,20))|>
  mutate(year=c("10 - 20 years"))
head(z)

z2<-s6 |> filter(between(time_y,90,100))|>
  mutate(year=c("90 - 100 years"))
head(z2)

z3<-rbind(z,z2)
head(z3)

p6z<-ggplot(z3) + 
  geom_line(aes(x = time_y, y = value, color = class), linewidth = 1)+
  labs(x="years", y= "population",
       title= 'Gaur population with FMD infection') +
  facet_wrap(year~., scale = "free_x") +
  
  scale_x_continuous(breaks=seq(0, (365*10), by = 2))+
  
  scale_color_manual( name = "population",
                      labels = c('S','E','I',"R",'M','total' ),
                      values = c('S'='seagreen4',
                                 'E'='darkorange2',
                                 'I'='firebrick',
                                 "R"='dodgerblue3',
                                 "M"='mediumorchid4',
                                 "N"='#153030'))+ 
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))

print(p6z)
#ref:https://stackoverflow.com/questions/56064042/using-facet-tags-and-strip-labels-together-in-ggplot2
tag_facet2 <- function(p, close = ")", tag_pool = LETTERS, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 1, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

p6z2<-tag_facet2(p6z)
print(p6z2)
ggsave("gaur_fmd_1sim_zoom10y.png",p6z2, width = 25, height = 15, units = 'cm', dpi = 600)

# Plot the average population change ##########
# Prepare the dataframe for boxplot 
# the % of population change in 100 for all the models

# arrange dataframe using melt
mx<-list()
for (i in 1:length(m)){
  
  mx[[i]] <- melt(m[[i]], id.vars = c('time_y','time_d','run','model',"Ndiff"),
                  value.name = 'value', variable.name = 'class')
}

####### summerize basic stat ####### 
# single run
for(i in 1:length(s)){
  print(s[[i]] |> group_by(class) |>
          dplyr::summarise(Mean = mean(value),
                           Max=max(value),
                           Min=min(value)))
}

# multiple runs
for (i in 1:length(mx)){
  print(mx[[i]] |> group_by(class) |>
          dplyr::summarise(#Median = median(value), 
            #Mean = mean(value),
            Max=max(value),
            Min=min(value)))
}

mx2 <-data.table::rbindlist(mx)

str(mx2)

#Prepare the dataframe for boxplot ########
#the % of population change in 100 for all the models
dft<-mx2 %>% 
  dplyr::select(Ndiff,run,model)%>% 
  drop_na()%>%
  distinct()

table(dft$Ndiff)
write_csv(dft,"df_ndiff.csv")

#load df back
#dft <- read_csv("df_ndiff.csv")

#calculating mean,median
dft |>
  group_by(model)|>
  summarise(Median = median(Ndiff),
            Mean = mean(Ndiff))

# Boxplot
box<-dft %>%
  ggplot(aes(x= reorder(model,-Ndiff), y = Ndiff,fill=model)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.6, alpha=0.4) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)) +
  ggtitle("Gaur population change by infectious disease models in 100 years") +
  xlab("") +
  ylab("Population change (%)")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))


print(box)
ggsave("gaur_ndiff_box1_100sim.png",box,width = 22, height = 15, units = 'cm', dpi = 600)

#reorder by population change max-min
dft$model <- recode_factor(dft$model, no_infection = "no infection" )
dft$model2 <- factor(dft$model, 
                    levels = c("no infection", 
                               "HS",
                               "LSD",
                               "Anthrax",
                               "FMD",
                               "bTB",
                               "Brucellosis"))
plt<-dft%>%ggbetweenstats(
  x=model2,
  y=Ndiff,
  k=0,
  plot.type = "boxviolin",
  pairwise.comparisons=F,
  bf.message = F,
  results.subtitle = FALSE,
  centrality.point.args = list(size = 2, color = "darkred"),
  title= "Gaur population change by infectious disease models in 100 years",
  xlab = "",
  ylab = "Population change (%)",
  package = "ggsci",
  palette = "default_jco")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))

print(plt)
ggsave("gaur_ndiff_boxviolin_max-min2_100sim.png",plt,width = 20, height = 15, units = 'cm', dpi = 600)

#plotting histogram N==0  for brucellosis #######
m7<-readRDS("df_m7_Brucellosis_100runs.rds")
mx7 <- melt(m7, id.vars = c('time_y','time_d','run','model',"Ndiff"),
                value.name = 'value', variable.name = 'class')
head(mx7)
df.agg <- aggregate(time_y ~ run + value +class, mx7, min)
head(df.agg)

df.ag <- (df.agg[df.agg$value==0,c('class','time_y','value','run')])
head(df.ag)
str(df.ag)
df.ag2 <- df.ag|>filter(class %in% c("N","I") )
head(df.ag2)
str(df.ag2)

# plotting N extinction 
pn<- df.ag|>filter(class %in% c("N"))|>
  
  ggplot(aes(x=time_y,fill=run)) + 
  geom_histogram(binwidth =2,color="black",linewidth=0.2,alpha=0.7)+
  scale_y_continuous(limits = c(0,3), breaks=seq(0, 3, by = 2)) +
  scale_x_continuous(limits = c(0,105),breaks=seq(0, 105, by = 20))+
  ggtitle("Extinction events of the gaur population from brucellosis") +
  xlab("years") +
  theme( plot.title = element_text(size = 9),
         axis.title.x = element_text(size = 7),
         axis.title.y = element_text(size = 7),
         legend.title=element_text(size=7),
         legend.text = element_text(size = 7),
         axis.text=element_text(size=7))+
  scale_fill_viridis_d(option="mako",direction = -1)+
  #scale_fill_brewer(palette = "BrBG")+
  theme(legend.position = "none")

print(pn)  

tiff("extinction_bru_n.png",width = 12, height = 10, units = 'cm', res = 600)
print(pn)
dev.off()

# plotting I&N not separate 
p<-ggplot(df.ag2,aes(x=time_y,fill=class)) + 
  geom_histogram(binwidth =2,color="black",linewidth=0.2,alpha=0.7)+
  #facet_wrap(class~., ncol = 2)+
  scale_y_continuous(limits = c(0,8), breaks=seq(0, 8, by = 2)) +
  scale_x_continuous(limits = c(0,105),breaks=seq(0, 105, by = 20))+
  ggtitle("Extinction events of the gaur population from brucellosis") +
  xlab("years") +
  theme( plot.title = element_text(size = 8),
         axis.title.x = element_text(size = 7),
         axis.title.y = element_text(size = 7),
         legend.title=element_text(size=7),
         legend.text = element_text(size = 7),
         axis.text=element_text(size=7))+
  #theme(legend.position = "none")+
  #scale_fill_brewer(palette = "BrBG")
  scale_fill_viridis_d(option="mako",direction = -1,
                       name = "population", labels = c("I","total"))+
  theme(legend.key.size = unit(0.3, 'cm'))+
  #scale_fill_discrete(name = "population", labels = c("I","total"))
  theme(legend.position=c(.91, .9),legend.box.background = element_rect(colour = "black"))
print(p)  

tiff("extinction_bru_ni.png",width = 12, height = 10, units = 'cm', res = 600)
print(p)
dev.off()

# plotting I&N separately using facet
#classified by simulation
df.ag2$class <- factor(df.ag2$class, levels = c("I", "N"), 
                       labels = c("I", "total"))

p2<-ggplot(df.ag2,aes(x=time_y,fill=run)) + 
  geom_histogram(binwidth =2,color="black",linewidth=0.2,alpha=0.7)+
  facet_wrap(class~., ncol = 2)+
  scale_y_continuous(limits = c(0,8), breaks=seq(0, 8, by = 2)) +
  scale_x_continuous(limits = c(0,105),breaks=seq(0, 105, by = 20))+
  ggtitle("Extinction events of the gaur population from brucellosis") +
  xlab("years") +
  theme( plot.title = element_text(size = 9),
         axis.title.x = element_text(size = 7),
         axis.title.y = element_text(size = 7),
         legend.title=element_text(size=7),
         legend.text = element_text(size = 7),
         axis.text=element_text(size=7))+
  theme(legend.position = "none")+
  #scale_fill_brewer(palette = "BrBG")
  scale_fill_viridis_d(option="mako",direction = -1)

print(p2)  

tiff("extinction_bru_facet.png",width = 12, height = 9, units = 'cm', res = 600)
print(p2)
dev.off()

