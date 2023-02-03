############## NOTES ##################
# Extra codes for running and plotting 
# - an option for modeling both frequency disease (FD) and density disease (DD) transmissions;    

# Working on:  R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
# RStudio 2022.07.1+554 "Spotted Wakerobin" Release (7872775ebddc40635780ca1ed238934c3345c5de, 2022-07-22) for macOS
# Mozilla/5.0 (Macintosh; Intel Mac OS X 13_0_1)

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

############### 1) Set up the infectious disease model's function ############### 
# > MODEL 1 - no infection   ############### 
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

#> Frequency dependent transmission  ############### 

#  MODEL 2 SI - Anthrax 
model2fd = 
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
# MODEL 3 SEI - Bovine tuberculosis

model3fd =
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
        rate[24] <- rho_a * gamma_a * Ia
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

# MODEL 4 SIRS - Hemorrhagic septicemia
model4fd =
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

# MODEL 5 SEIRS - Lumpy skin disease
model5fd =
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

# MODEL 6 SEIRMS/E - Foot and mouth disease
model6fd =
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

# MODEL 7 SEIRMS/E - Bovine brucellosis
model7fd =
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

#> Density dependent transmission ############### 
# MODEL 2 SI - Anthrax
model2dd = 
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
        rate[2] <- beta_c * Sc * (Ic+Isa+Ia)
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
        rate[7] <- beta_sa * Ssa * (Ic+Isa+Ia)
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
        rate[12] <- beta_a * Sa * (Ic+Isa+Ia)
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
# MODEL 3 SEI - Bovine tuberculosis
model3dd =
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
        rate[2] <- beta_c * Sc * (Ic+Isa+Ia)
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
        rate[12] <- beta_sa * Ssa * (Ic+Isa+Ia)
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
        rate[22] <- beta_a * Sa * (Ic+Isa+Ia)
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

# MODEL 4 SIRS - Hemorrhagic septicemia
model4dd =
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
        rate[2] <- beta_c * Sc * (Ic+Isa+Ia)
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
        rate[13] <- beta_sa * Ssa * (Ic+Isa+Ia)
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
        rate[24] <- beta_a * Sa * (Ic+Isa+Ia)
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

# MODEL 5 SEIRS - Lumpy skin disease
model5dd =
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
        rate[3] <- beta_c * Sc * (Ic+Isa+Ia)
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
        rate[17] <- beta_sa * Ssa * (Ic+Isa+Ia)
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
        rate[32] <- beta_a * Sa * (Ic+Isa+Ia)
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

# MODEL 6 SEIRMS/E - Foot and mouth disease
model6dd = 
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
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc,Sm)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))%>%
      dplyr::mutate(N = rowSums(across(c(S,E,I,R,M), na.rm=TRUE)))
    
    return (list(pars = pars, init = init2, time = time, results = results))
    
  }

# MODEL 7 SEIRMS/E - Bovine brucellosis
model7dd=
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
                        Sc, Ec, Ic, Rc, Ssa, Esa, Isa, Rsa, Sa, Ea, Ia, Ra,  M, Sm )%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc,Sm)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))%>%
      dplyr::mutate(N = rowSums(across(c(S,E,I,R,M), na.rm=TRUE)))
    
    return (list(pars = pars, init = init2, time = time, results = results))
    
  }

#############  2) Set up parameters before running #############  
#> set end.time and number of replications (n_rep)  #############
end.time <- 100 * 365
n_rep <- 100

#> population parameters #############
# gaur population   #########
N = 300 
#estimate the age structure proportion
#calf:subadult:adult ratio
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3, 0) 
a = round((N/rat)*1.5, 0) 

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
# FD
parameters_m7fd <- c( 
  beta_c = 2/365, beta_sa = 2/365, beta_a = 2/365, #fequency: bata = 2
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

# DD
parameters_m7dd <- c( 
  beta_c = 0.002/365, beta_sa = 0.002/365, beta_a = 0.002/365, #density: bata = 0.002
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

res_model1<- model1(pars = parameters_m1, init = initials_m1,
                       end.time = end.time)
#PlotMods(res_model1)

# FD: run model 2 - 7 ##########
res_model2fd <-model2fd(pars = parameters_m2, init = initials_m2,
                   end.time = end.time)
#PlotMods(res_model2fd)

res_model3fd<-model3fd(pars = parameters_m3, init = initials_m3,
                   end.time = end.time)
#PlotMods(res_model3fd)

res_model4fd<-model4fd(pars = parameters_m4, init = initials_m4,
                   end.time = end.time)
#PlotMods(res_model4fd)

res_model5fd<-model5fd(pars = parameters_m5, init = initials_m5,
                   end.time = end.time)
#PlotMods(res_model5fd)

res_model6fd<-model6fd(pars = parameters_m6, init = initials_m6,
                   end.time = end.time)
res_model6fd$results<-res_model6fd$results %>% 
  relocate(M,.after = R) 
#PlotMods(res_model6fd)

res_model7fd<-model7fd(pars = parameters_m7fd, init = initials_m7,
                   end.time = end.time)
res_model7fd$results<-res_model7fd$results %>% 
  relocate(M,.after = R) 
#PlotMods(res_model7fd)

# DD : run model 2 - 7 ##########
res_model2dd <-model2dd(pars = parameters_m2, init = initials_m2,
                        end.time = end.time)
#PlotMods(res_model2dd)

res_model3dd<-model3dd(pars = parameters_m3, init = initials_m3,
                       end.time = end.time)
#PlotMods(res_model3dd)

res_model4dd<-model4dd(pars = parameters_m4, init = initials_m4,
                       end.time = end.time)
#PlotMods(res_model4dd)

res_model5dd<-model5dd(pars = parameters_m5, init = initials_m5,
                       end.time = end.time)
#PlotMods(res_model5dd)

res_model6dd<-model6dd(pars = parameters_m6, init = initials_m6,
                       end.time = end.time)
res_model6dd$results<-res_model6fd$results %>% 
  relocate(M,.after = R) 
#PlotMods(res_model6dd)

res_model7dd<-model7dd(pars = parameters_m7dd, init = initials_m7,
                       end.time = end.time)
res_model7dd$results<-res_model7dd$results %>% 
  relocate(M,.after = R) 
#PlotMods(res_model7dd)

############## 4) MULTIPLE RUNS - ALL MODELS ################################
# use the same parameters as the single run
#> MX RUNS MODEL 1 - NO infection #####
sim_rep_m1<-replicate(n_rep,(model1(pars = parameters_m1, init = initials_m1,
                                    end.time = end.time)))

#> MX RUNS MODEL 2 SI - Anthrax #####
sim_rep_m2fd<-replicate(n_rep,(model2fd(pars = parameters_m2, init = initials_m2,
                                    end.time = end.time)))
sim_rep_m2dd<-replicate(n_rep,(model2dd(pars = parameters_m2, init = initials_m2,
                                        end.time = end.time)))
#> MX RUNS MODEL 3 SEI - Bovine tuberculosis #####
sim_rep_m3fd<-replicate(n_rep,(model3fd(pars = parameters_m3, init = initials_m3,
                                    end.time = end.time)))
sim_rep_m3dd<-replicate(n_rep,(model3dd(pars = parameters_m3, init = initials_m3,
                                        end.time = end.time)))

#> MX RUNS MODEL 4 SIRS - Hemorrhagic septicemia #####
sim_rep_m4fd<-replicate(n_rep,(model4fd(pars = parameters_m4, init = initials_m4,
                                    end.time = end.time)))
sim_rep_m4dd<-replicate(n_rep,(model4fd(pars = parameters_m4, init = initials_m4,
                                        end.time = end.time)))
#> MX RUNS MODEL 5 SEIRS - Lumpy skin disease #####
sim_rep_m5fd<-replicate(n_rep,(model5fd(pars = parameters_m5, init = initials_m5,
                                    end.time = end.time)))
sim_rep_m5dd<-replicate(n_rep,(model5dd(pars = parameters_m5, init = initials_m5,
                                        end.time = end.time)))
#> MX RUNS MODEL 6 SEIRMS/E - Foot and mouth disease #####
sim_rep_m6fd<-replicate(n_rep,(model6fd(pars = parameters_m6, init = initials_m6,
                                    end.time = end.time)))
sim_rep_m6dd<-replicate(n_rep,(model6fd(pars = parameters_m6, init = initials_m6,
                                        end.time = end.time)))
#> MX RUNS MODEL 7 SEIRMS/E - Brucellosis #####
sim_rep_m7fd<-replicate(n_rep,(model7fd(pars = parameters_m7fd, init = initials_m7,
                                    end.time = end.time)))
sim_rep_m7dd<-replicate(n_rep,(model7dd(pars = parameters_m7dd, init = initials_m7,
                                        end.time = end.time)))

######## 5) PLOTTING Single and Multiple runs ######## 

# SET UP single run df&plots ########
# convert res_model to data.frame, change days -> years
s_fd<-list(res_model1,
        res_model2fd,
        res_model3fd,
        res_model4fd,
        res_model5fd,
        res_model6fd,
        res_model7fd)
s_dd<-list(res_model1,
           res_model2dd,
           res_model3dd,
           res_model4dd,
           res_model5dd,
           res_model6dd,
           res_model7dd)

nam<-c('no_infection',
       'Anthrax',
       'bTB',
       'HS',
       'LSD',
       'FMD',
       'Brucellosis')

# arrange dataframe using melt function 
for(i in 1:length(s_fd)){
  s_fd[[i]]<-s_fd[[i]]$results %>%
    mutate(time_y = time/365) %>% #convert day to year for plotting
    melt(id.vars = c('time','time_y'),
         value.name = 'value', variable.name = 'class')
  
  s_fd[[i]]$model <- paste0(nam[[i]]) # adding name
#}  
#for(i in 1:length(s_dd)){
  s_dd[[i]]<-s_dd[[i]]$results %>%
    mutate(time_y = time/365) %>% #convert day to year for plotting
    melt(id.vars = c('time','time_y'),
         value.name = 'value', variable.name = 'class')
  
  s_dd[[i]]$model <- paste0(nam[[i]]) # adding name
}
s_fd
s_dd
# save the data frame  (.rds) for working next time
for (i in 1:length(s_fd)) {
  saveRDS(s_fd[[i]], file = paste0("df_m",i,"_",nam[[i]],"_1run_fd.rds")) }
# save the data frame  (.rds) for working next time
for (i in 1:length(s_dd)) {
  saveRDS(s_dd[[i]], file = paste0("df_m",i,"_",nam[[i]],"_1run_dd.rds")) }

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
sim_rep_m_fd<-list(sim_rep_m1,
                sim_rep_m2fd,
                sim_rep_m3fd,
                sim_rep_m4fd,
                sim_rep_m5fd,
                sim_rep_m6fd,
                sim_rep_m7fd)
sim_rep_m_dd<-list(sim_rep_m1,
                   sim_rep_m2dd,
                   sim_rep_m3dd,
                   sim_rep_m4dd,
                   sim_rep_m5dd,
                   sim_rep_m6dd,
                   sim_rep_m7dd)

# > rearrange df and calculate total population change (%) loop--------
mfd<-list()
mdd<-list()
for (i in 1:length(sim_rep_m_fd)) {
  mfd[[i]]<- pop_sim_prep(x = sim_rep_m_fd[[i]], n_rep=n_rep, end.time= end.time, melt = F)
  mfd[[i]]<- mfd[[i]]%>%
    group_by(run) %>%
    mutate(Ndiff = ((last(N)-first(N))/first(N))*100)%>% #calculate change in the total population at year100, and year0
    mutate(time_y = time_d/365) %>% #convert day to year for plotting
    as.data.frame()
  mfd[[i]]$model <- paste0(nam[[i]])
}

for (i in 1:length(sim_rep_m_dd)) {
  mdd[[i]]<- pop_sim_prep(x = sim_rep_m_dd[[i]], n_rep=n_rep, end.time= end.time, melt = F)
  mdd[[i]]<- mdd[[i]]%>%
    group_by(run) %>%
    mutate(Ndiff = ((last(N)-first(N))/first(N))*100)%>% #calculate change in the total population at year100, and year0
    mutate(time_y = time_d/365) %>% #convert day to year for plotting
    as.data.frame()
  mdd[[i]]$model <- paste0(nam[[i]])
}

mfd[[6]]<-mfd[[6]]%>% 
  relocate(M,.after = R) 
mfd[[7]]<-mfd[[7]]%>% 
  relocate(M,.after = R) 
mdd[[6]]<-mdd[[6]]%>% 
  relocate(M,.after = R) 
mdd[[7]]<-mdd[[7]]%>% 
  relocate(M,.after = R) 
# check columns before plotting: the population class should be in an order like SIRN, SIN, SIERN, if not back to the relocate
for (i in 1:length(mfd)){
  names(mfd)[i]
  print(head(mfd[[i]]))
}
for (i in 1:length(mdd)){
  names(mdd)[i]
  print(head(mdd[[i]]))
}
# save the data frame  (.rds) for working next time
for (i in 1:length(mfd)) {
  saveRDS(mfd[[i]], file = paste0("df_m",i,"_",nam[[i]],"_mruns_fd.rds")) }

for (i in 1:length(mdd)) {
  saveRDS(mdd[[i]], file = paste0("df_m",i,"_",nam[[i]],"_mruns_dd.rds")) }

# # # #

# PLOTTING models #######
####### this can skip #######
# Load single runs .rds files and create a list
slfd <- list.files(path = getwd(), pattern = "_1run_fd.rds")
slfd
sldd <- list.files(path = getwd(), pattern = "_1run_dd.rds")
sldd

sfd = lapply(slfd, readRDS)
sdd = lapply(sldd, readRDS)

str(sfd)
str(sdd)

####### this can skip #######
# Load multiple runs .rds files and create a list
lfd <- list.files(path = getwd(), pattern = "mruns_fd.rds")
lfd

ldd <- list.files(path = getwd(), pattern =  "_mruns_dd.rds")
ldd

mfd = lapply(lfd, readRDS)
mdd = lapply(ldd, readRDS)

str(mfd)
str(mdd)

# # # #
# name for plotting
nam2<-list('No',
           'Anthrax',
           'bTB',
           'HS',
           'LSD',
           'FMD',
           'Brucellosis')
#1) plot single run ######
# FD plots  ######
pfd<-list()

for (i in 1:length(sfd)) {
  # Condition
  model <- sfd[[i]]$model
  
  #> p1 - plot 1 run no infection ########
  if (any(model == "no_infection")){
    pfd[[i]]<-ggplot(sfd[[i]]) + 
      geom_line(aes(x = time_y ,y = value,  color = class))  +
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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

    #ggsave("gaur_m1_noinf_1run.png",p1,width = 22, height = 15, units = 'cm', dpi = 600) 
  } 
  else if (any(model == "Anthrax")) {
    #> p2 - plot 1 run SI Anthrax ######
    pfd[[i]]<-sfd[[i]] %>%  filter(class %in% c("S","I","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y ,y = value,  color = class))  +
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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
    
    #ggsave("gaur_m2_anthrax_1run_fd.png",p2, width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if(any(model == "bTB")) {
    #> p3 - plot 1 run SEI bTB ######
    pfd[[i]] <- sfd[[i]] %>% filter(class %in% c("S","E","I", "N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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
    #ggsave("gaur_m3_btb_1run_fd.png",p3, width = 22, height = 15, units = 'cm', dpi = 600)
  } 
  else if (any(model == "HS")) {
    #> p4 plot 1 run SIRS HS ######
    pfd[[i]] <- sfd[[i]] %>% filter(class %in% c("S","I","R","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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
      #ggsave("gaur_m4_HS_1run_fd.png",p4, width = 22, height = 15, units = 'cm', dpi = 600)
    }

  else if  (any(model =="LSD")) {
    #> p5 plot 1 run SEIR LSD ######
    pfd[[i]]<- sfd[[i]] %>% filter(class %in% c("S","E","I","R","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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
  
    #ggsave("gaur_m5_LSD_1run_fd.png",p5, width = 22, height = 15, units = 'cm', dpi = 600)
  } 
  else if  (any(model == "FMD")) {
    #> p6 plot 1 run SEIRMS/E FMD ######
    pfd[[i]]<-sfd[[i]] |>filter(class %in% c("S","E","I","R","M","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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
    #ggsave("gaur_m6_fmd_1run_fd.png",p6, width = 22, height = 15, units = 'cm', dpi = 600)
  
    } else {
    #> p7 plot 1 run SEIRMS/E Brucellosis ######
    pfd[[i]]<-sfd[[i]]%>% filter(class %in% c("S","E","I","R","M","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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
    #ggsave("gaur_m7_bru_1run_fd.png",p7, width = 22, height = 15, units = 'cm', dpi = 600)
  }
}

print(pfd)

# DD plots   ######
pdd<-list()

for (i in 1:length(sdd)) {
  # Condition
  model <- sdd[[i]]$model
  
  #> p1 - plot 1 run no infection ########
  if (any(model == "no_infection")){
    pdd[[i]]<-ggplot(sdd[[i]]) + 
      geom_line(aes(x = time_y ,y = value,  color = class))  +
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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
  } 
  else if (any(model == "Anthrax")) {
    #> p2 - plot 1 run SI Anthrax ######
    pdd[[i]]<-sdd[[i]] %>%  filter(class %in% c("S","I","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y ,y = value,  color = class))  +
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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

    #ggsave("gaur_m2_anthrax_1run_dd.png",p2, width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if(any(model == "bTB")) {
    #> p3 - plot 1 run SEI bTB ######
    pdd[[i]] <- sdd[[i]] %>% filter(class %in% c("S","E","I", "N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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
  
    #ggsave("gaur_m3_btb_1run_dd.png",p3, width = 22, height = 15, units = 'cm', dpi = 600)
  } 
  else if (any(model == "HS")) {
    #> p4 plot 1 run SIRS HS ######
    pdd[[i]] <- sdd[[i]] %>% filter(class %in% c("S","I","R","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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
      #ggsave("gaur_m4_HS_1run_dd.png",p4, width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if  (any(model =="LSD")) {
    #> p5 plot 1 run SEIR LSD ######
    pdd[[i]]<- sdd[[i]] %>% filter(class %in% c("S","E","I","R","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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

    #ggsave("gaur_m5_LSD_1run_dd.png",p5, width = 22, height = 15, units = 'cm', dpi = 600)
  } 
  else if  (any(model == "FMD")) {
    #> p6 plot 1 run SEIRMS/E FMD ######
    pdd[[i]]<-sdd[[i]] |>filter(class %in% c("S","E","I","R","M","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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
    #ggsave("gaur_m6_fmd_1run_dd.png",p6, width = 22, height = 15, units = 'cm', dpi = 600)
  } else {
    #> p7 plot 1 run SEIRMS/E Brucellosis ######
    pdd[[i]]<-sdd[[i]]%>% filter(class %in% c("S","E","I","R","M","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="years", y= "population",
           title= paste0(LETTERS[[i]],") ", nam2[[i]], " infection")) +
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
    #ggsave("gaur_m7_bru_1run_dd.png",p7, width = 22, height = 15, units = 'cm', dpi = 600)
  }
}

print(pdd)
# # # # #

#2) multiple run plots ###### 
# FD plots  ######
plfd<-list()
for (i in 1:length(mfd)) {
  # Condition
  model <- mfd[[i]]$model
  
  #  > pl1 - plot 100 runs no infection ########
  if (any(mfd[[i]]$model == "no_infection")) {
    plfd[[i]]<-ggplot(mfd[[i]]) + 
      geom_line(aes(x = time_y ,y = a,  group = run, color = 'adult'), linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = sa, group = run, color = 'subadult'),linewidth = 0.1, alpha = 0.12) +
      geom_line(aes(x = time_y, y = c,  group = run, color = 'calf' ), linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = N,  group = run, color = 'total'), linewidth = 0.1, alpha = 0.12) +
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "population",
                          labels = c('adult','subadult','calf','total'),
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
             axis.text=element_text(size=11)) +
      
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth =0.7 )))+
      
      stat_summary(mfd[[i]], mapping =aes( x = time_y, y = a, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping =aes( x = time_y, y = sa, group = 1), fun=mean, geom="line", colour= 'firebrick',linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping =aes( x = time_y, y = c, group = 1), fun=mean, geom="line", colour= 'dodgerblue3',linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping =aes( x = time_y, y = N, group = 1), fun=mean, geom="line", colour= '#153030',linewidth = 0.5)
    
    #print(pl1)
    #ggsave("gaur_m1_noinf_100runs.png",pl1,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if (any(model == "Anthrax")) {       
    #> pl2 - plot 100 runs SI Anthrax ######
    plfd[[i]] <-ggplot(mfd[[i]])+ 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'), linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'), linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ), linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom ="line", colour="firebrick",linewidth = 0.5)
    
    #print(pl2)
    #ggsave("gaur_m2_anthrax_100runs_fd.png",pl2,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if(any(model == "bTB")) {
    #> pl3 - plot 100 runs SEI bTB ######
    plfd[[i]] <-ggplot(mfd[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mfd[[i]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5) #blackgreen
  
    #ggsave("gaur_m3_bTB_100runs_fd.png",pl3,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if(any(model == "HS")) {
    #> pl4 - plot 100 runs SIRS HS ######
    plfd[[i]] <-ggplot(mfd[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mfd[[4]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(mfd[[4]], mapping = aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mfd[[4]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(mfd[[4]], mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)
    
    #print(pl4)
    #ggsave("gaur_m4_HS_100runs_fd.png",pl4,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if(any(model == "LSD")){
    #> pl5 - plot 100 runs SEIRS LSD ######
    plfd[[i]] <-ggplot(mfd[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)#blackgreen
    
    #print(pl5)
    #ggsave("gaur_m5_LSD_100runs_fd.png",pl5,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if(any(model == "FMD"))  {
    #> pl6 - plot 100 runs SEIRMS/E FMD ######
    plfd[[i]] <-ggplot(mfd[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)
    
    #print(pl6)
    #ggsave("gaur_m6_FMD_100runs_fd.png",pl6,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else {
    #> pl7 - plot 100 runs SEIRMS/E Brucellosis ######
    plfd[[i]] <-ggplot(mfd[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
      stat_summary(mfd[[i]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)
    #ggsave("gaur_m7_Brucellosis_100runs_fd.png",pl7,width = 22, height = 15, units = 'cm', dpi = 600)
  }
}   

print(plfd)

# DD plots  ######
pldd<-list()

for (i in 1:length(mdd)) {
  # Condition
  model <- mdd[[i]]$model
  
  #  > pl1 - plot 100 runs no infection ########
  if (any(mdd[[i]]$model == "no_infection")) {
    pldd[[i]]<-ggplot(mdd[[i]]) + 
      geom_line(aes(x = time_y ,y = a,  group = run, color = 'adult'), linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = sa, group = run, color = 'subadult'),linewidth = 0.1, alpha = 0.12) +
      geom_line(aes(x = time_y, y = c,  group = run, color = 'calf' ), linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = N,  group = run, color = 'total'), linewidth = 0.1, alpha = 0.12) +
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      
      stat_summary(mdd[[i]], mapping =aes( x = time_y, y = a, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping =aes( x = time_y, y = sa, group = 1), fun=mean, geom="line", colour= 'firebrick',linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping =aes( x = time_y, y = c, group = 1), fun=mean, geom="line", colour= 'dodgerblue3',linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping =aes( x = time_y, y = N, group = 1), fun=mean, geom="line", colour= '#153030',linewidth = 0.5)
    
    #ggsave("gaur_m1_noinf_100runs.png",pl1,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if (any(model == "Anthrax")) {       
    #> pl2 - plot 100 runs SI Anthrax ######
    pldd[[i]] <-ggplot(mdd[[i]])+ 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'), linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'), linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ), linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom ="line", colour="firebrick",linewidth = 0.5)
    
    #ggsave("gaur_m2_anthrax_100runs_dd.png",pl2,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if(any(model == "bTB")) {
    #> pl3 - plot 100 runs SEI bTB ######
    pldd[[i]] <-ggplot(mdd[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mdd[[i]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5) #blackgreen
    
    #print(pl3)
    #ggsave("gaur_m3_bTB_100runs.png",pl3,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if(any(model == "HS")) {
    #> pl4 - plot 100 runs SIRS HS ######
    pldd[[i]] <-ggplot(mdd[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mdd[[4]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(mdd[[4]], mapping = aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mdd[[4]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(mdd[[4]], mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)
    
      #ggsave("gaur_m4_HS_100runs_dd.png",pl4,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if(any(model == "LSD")){
    #> pl5 - plot 100 runs SEIRS LSD ######
    pldd[[i]] <-ggplot(mdd[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)#blackgreen
    
    #print(pl5)
    #ggsave("gaur_m5_LSD_100runs_dd.png",pl5,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else if(any(model == "FMD"))  {
    #> pl6 - plot 100 runs SEIRMS/E FMD ######
    pldd[[i]] <-ggplot(mdd[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title= paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)
    
    #print(pl6)
    #ggsave("gaur_m6_FMD_100runs_fd.png",pl6,width = 22, height = 15, units = 'cm', dpi = 600)
  }
  else {
    #> pl7 - plot 100 runs SEIRMS/E Brucellosis ######
    pldd[[i]] <-ggplot(mdd[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", title=paste0(LETTERS[8:14][i],") ", nam2[[i]], " infection")) +
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
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
      stat_summary(mdd[[i]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)

    #ggsave("gaur_m7_Brucellosis_100runs_dd.png",pl7,width = 22, height = 15, units = 'cm', dpi = 600)
  }
}   

print(pldd)

# Combine plots using Patchwork #######
library(patchwork)
library(grid)

# FD plots 
ppfd1 <- wrap_plots(pfd,ncol = 1) 
ppfd2 <- wrap_plots (plfd, ncol=1) 
ppfd3 <- (ppfd1|ppfd2) & plot_annotation(title = "Frequency-dependent transmission")
print(ppfd3)
ggsave("patchwork_fd_all_test.png",ppfd3,width = 30, height = 40, units = 'cm', dpi = 600)

# dd plots 
ppdd1 <- wrap_plots(pdd,ncol=1) 
ppdd2 <- wrap_plots (pldd, ncol=1) 
ppdd3 <- (ppdd1|ppdd2) & plot_annotation(title = "Density-dependent transmission")
ggsave("patchwork_dd_all_test.png",ppdd3,width = 30, height = 40, units = 'cm', dpi = 600)

# Plot the average % of the population change in 100 years ##########
# Prepare the dataframe for boxplot using melt
mxfd<-list()
for (i in 1:length(mfd)){
  
  mxfd[[i]] <- melt(mfd[[i]], id.vars = c('time_y','time_d','run','model',"Ndiff"),
                  value.name = 'value', variable.name = 'class')
}
mxdd<-list()
for (i in 1:length(mdd)){
  
  mxdd[[i]] <- melt(mdd[[i]], id.vars = c('time_y','time_d','run','model',"Ndiff"),
                    value.name = 'value', variable.name = 'class')
}

####### summerize basic stat ####### 
# single run
for(i in 1:length(sfd)){
  print(sfd[[i]] |> group_by(class) |>
          dplyr::summarise(
            Max=max(value),
            Min=min(value)))
}
for(i in 1:length(sdd)){
  print(sdd[[i]] |> group_by(class) |>
          dplyr::summarise(
            Max=max(value),
            Min=min(value)))
}

# multiple runs
for (i in 1:length(mxfd)){
  print(mxfd[[i]] |> group_by(class) |>
          dplyr::summarise(Mean = mean(value),Max=max(value), Min=min(value)))
}

for (i in 1:length(mxdd)){
  print(mxdd[[i]] |> group_by(class) |>
          dplyr::summarise(Mean = mean(value),Max=max(value), Min=min(value)))
}
# combining multiple runs df list
mxfd2 <-data.table::rbindlist(mxfd)
mxdd2 <-data.table::rbindlist(mxdd)

#select some columns: Ndiff,run,model
dft_fd<-mxfd2 %>% 
  dplyr::select(Ndiff,run,model) %>% 
  drop_na() %>%
  distinct()

dft_fd$transmit<-c("FD")

dft_dd<-mxdd2 %>% 
  dplyr::select(Ndiff,run,model)%>% 
  drop_na()%>%
  distinct()

dft_dd$transmit<-c("DD")

#calculating mean,median by models
dft_fd |>
  group_by(model)|>
  summarise(Median = median(Ndiff),
            Mean = mean(Ndiff))

dft_dd |>
  group_by(model)|>
  summarise(Median = median(Ndiff),
            Mean = mean(Ndiff))

#recoding
dft_fd$model <- recode_factor(dft_fd$model, no_infection = "no infection" )
dft_dd$model <- recode_factor(dft_dd$model, no_infection = "no infection" )

# write .csv file
#write_csv(dft_fd,"df_ndiff_fd.csv")
#write_csv(dft_dd,"df_ndiff_dd.csv")

# load df 
# dft_dd <- read.csv("df_ndiff_dd.csv")
# dft_fd <- read.csv("df_ndiff_fd.csv")


# Boxplot #####
#1) FD
box_fd<-dft_fd %>%
  ggplot(aes(x= reorder(model,-Ndiff), y = Ndiff,fill=model)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.6, alpha=0.4) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)) +
  ggtitle("Gaur population change by FD infectious diseases in 100 years") +
  xlab("") +
  ylab("Population change (%)")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))
print(box_fd)
#ggsave("gaur_ndiff_box1_100sim_FD.png",box_fd,width = 22, height = 15, units = 'cm', dpi = 600)

#reorder by population change max-min
dft_fd$model2 <- factor(dft_fd$model, 
                        levels = c("no infection", 
                                   "FMD",
                                   "HS",
                                   "LSD",
                                   "Anthrax",
                                   "bTB",
                                   "Brucellosis"))
plt_fd<-dft_fd%>%ggbetweenstats(
  x=model2,
  y=Ndiff,
  k=0,
  plot.type = "boxviolin",
  pairwise.comparisons=F,
  bf.message = F,
  results.subtitle = FALSE,
  centrality.point.args = list(size = 2, color = "darkred"),
  title= "Gaur population change by FD infectious disease in 100 years",
  xlab = "",
  ylab = "Population change (%)",
  package = "ggsci",
  palette = "default_jco")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))

print(plt_fd)
#ggsave("gaur_ndiff_boxviolin_100sim_fd.png",plt_fd,width = 20, height = 15, units = 'cm', dpi = 600)

# 2) DD
box_dd<-dft_dd %>%
  ggplot(aes(x= reorder(model,-Ndiff), y = Ndiff,fill=model)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.6, alpha=0.4) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)) +
  ggtitle("Gaur population change by DD infectious disease models in 100 years") +
  xlab("") +
  ylab("Population change (%)")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))

print(box_dd)
#ggsave("gaur_ndiff_box_100sim_dd.png",box_dd,width = 22, height = 15, units = 'cm', dpi = 600)

#reorder by population change max-min
dft_dd$model2 <- factor(dft_dd$model, 
                        levels = c("no infection", 
                                   "Anthrax",
                                   "FMD",
                                   "Brucellosis",
                                   "HS",
                                   "LSD",
                                   "bTB"))

plt_dd<-dft_dd%>%ggbetweenstats(
  x=model2,
  y=Ndiff,
  k=0,
  plot.type = "boxviolin",
  pairwise.comparisons=F,
  bf.message = F,
  results.subtitle = FALSE,
  centrality.point.args = list(size = 2, color = "darkred"),
  title= "Gaur population change by DD infectious disease models in 100 years",
  xlab = "",
  ylab = "Population change (%)",
  package = "ggsci",
  palette = "default_jco")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))

print(plt_dd)
#ggsave("gaur_ndiff_boxviolin_100sim_dd.png",plt,width = 20, height = 15, units = 'cm', dpi = 600)

#compare boxplot of dd & fd disease transmission models
dft_all<-full_join(dft_dd,dft_fd)
View(dft_all)

#reorder by population change max-min
dft_all$model2 <- factor(dft_all$model, 
                     levels = c("no infection", 
                                "Anthrax",
                                "FMD",
                                "Brucellosis",
                                "HS",
                                "LSD",
                                "bTB"))

pltall<-grouped_ggbetweenstats(
  data =dft_all,
  x=model2,
  y=Ndiff,
  grouping.var = transmit,
  k=0,
  plot.type = "boxviolin",
  pairwise.comparisons=F,
  bf.message = F,
  results.subtitle = FALSE,
  centrality.label.args = list(size  = 3.5),
  centrality.point.args = list(size = 2, color = "darkred"),
  annotation.args  = list(title = "Gaur population change by infectious disease models in 100 years"), 
  plotgrid.args = list(nrow = 2),
  xlab = "",
  ylab = "Population change (%)",
  package = "ggsci",
  palette = "default_jco")+
  
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))

print(pltall)
#ggsave("gaur_ndiff_boxviolin_100sim_all.png",pltall,width = 22, height = 25, units = 'cm', dpi = 600)

# another boxplot style
dft_all$transmit[dft_all$model == "no infection" ] <- "no"
# Boxplot
boxall<-dft_all %>%
  ggplot(aes(x= model2, y = Ndiff, fill=transmit)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  # geom_jitter(color="black", size=0.6, alpha=0.4) +
  theme_ipsum() +
  theme(
    #legend.position="none",
    plot.title = element_text(size=11)) +
  ggtitle("Gaur population change by infectious disease models in 100 years") +
  xlab("") +
  ylab("Population change (%)")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))

print(boxall)
#ggsave("gaur_ndiff_box_100sim_all.png",boxall,width = 22, height = 15, units = 'cm', dpi = 600)

# DONE :) #
