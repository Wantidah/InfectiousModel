############## NOTES ##################
# Infectious disease modelling code for wild Bovidae populations 
# - 18 Nov 2022-

# # Code adapted from: https://github.com/dtsh2/ebola_model; 
# # Reference Article:https://doi.org/10.1098/rsif.2021.0638  (Transmission models indicate Ebola virus persistence in non-human primate populations is unlikely)

# Working on:  R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
# RStudio 2022.07.1+554 "Spotted Wakerobin" Release (7872775ebddc40635780ca1ed238934c3345c5de, 2022-07-22) for macOS
# Mozilla/5.0 (Macintosh; Intel Mac OS X 13_0_1)

# Note:this models are based on
# - Single, close population, no migration
# - Number of the population == 300
# - All parameters can be adapted if needed
# - All models have 3 age classes (c = calf, sa = Subadult, a = adult)
# unit == day (per day)

# Script Outline:
# 1) Set up the infectious disease model's function
    # - FD = Frequency-dependent model
    # - DD = Density-dependent model
# 2) Input parameters
# 3) Single run for all the models
# 4) Multiple run for all the models 
# 5) Plotting single and multiple runs (e.g. SIR plot, boxplot)

#rm(list=ls())
set.seed(111)

options(digits = 3, scipen = 999)

# set library path
# path <- '...'

#install.packages(c('EpiDynamics','dplyr','tidyverse','reshape2',
#                   'stringr','hrbrthemes','ggstatsplot','patchwork','grid',
#                    'factoextra','FactoMineR','readxl','ggrepel','PCAtools','data.table','wesanderson'), 
#                 lib = path)
#devtools::install_github('kevinblighe/PCAtools')
#if(!require(ggstatsplot)) install.packages('ggstatsplot', dep = TRUE, quiet = TRUE)

############## LOAD PACKAGES ##########

library(EpiDynamics)
library(dplyr)   
library(tidyverse)
library(reshape2) 
library(stringr)
library(ggplot2); theme_set(theme_bw())
library(hrbrthemes)
library(ggstatsplot)
library(patchwork)
library(grid)
library(factoextra)
library(FactoMineR)
library(readxl)
library(ggrepel)
library(PCAtools)
library(data.table)
library(wesanderson)

my_theme <- theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))

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
          num.min <- min(num, init[which(change[i, ] < 0)])
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
    results <- data.frame(time, a,  sa,  c)%>% 
      dplyr::mutate(N = rowSums(across(-c(time)), na.rm=TRUE))
    return(list(pars = pars, init = init2, time = time, results = results))
    
  }

############## MODEL 2 SI - Anthrax DD #####
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
        rate[3] <- Ic * rho_c
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
        rate[8] <- Isa * rho_sa
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
        rate[13] <- Ia *rho_a
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
      dplyr::mutate(N = rowSums(across(-c(time,S,I)), na.rm=TRUE))
    return(list(pars = pars, init = init2, time = time, results = results))
    
  }

############## MODEL 2 SI - Anthrax FD #####
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
        rate[3] <- Ic * rho_c
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
        rate[8] <- Isa * rho_sa
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
        rate[13] <-  Ia *rho_a
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
    results <- data.frame(time, 
                        Sc,  Ic,  Ssa, Isa, Sa, Ia)%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(N = rowSums(across(-c(time,S,I)), na.rm=TRUE))
    return(list(pars = pars, init = init2, time = time, results = results))
  }


############## MODEL 3 SEI - Bovine tuberculosis DD #####
model3dd =
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 29)
        change <- matrix(0, nrow = 29, ncol = 9)
        
        N <- Sc + Ec + Ic + Ssa + Esa + Isa + Sa + Ea + Ia 
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
        
        # subadult
        rate[12] <- beta_sa * Ssa * (Ic+Isa+Ia)
        change[12, ] <- c(0, 0, 0, -1, 1, 0, 0, 0, 0)
        rate[13] <- phi_sa * Esa 
        change[13, ] <- c(0, 0, 0, 0,-1, 1, 0, 0, 0)
        rate[14] <- rho_sa * gamma_sa * Isa
        change[14, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0)
        rate[15] <- delta_sa * Ssa
        change[15, ] <- c(0, 0, 0, -1, 0, 0, 1, 0, 0)  
        rate[16] <- delta_sa * Esa
        change[16, ] <- c(0, 0, 0, 0, -1, 0, 0, 1, 0) 
        rate[17] <-  delta_sa *  Isa
        change[17, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 1) 
        rate[18] <- mu_sa * Ssa
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
      dplyr::mutate(N = rowSums(across(c(S,E,I)), na.rm=TRUE))
    
    return(list(pars = pars, init = init2, time = time, results = results))
  }

############## MODEL 3 SEI - Bovine tuberculosis FD #####
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
        change[14, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0)
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
      dplyr::mutate(N = rowSums(across(c(S,E,I)), na.rm=TRUE)) 
    
    return(list(pars = pars, init = init2, time = time, results = results))
  }

############## MODEL 4 SIRS - Hemorrhagic septicemia DD #####
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
        rate[21] <-  mu_sa * Isa
        change[21, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 0)
        rate[22] <-  mu_sa * Rsa
        change[22, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0)
        rate[23] <-  epsilon * Ssa
        change[23, ] <- c(0, 0, 0, -1, 1, 0, 0, 0, 0)
        
        #adult
        rate[24] <- beta_a * Sa * (Ic+Isa+Ia)
        change[24, ] <- c(0, 0, 0, 0, 0, 0, -1, 1, 0)
        rate[25] <- (1-rho_a) * gamma_a * Ia
        change[25, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 1)
        rate[26] <-  rho_a * gamma_a * Ia
        change[26, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0)
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
          num.min <- min(num, init[which(change[i, ] < 0)])
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
    results <- data.frame(time, 
                        Sc,  Ic, Rc, Ssa, Isa, Rsa, Sa, Ia, Ra)%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))%>%
      dplyr::mutate(N = rowSums(across(c(S,I,R)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
  }

############## MODEL 4 SIRS - Hemorrhagic septicemia FD #####
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
        change[26, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0)
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
          num.min <- min(num, init[which(change[i, ] < 0)])
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
      dplyr::mutate(N = rowSums(across(c (S,I,R)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
  }


############## MODEL 5 SEIRS - Lumpy skin disease DD #####
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
          num.min <- min(num, init[which(change[i, ] < 0)])
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
      dplyr::mutate(N = rowSums(across(c(S,E,I,R)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
  }

############## MODEL 5 SEIRS - Lumpy skin disease FD #####
model5fd =
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 40)
        change <- matrix(0, nrow = 40, ncol = 12)
        
        N <- Sc + Ec + Ic + Rc + Ssa + Esa + Isa + Rsa + Sa + Ea + Ia + Ra
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
        rate[11] <- delta_c * Rc
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
        rate[34] <- (1- rho_a) * gamma_a * Ia
        change[34, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1)
        rate[35] <- rho_a * gamma_a * Ia
        change[35, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0)
        rate[36] <- omega_a *  Ra
        change[36, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1)
        rate[37] <-  mu_a * Sa
        change[37, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0) 
        rate[38] <- mu_a * Ea
        change[38, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0) 
        rate[39] <- mu_a * Ia
        change[39, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0) 
        rate[40] <- mu_a * Ra
        change[40, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1) 
        
        init <- c(Sc = Sc, Ec = Ec, Ic = Ic, Rc = Rc, 
                  Ssa = Ssa, Esa = Esa, Isa = Isa, Rsa = Rsa,
                  Sa = Sa, Ea = Ea, Ia = Ia, Ra = Ra)
        for (i in 1:40) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 0)])
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
    results <- data.frame(time, 
                        Sc, Ec, Ic, Rc ,Ssa, Esa, Isa, Rsa, Sa, Ea, Ia, Ra )%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))%>% 
      dplyr::mutate(N = rowSums(across(c(S,E,I,R)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
  }

############## MODEL 6 SEIRMS/E - Foot and mouth disease DD #####
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
        #Ia produce Ic
        rate[16] <- alpha * mu_bI * Ia
        change[16, ] <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        #Ra produce M (maternal derived immunity calf)
        rate[17] <- mu_b * Ra
        change[17, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
        rate[18] <- omega_m * M
        change[18, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1)
        rate[19] <- mu_c * M
        change[19, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0)
        #Sm go back to Ssa
        rate[20] <- (delta_c + omega_m) * Sm  
        change[20, ] <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        #Sm go back to Ec
        rate[21] <-  beta_c * Sm * (Ic+Isa+Ia)
        change[21, ] <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[22] <- mu_c * Sm
        change[22, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[23] <- epsilon * Sc  
        change[23, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        
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
        change[38, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0) 
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
          num.min <- min(num, init[which(change[i, ] < 0)])
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
      dplyr::mutate(N = rowSums(across(c(S,E,I,R,M)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
    
  }
############## MODEL 6 SEIRMS/E - Foot and mouth disease FD #####
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
        #Ia produce Ic
        rate[16] <- alpha * mu_bI * Ia
        change[16, ] <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        #Ra produce M (maternal derived immunity calf)
        rate[17] <- mu_b * Ra
        change[17, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
        rate[18] <- omega_m * M
        change[18, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1)
        rate[19] <- mu_c * M
        change[19, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0)
        #Sm go back to Ssa
        rate[20] <- (delta_c + omega_m) * Sm  
        change[20, ] <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        #Sm go back to Ec
        rate[21] <-  beta_c * Sm * (Ic+Isa+Ia)/N
        change[21, ] <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[22] <- mu_c * Sm
        change[22, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[23] <- epsilon * Sc  
        change[23, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        
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
          num.min <- min(num, init[which(change[i, ] < 0)])
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
      dplyr::mutate(N = rowSums(across(c(S,E,I,R,M)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
    
  }

############## MODEL 7 SEIRMS/E - Bovine brucellosis DD #####
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
        #Ia produce Ic
        rate[16] <- alpha * mu_bI * Ia
        change[16, ] <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        #Ra produce M (maternal derived immunity calf)
        rate[17] <- mu_b * Ra
        change[17, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
        rate[18] <- omega_m * M
        change[18, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1)
        rate[19] <- mu_c * M
        change[19, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0)
        #Sm go back to Ssa
        rate[20] <- (delta_c + omega_m) * Sm  
        change[20, ] <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        #Sm go back to Ec
        rate[21] <-  beta_c * Sm * (Ic+Isa+Ia)
        change[21, ] <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[22] <- mu_c * Sm
        change[22, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[23] <- epsilon * Sc  
        change[23, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        
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
          num.min <- min(num, init[which(change[i, ] < 0)])
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
      dplyr::mutate(N = rowSums(across(c(S,E,I,R,M)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
    
  }

############## MODEL 7 SEIRMS/E - Bovine brucellosis FD #####
model7fd=
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
        rate[5] <- (1-rho_c) * gamma_c * Ic
        change[5, ] <- c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[6] <- rho_c * gamma_c * Ic
        change[6, ] <- c(0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[7] <- omega_c *  Rc
        change[7, ] <- c(1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[8] <-  delta_c * Sc
        change[8, ] <- c(-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)  
        rate[9] <- delta_c * Ec
        change[9, ] <- c(0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0) 
        rate[10] <- delta_c * Ic
        change[10, ] <- c(0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0) 
        rate[11] <- delta_c * Rc
        change[11, ] <- c(0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0) 
        rate[12] <- mu_c * Sc
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
        rate[20] <- (delta_c + omega_m) * Sm  
        change[20, ] <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        #Sm go back to Ec
        rate[21] <-  beta_c * Sm * (Ic+Isa+Ia)/N
        change[21, ] <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[22] <- mu_c * Sm
        change[22, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
        rate[23] <- epsilon * Sc  
        change[23, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        
        # saubadult
        rate[24] <- beta_sa * Ssa * (Ic+Isa+Ia)/N
        change[24, ] <- c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
        rate[25] <- phi_sa * Esa 
        change[25, ] <- c(0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0)
        rate[26] <- (1-rho_sa) * gamma_sa * Isa
        change[26, ] <- c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0)
        rate[27] <-  rho_sa * gamma_sa * Isa
        change[27, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0)
        rate[28] <- omega_sa * Rsa
        change[28, ] <- c(0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0)
        rate[29] <- delta_sa * Ssa
        change[29, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0)  
        rate[30] <- delta_sa * Esa
        change[30, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0) 
        rate[31] <-  delta_sa *  Isa
        change[31, ] <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0) 
        rate[32] <-  delta_sa *  Rsa
        change[32, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0) 
        rate[33] <- mu_sa * Ssa
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
        rate[41] <- (1- rho_a) * gamma_a * Ia
        change[41, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0)
        rate[42] <- rho_a * gamma_a * Ia
        change[42, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0)
        rate[43] <- omega_a *  Ra
        change[43, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0)
        rate[44] <-  mu_a * Sa
        change[44, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0) 
        rate[45] <- mu_a * Ea
        change[45, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0) 
        rate[46] <- mu_a * Ia
        change[46, ] <-  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0) 
        rate[47] <- mu_a * Ra
        change[47, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0) 
        
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
    results <- data.frame(time, 
                        Sc, Ec, Ic, Rc, Ssa, Esa, Isa, Rsa, Sa, Ea, Ia, Ra,  M, Sm )%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc,Sm)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))%>%
      dplyr::mutate(N = rowSums(across(c(S,E,I,R,M)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
    
  }

#############  2) Set up parameters before running ##################  
# setting end time and number of replications (n_rep)  
# in the main text we used end_time = 100 * 365, n_rep = 100
end.time <- 100 * 365
n_rep <- 100

# population parameters #############
# gaur population 
N = 300 

#estimate the age structure proportion
#calf (c) :subadult (sa) :adult (a) ratio
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

# Disease parameters ####
# MODEL 1 - NO infection #####
initials_m1 <- c(c = c, sa = sa, a = a )
parameters_m1 <- c(mu_b  = 0.34/365, 
                   mu_c  = 0.27/365, 
                   mu_sa = 0.15/365,
                   mu_a  = 0.165/365,
                   delta_c  = 1/365,
                   delta_sa = 1/(3*365),
                   N = sum(initials_m1), 
                   tau = 1)
t(t(parameters_m1))

# MODEL 2 SI - Anthrax ####
initials_m2 <- c(Sc = c, Ic = 0, Ssa = sa, Isa = 0, Sa = (a-1), Ia = 1 )

#rho = 1;  fatality 100%
## m2fd : beta 0.01 ####
parameters_m2fd <- c(
  beta_c   = 0.01, 
  beta_sa  = 0.01, 
  beta_a   = 0.01,
  gamma_c  = 1, 
  gamma_sa = 1, 
  gamma_a  = 1, 
  rho_c    = 1, 
  rho_sa   = 1,  
  rho_a    = 1,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365, 
  mu_a     = 0.165/365,
  delta_c  = 1/365, 
  delta_sa = 1/(3*365),
  N = sum(initials_m2),
  tau = 1)
t(t(parameters_m2fd))

## m2dd* : beta 0.01/N ####
parameters_m2dd <- c(
  beta_c   = 0.01/N, 
  beta_sa  = 0.01/N, 
  beta_a   = 0.01/N,
  gamma_c  = 1, 
  gamma_sa = 1, 
  gamma_a  = 1, 
  rho_c    = 1, 
  rho_sa   = 1,  
  rho_a    = 1,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365, 
  mu_a     = 0.165/365,
  delta_c  = 1/365, 
  delta_sa = 1/(3*365),
  N = sum(initials_m2),
  tau = 1)
t(t(parameters_m2dd))

# MODEL 3 SEI - Bovine tuberculosis #####
initials_m3 <- c(Sc = c, Ec = 0, Ic = 0, Ssa = sa, Esa = 0, Isa = 0, Sa = (a-1), Ea = 0, Ia = 1)

# Lifelong infection
# Ia birth rate reduce by = 27%

## m3dd : beta 0.00143 ####
parameters_m3dd <- c( 
  beta_c   = 0.00143, 
  beta_sa  = 0.00143, 
  beta_a   = 0.00143,
  phi_c    = 0.21/30, 
  phi_sa   = 0.21/30, 
  phi_a    = 0.21/30,
  gamma_c  = 0, 
  gamma_sa = 0, 
  gamma_a  = 0,
  rho_c    = 0, 
  rho_sa   = 0, 
  rho_a    = 0.1, 
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(1-0.27),    
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m3), 
  tau=1)
t(t(parameters_m3dd))
 
## m3fd* : beta 0.00143*N  ####
parameters_m3fd <- c( 
  beta_c   = 0.00143*N, 
  beta_sa  = 0.00143*N, 
  beta_a   = 0.00143*N,
  phi_c    = 0.21/30, 
  phi_sa   = 0.21/30, 
  phi_a    = 0.21/30,
  gamma_c  = 0, 
  gamma_sa = 0, 
  gamma_a  = 0,
  rho_c    = 0, 
  rho_sa   = 0, 
  rho_a    = 0.1, 
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(1-0.27),    
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m3), 
  tau=1)
t(t(parameters_m3fd))

## m3dd2 : beta 2.7*10^-5 ####
parameters_m3dd2 <- c( 
  beta_c   = (2.7*10^-5), 
  beta_sa  = (2.7*10^-5), 
  beta_a   = (2.7*10^-5),
  phi_c    = 0.21/30, 
  phi_sa   = 0.21/30, 
  phi_a    = 0.21/30,
  gamma_c  = 0, 
  gamma_sa = 0, 
  gamma_a  = 0,
  rho_c    = 0, 
  rho_sa   = 0, 
  rho_a    = 0.1, 
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(1-0.27),    
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m3), 
  tau=1)
t(t(parameters_m3dd2))

## m3fd2* : beta 2.7*10^-5*N ####
parameters_m3fd2 <- c( 
  beta_c   = (2.7*10^-5)*N, 
  beta_sa  = (2.7*10^-5)*N,
  beta_a   = (2.7*10^-5)*N,
  phi_c    = 0.21/30, 
  phi_sa   = 0.21/30, 
  phi_a    = 0.21/30,
  gamma_c  = 0, 
  gamma_sa = 0, 
  gamma_a  = 0,
  rho_c    = 0, 
  rho_sa   = 0, 
  rho_a    = 0.1, 
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(1-0.27),    
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m3), 
  tau=1)
t(t(parameters_m3fd2))

## m3fd3 : beta 0.00143 ####
parameters_m3fd3 <- c( 
  beta_c   = 0.00143, 
  beta_sa  = 0.00143, 
  beta_a   = 0.00143,
  phi_c    = 0.21/30, 
  phi_sa   = 0.21/30, 
  phi_a    = 0.21/30,
  gamma_c  = 0, 
  gamma_sa = 0, 
  gamma_a  = 0,
  rho_c    = 0, 
  rho_sa   = 0, 
  rho_a    = 0.1, 
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(1-0.27),    
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m3), 
  tau=1)
t(t(parameters_m3fd3))

## m3dd3* : beta 0.00143/N ####
parameters_m3dd3 <- c( 
  beta_c   = 0.00143/N, 
  beta_sa  = 0.00143/N, 
  beta_a   = 0.00143/N,
  phi_c    = 0.21/30, 
  phi_sa   = 0.21/30, 
  phi_a    = 0.21/30,
  gamma_c  = 0, 
  gamma_sa = 0, 
  gamma_a  = 0,
  rho_c    = 0, 
  rho_sa   = 0, 
  rho_a    = 0.1, 
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(1-0.27),    
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m3), 
  tau=1)
t(t(parameters_m3dd3))

## m3fd4 : beta 0.0063 ####
parameters_m3fd4 <- c( 
  beta_c   = 0.0063, 
  beta_sa  = 0.0063, 
  beta_a   = 0.0063,
  phi_c    = 0.21/30, 
  phi_sa   = 0.21/30, 
  phi_a    = 0.21/30,
  gamma_c  = 0, 
  gamma_sa = 0, 
  gamma_a  = 0,
  rho_c    = 0, 
  rho_sa   = 0, 
  rho_a    = 0.1, 
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(1-0.27),    
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m3), 
  tau=1)
t(t(parameters_m3fd4))

## m3dd4* : beta 0.0063/N ####
parameters_m3dd4 <- c( 
  beta_c   = 0.0063/N, 
  beta_sa  = 0.0063/N, 
  beta_a   = 0.0063/N,
  phi_c    = 0.21/30, 
  phi_sa   = 0.21/30, 
  phi_a    = 0.21/30,
  gamma_c  = 0, 
  gamma_sa = 0, 
  gamma_a  = 0,
  rho_c    = 0, 
  rho_sa   = 0, 
  rho_a    = 0.1, 
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(1-0.27),    
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m3), 
  tau=1)
t(t(parameters_m3dd4))

# MODEL 4 SIRS - Hemorrhagic septicemia #####
initials_m4 <- c(Sc = c, Ic = 0, Rc = 0, Ssa = sa, Isa = 0, Rsa = 0, Sa = (a-1), Ia = 1, Ra = 0)

## m4dd : beta 0.335, rho 5.26E-03 ####
parameters_m4dd <- c( 
  beta_c   = 0.335, 
  beta_sa  = 0.335, 
  beta_a   = 0.335,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 5.26E-03,
  rho_sa   = 5.26E-03,
  rho_a    = 5.26E-03,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4dd))

## m4fd* : beta 0.335*N, rho 5.26E-03 ####
parameters_m4fd <- c( 
  beta_c   = 0.335*N, 
  beta_sa  = 0.335*N, 
  beta_a   = 0.335*N,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 5.26E-03,
  rho_sa   = 5.26E-03,
  rho_a    = 5.26E-03,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4fd))

## m4dd2 : beta 0.335, rho 0.058 ####
parameters_m4dd2 <- c( 
  beta_c   = 0.335,
  beta_sa  = 0.335,
  beta_a   = 0.335,
  gamma_c  = 1/3,
  gamma_sa = 1/3,
  gamma_a  = 1/3,
  rho_c    = 0.058,
  rho_sa   = 0.058,
  rho_a    = 0.058,
  omega_c  = 1/180,
  omega_sa = 1/180,
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4dd2))

## m4fd2* : beta 0.335*N, rho 0.058 ####
parameters_m4fd2 <- c( 
  beta_c   = 0.335*N,
  beta_sa  = 0.335*N,
  beta_a   = 0.335*N,
  gamma_c  = 1/3,
  gamma_sa = 1/3,
  gamma_a  = 1/3,
  rho_c    = 0.058,
  rho_sa   = 0.058,
  rho_a    = 0.058,
  omega_c  = 1/180,
  omega_sa = 1/180,
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4fd2))

## m4fd3 : beta 0.335, rho 5.26E-03 ####
parameters_m4fd3 <- c( 
  beta_c   = 0.335,
  beta_sa  = 0.335,
  beta_a   = 0.335,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 5.26E-03,
  rho_sa   = 5.26E-03,
  rho_a    = 5.26E-03,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4fd3))

## m4dd3* : beta 0.335/N, rho 5.26E-03 ####
parameters_m4dd3 <- c(
  beta_c   = 0.335/N,
  beta_sa  = 0.335/N,
  beta_a   = 0.335/N,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 5.26E-03,
  rho_sa   = 5.26E-03,
  rho_a    = 5.26E-03,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4dd3))

## m4fd4 : beta 0.335, rho 0.058 ####
parameters_m4fd4 <- c( 
  beta_c   = 0.335,
  beta_sa  = 0.335,
  beta_a   = 0.335,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 0.058,
  rho_sa   = 0.058,
  rho_a    = 0.058,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4fd4))

## m4dd4* : beta 0.335/N, rho 0.058 ####
parameters_m4dd4 <- c(
  beta_c   = 0.335/N, 
  beta_sa  = 0.335/N,
  beta_a   = 0.335/N,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 0.058,
  rho_sa   = 0.058,
  rho_a    = 0.058,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4dd4))

## m4dd5 : beta 0.55, rho 5.26E-03 ####
parameters_m4dd5 <- c( 
  beta_c   = 0.55,
  beta_sa  = 0.55,
  beta_a   = 0.55,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 5.26E-03,
  rho_sa   = 5.26E-03,
  rho_a    = 5.26E-03,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4dd5))

## m4fd5* : beta 0.55*N, rho 5.26E-03 ####
parameters_m4fd5 <- c( 
  beta_c   = 0.55*N,
  beta_sa  = 0.55*N,
  beta_a   = 0.55*N,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 5.26E-03,
  rho_sa   = 5.26E-03,
  rho_a    = 5.26E-03,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4fd5))

## m4dd6 : beta 0.55, rho 0.058 ####
parameters_m4dd6 <- c( 
  beta_c   = 0.55,
  beta_sa  = 0.55,
  beta_a   = 0.55,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 0.058,
  rho_sa   = 0.058,
  rho_a    = 0.058,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4dd6))

## m4fd6* : beta 0.55*N, rho 0.058 ####
parameters_m4fd6 <- c( 
  beta_c   = 0.55*N,
  beta_sa  = 0.55*N,
  beta_a   = 0.55*N,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 0.058,
  rho_sa   = 0.058,
  rho_a    = 0.058,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4fd6))

## m4fd7: beta 0.55, rho 5.26E-03 ####
parameters_m4fd7 <- c( 
  beta_c   = 0.55,
  beta_sa  = 0.55,
  beta_a   = 0.55,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 5.26E-03,
  rho_sa   = 5.26E-03,
  rho_a    = 5.26E-03,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4fd7))

## m4dd7*: beta 0.55/N, rho 5.26E-03 ####
parameters_m4dd7 <- c(
  beta_c   = 0.55/N,
  beta_sa  = 0.55/N,
  beta_a   = 0.55/N,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 5.26E-03,
  rho_sa   = 5.26E-03,
  rho_a    = 5.26E-03,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4dd7))

## m4fd8: beta 0.55, rho 0.058 ####
parameters_m4fd8 <- c( 
  beta_c   = 0.55,
  beta_sa  = 0.55,
  beta_a   = 0.55,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 0.058,
  rho_sa   = 0.058,
  rho_a    = 0.058,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4fd8))

## m4dd8*: beta 0.55/N, rho 0.058 ####
parameters_m4dd8 <- c(
  beta_c   = 0.55/N,
  beta_sa  = 0.55/N,
  beta_a   = 0.55/N,
  gamma_c  = 1/3, 
  gamma_sa = 1/3, 
  gamma_a  = 1/3,
  rho_c    = 0.058,
  rho_sa   = 0.058,
  rho_a    = 0.058,
  omega_c  = 1/180,
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365,
  mu_a     = 0.165/365,
  delta_c  = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials_m4),
  tau=1)
t(t(parameters_m4dd8))

# MODEL 5 SEIRS - Lumpy skin disease #####
initials_m5 <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = (a-1), Ea = 0, Ia = 1, Ra = 0)

## m5dd: beta 0.008 ####
parameters_m5dd <- c( 
  beta_c   = 0.008, 
  beta_sa  = 0.008, 
  beta_a   = 0.008, 
  phi_c    = 1/7,
  phi_sa   = 1/7,
  phi_a    = 1/7,
  gamma_c  = 1/35, 
  gamma_sa = 1/35,
  gamma_a  = 1/35,
  rho_c    = 0.05,
  rho_sa   = 0.03, 
  rho_a    = 0.01, 
  omega_c  = 1/180, 
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(0.9), 
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365, 
  mu_a     = 0.165/365,
  delta_c  = 1/365, 
  delta_sa = 1/(3*365),
  N = sum(initials_m5),
  tau=1)
t(t(parameters_m5dd))

## m5fd*: beta 0.008*N  ####
parameters_m5fd <- c( 
  beta_c   = 0.008*N, 
  beta_sa  = 0.008*N, 
  beta_a   = 0.008*N, 
  phi_c    = 1/7,
  phi_sa   = 1/7,
  phi_a    = 1/7,
  gamma_c  = 1/35, 
  gamma_sa = 1/35,
  gamma_a  = 1/35,
  rho_c    = 0.05,
  rho_sa   = 0.03, 
  rho_a    = 0.01, 
  omega_c  = 1/180, 
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(0.9), 
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365, 
  mu_a     = 0.165/365,
  delta_c  = 1/365, 
  delta_sa = 1/(3*365),
  N = sum(initials_m5),
  tau=1)
t(t(parameters_m5fd))

## m5fd2: beta 0.032  ####
parameters_m5fd2 <- c( 
  beta_c   = 0.032, 
  beta_sa  = 0.032, 
  beta_a   = 0.032, 
  phi_c    = 1/7,
  phi_sa   = 1/7,
  phi_a    = 1/7,
  gamma_c  = 1/35, 
  gamma_sa = 1/35,
  gamma_a  = 1/35,
  rho_c    = 0.05,
  rho_sa   = 0.03, 
  rho_a    = 0.01, 
  omega_c  = 1/180, 
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(0.9), 
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365, 
  mu_a     = 0.165/365,
  delta_c  = 1/365, 
  delta_sa = 1/(3*365),
  N = sum(initials_m5),
  tau=1)
t(t(parameters_m5fd2))

## m5dd2: beta 0.032/N  ####
parameters_m5dd2 <- c( 
  beta_c   = 0.032/N, 
  beta_sa  = 0.032/N, 
  beta_a   = 0.032/N, 
  phi_c    = 1/7,
  phi_sa   = 1/7,
  phi_a    = 1/7,
  gamma_c  = 1/35, 
  gamma_sa = 1/35,
  gamma_a  = 1/35,
  rho_c    = 0.05,
  rho_sa   = 0.03, 
  rho_a    = 0.01, 
  omega_c  = 1/180, 
  omega_sa = 1/180, 
  omega_a  = 1/180,
  epsilon  = 2e-5,
  mu_b     = 0.34/365,  
  mu_bI    = (0.34/365)*(0.9), 
  mu_c     = 0.27/365,  
  mu_sa    = 0.15/365, 
  mu_a     = 0.165/365,
  delta_c  = 1/365, 
  delta_sa = 1/(3*365),
  N = sum(initials_m5),
  tau=1)
t(t(parameters_m5dd2))

# MODEL 6 SEIRMS/E - Foot and mouth disease  #####
initials_m6 <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, 
                 Sa = (a-1), Ea = 0, Ia = 1, Ra = 0)

## m6dd : beta 0.026 ####
parameters_m6dd <- c( 
  beta_c   = 0.026, 
  beta_sa  = 0.026, 
  beta_a   = 0.026,
  phi_c    = 1/8, 
  phi_sa   = 1/6, 
  phi_a    = 1/6,
  gamma_c  = 1/5, 
  gamma_sa = 1/5, 
  gamma_a  = 1/5,
  rho_c    = 0.1, 
  rho_sa   = 0.05, 
  rho_a    = 0.03, 
  alpha    = 0.5,
  omega_c  = (1/120), 
  omega_sa =  (1/120), 
  omega_a  = (1/565), 
  omega_m  = (1/144),
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_bI    = (0.34/365)*(0.9),
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365, 
  mu_a     = 0.165/365,
  delta_c  = 1/365, 
  delta_sa = 1/(3*365),
  N = sum(initials_m6),
  tau=1)
t(t(parameters_m6dd))

## m6fd* : beta 0.026*N ####
parameters_m6fd<- c( 
    beta_c   = 0.026*N, 
    beta_sa  = 0.026*N, 
    beta_a   = 0.026*N,
    phi_c    = 1/8, 
    phi_sa   = 1/6, 
    phi_a    = 1/6,
    gamma_c  = 1/5, 
    gamma_sa = 1/5, 
    gamma_a  = 1/5,
    rho_c    = 0.1, 
    rho_sa   = 0.05, 
    rho_a    = 0.03, 
    alpha    = 0.5,
    omega_c  = (1/120), 
    omega_sa =  (1/120), 
    omega_a  = (1/565), 
    omega_m  = (1/144),
    epsilon  = 2e-5,
    mu_b     = 0.34/365, 
    mu_bI    = (0.34/365)*(0.9),
    mu_c     = 0.27/365, 
    mu_sa    = 0.15/365, 
    mu_a     = 0.165/365,
    delta_c  = 1/365, 
    delta_sa = 1/(3*365),
    N = sum(initials_m6),
    tau=1)
t(t(parameters_m6fd))

## m6dd2 : beta 21 ####
parameters_m6dd2 <- c( 
  beta_c   = 21.84, 
  beta_sa  = 21.84, 
  beta_a   = 21.84,
  phi_c    = 1/8, 
  phi_sa   = 1/6, 
  phi_a    = 1/6,
  gamma_c  = 1/5, 
  gamma_sa = 1/5, 
  gamma_a  = 1/5,
  rho_c    = 0.1, 
  rho_sa   = 0.05, 
  rho_a    = 0.03, 
  alpha    = 0.5,
  omega_c  = (1/120), 
  omega_sa =  (1/120), 
  omega_a  = (1/565), 
  omega_m  = (1/144),
  epsilon  = 2e-5,
  mu_b     = 0.34/365, 
  mu_bI    = (0.34/365)*(0.9),
  mu_c     = 0.27/365, 
  mu_sa    = 0.15/365, 
  mu_a     = 0.165/365,
  delta_c  = 1/365, 
  delta_sa = 1/(3*365),
  N = sum(initials_m6),
  tau=1)
t(t(parameters_m6dd2))

## m6fd2* : beta 21*N ####
parameters_m6fd2 <- c( 
    beta_c   = 21.84*N, 
    beta_sa  = 21.84*N, 
    beta_a   = 21.84*N,
    phi_c    = 1/8, 
    phi_sa   = 1/6, 
    phi_a    = 1/6,
    gamma_c  = 1/5, 
    gamma_sa = 1/5, 
    gamma_a  = 1/5,
    rho_c    = 0.1, 
    rho_sa   = 0.05, 
    rho_a    = 0.03, 
    alpha    = 0.5,
    omega_c  = (1/120), 
    omega_sa =  (1/120), 
    omega_a  = (1/565), 
    omega_m  = (1/144),
    epsilon  = 2e-5,
    mu_b     = 0.34/365, 
    mu_bI    = (0.34/365)*(0.9),
    mu_c     = 0.27/365, 
    mu_sa    = 0.15/365, 
    mu_a     = 0.165/365,
    delta_c  = 1/365, 
    delta_sa = 1/(3*365),
    N = sum(initials_m6),
    tau=1)
t(t(parameters_m6fd2))

## m6fd3 : beta 0.115 ####
parameters_m6fd3 <- c( 
    beta_c   = 0.115, 
    beta_sa  = 0.115, 
    beta_a   = 0.115,
    phi_c    = 1/8, 
    phi_sa   = 1/6, 
    phi_a    = 1/6,
    gamma_c  = 1/5, 
    gamma_sa = 1/5, 
    gamma_a  = 1/5,
    rho_c    = 0.1, 
    rho_sa   = 0.05, 
    rho_a    = 0.03, 
    alpha    = 0.5,
    omega_c  = (1/120), 
    omega_sa =  (1/120), 
    omega_a  = (1/565), 
    omega_m  = (1/144),
    epsilon  = 2e-5,
    mu_b     = 0.34/365, 
    mu_bI    = (0.34/365)*(0.9),
    mu_c     = 0.27/365, 
    mu_sa    = 0.15/365, 
    mu_a     = 0.165/365,
    delta_c  = 1/365, 
    delta_sa = 1/(3*365),
    N = sum(initials_m6),
    tau=1)
t(t(parameters_m6fd3))

## m6dd3* : beta 0.115/N ####
parameters_m6dd3 <- c( 
    beta_c   = 0.115/N, 
    beta_sa  = 0.115/N, 
    beta_a   = 0.115/N,
    phi_c    = 1/8, 
    phi_sa   = 1/6, 
    phi_a    = 1/6,
    gamma_c  = 1/5, 
    gamma_sa = 1/5, 
    gamma_a  = 1/5,
    rho_c    = 0.1, 
    rho_sa   = 0.05, 
    rho_a    = 0.03, 
    alpha    = 0.5,
    omega_c  = (1/120), 
    omega_sa =  (1/120), 
    omega_a  = (1/565), 
    omega_m  = (1/144),
    epsilon  = 2e-5,
    mu_b     = 0.34/365, 
    mu_bI    = (0.34/365)*(0.9),
    mu_c     = 0.27/365, 
    mu_sa    = 0.15/365, 
    mu_a     = 0.165/365,
    delta_c  = 1/365, 
    delta_sa = 1/(3*365),
    N = sum(initials_m6),
    tau=1)
t(t(parameters_m6dd3))

## m6fd4 : beta 21.84 ####
parameters_m6fd4 <- c( 
    beta_c   = 21.84, 
    beta_sa  = 21.84, 
    beta_a   = 21.84,
    phi_c    = 1/8, 
    phi_sa   = 1/6, 
    phi_a    = 1/6,
    gamma_c  = 1/5, 
    gamma_sa = 1/5, 
    gamma_a  = 1/5,
    rho_c    = 0.1, 
    rho_sa   = 0.05, 
    rho_a    = 0.03, 
    alpha    = 0.5,
    omega_c  = (1/120), 
    omega_sa =  (1/120), 
    omega_a  = (1/565), 
    omega_m  = (1/144),
    epsilon  = 2e-5,
    mu_b     = 0.34/365, 
    mu_bI    = (0.34/365)*(0.9),
    mu_c     = 0.27/365, 
    mu_sa    = 0.15/365, 
    mu_a     = 0.165/365,
    delta_c  = 1/365, 
    delta_sa = 1/(3*365),
    N = sum(initials_m6),
    tau=1)
t(t(parameters_m6fd4))

## m6dd4 : beta 21.84/N ####
parameters_m6dd4 <- c( 
    beta_c   = 21.84/N, 
    beta_sa  = 21.84/N, 
    beta_a   = 21.84/N,
    phi_c    = 1/8, 
    phi_sa   = 1/6, 
    phi_a    = 1/6,
    gamma_c  = 1/5, 
    gamma_sa = 1/5, 
    gamma_a  = 1/5,
    rho_c    = 0.1, 
    rho_sa   = 0.05, 
    rho_a    = 0.03, 
    alpha    = 0.5,
    omega_c  = (1/120), 
    omega_sa = (1/120), 
    omega_a  = (1/565), 
    omega_m  = (1/144),
    epsilon  = 2e-5,
    mu_b     = 0.34/365, 
    mu_bI    = (0.34/365)*(0.9),
    mu_c     = 0.27/365, 
    mu_sa    = 0.15/365, 
    mu_a     = 0.165/365,
    delta_c  = 1/365, 
    delta_sa = 1/(3*365),
    N = sum(initials_m6),
    tau=1)
t(t(parameters_m6dd4))


# MODEL 7 SEIRMS/E - Brucellosis  #####
initials_m7 <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, 
                 Sa = (a-1), Ea = 0, Ia = 1, Ra = 0)

## m7dd : beta 0.002/365 ####
parameters_m7dd <- c( 
    beta_c   = 0.002/365, 
    beta_sa  = 0.002/365, 
    beta_a   = 0.002/365, 
    phi_c    = 1/14, 
    phi_sa   = 1/14, 
    phi_a    = 1/14,
    gamma_c  = 1/(2*365), 
    gamma_sa = 1/(2*365), 
    gamma_a  = 1/(2*365),
    rho_c    = 0.1, 
    rho_sa   = 0.05, 
    rho_a    = 0.03, 
    alpha    = 0.9,   
    omega_c  = 1/180, 
    omega_sa =  1/180, 
    omega_a  = 1/180, 
    omega_m  = 1/180,
    epsilon  = 2e-5,
    mu_b     = 0.34/365, 
    mu_bI    = (0.34/365)*(0.5), 
    mu_c     = 0.27/365, 
    mu_sa    = 0.15/365,
    mu_a     = 0.165/365,
    delta_c  = 1/365,
    delta_sa = 1/(3*365),
    N = sum(initials_m7),
    tau=1)
t(t(parameters_m7dd))

## m7fd* : beta (0.002/365)*N ####
parameters_m7fd <- c( 
    beta_c   = (0.002/365)*N, 
    beta_sa  = (0.002/365)*N, 
    beta_a   = (0.002/365)*N,
    phi_c    = 1/14, 
    phi_sa   = 1/14, 
    phi_a    = 1/14,
    gamma_c  = 1/(2*365), 
    gamma_sa = 1/(2*365), 
    gamma_a  = 1/(2*365),
    rho_c    = 0.1, 
    rho_sa   = 0.05, 
    rho_a    = 0.03, 
    alpha    = 0.9,   
    omega_c  = 1/180, 
    omega_sa =  1/180, 
    omega_a  = 1/180, 
    omega_m  = 1/180,
    epsilon  = 2e-5,
    mu_b     = 0.34/365, 
    mu_bI    = (0.34/365)*(0.5), 
    mu_c     = 0.27/365, 
    mu_sa    = 0.15/365,
    mu_a     = 0.165/365,
    delta_c  = 1/365,
    delta_sa = 1/(3*365),
    N = sum(initials_m7),
    tau=1)
t(t(parameters_m7fd))

## m7fd2 : beta 2/365 ####
parameters_m7fd2 <- c( 
    beta_c   = 2/365, 
    beta_sa  = 2/365, 
    beta_a   = 2/365,
    phi_c    = 1/14, 
    phi_sa   = 1/14, 
    phi_a    = 1/14,
    gamma_c  = 1/(2*365), 
    gamma_sa = 1/(2*365), 
    gamma_a  = 1/(2*365),
    rho_c    = 0.1, 
    rho_sa   = 0.05, 
    rho_a    = 0.03, 
    alpha    = 0.9,   
    omega_c  = 1/180, 
    omega_sa = 1/180, 
    omega_a  = 1/180, 
    omega_m  = 1/180,
    epsilon  = 2e-5,
    mu_b     = 0.34/365, 
    mu_bI    = (0.34/365)*(0.5), 
    mu_c     = 0.27/365, 
    mu_sa    = 0.15/365,
    mu_a     = 0.165/365,
    delta_c  = 1/365,
    delta_sa = 1/(3*365),
    N = sum(initials_m7),
    tau=1)
t(t(parameters_m7fd2))

## m7dd2 : beta (2/365)/N ####
parameters_m7dd2 <- c( 
    beta_c   = (2/365)/N, 
    beta_sa  = (2/365)/N, 
    beta_a   = (2/365)/N,
    phi_c    = 1/14, 
    phi_sa   = 1/14, 
    phi_a    = 1/14,
    gamma_c  = 1/(2*365), 
    gamma_sa = 1/(2*365), 
    gamma_a  = 1/(2*365),
    rho_c    = 0.1, 
    rho_sa   = 0.05, 
    rho_a    = 0.03, 
    alpha    = 0.9,   
    omega_c  = 1/180, 
    omega_sa = 1/180, 
    omega_a  = 1/180, 
    omega_m  = 1/180,
    epsilon  = 2e-5,
    mu_b     = 0.34/365, 
    mu_bI    = (0.34/365)*(0.5), 
    mu_c     = 0.27/365, 
    mu_sa    = 0.15/365,
    mu_a     = 0.165/365,
    delta_c  = 1/365,
    delta_sa = 1/(3*365),
    N = sum(initials_m7),
    tau=1)
t(t(parameters_m7dd2))

# 3) SINGLE RUNS - ALL MODELS  ###############################
#> m1 : No infection 
res_model1 <- model1(pars = parameters_m1, init = initials_m1,
                end.time = end.time)

#> m2: Anthrax 
res_model2dd<-model2dd(pars = parameters_m2dd, init = initials_m2,
                       end.time = end.time)
res_model2fd<-model2fd(pars = parameters_m2fd, init = initials_m2,
                       end.time = end.time)

#> m3: bTB
res_model3dd<-model3dd(pars = parameters_m3dd, init = initials_m3,
                        end.time = end.time)
res_model3dd2<-model3dd(pars = parameters_m3dd2, init = initials_m3,
                        end.time = end.time)
res_model3dd3<-model3dd(pars = parameters_m3dd3, init = initials_m3,
                        end.time = end.time)
res_model3dd4<-model3dd(pars = parameters_m3dd4, init = initials_m3,
                        end.time = end.time)

res_model3fd<-model3fd(pars = parameters_m3fd, init = initials_m3,
                        end.time = end.time)
res_model3fd2<-model3fd(pars = parameters_m3fd2, init = initials_m3,
                        end.time = end.time)
res_model3fd3<-model3fd(pars = parameters_m3fd3, init = initials_m3,
                        end.time = end.time)
res_model3fd4<-model3fd(pars = parameters_m3fd4, init = initials_m3,
                        end.time = end.time)

#> m4: HS
res_model4dd<-model4dd(pars = parameters_m4dd, init = initials_m4,
                       end.time = end.time)
res_model4dd2<-model4dd(pars = parameters_m4dd2, init = initials_m4,
                        end.time = end.time)
res_model4dd3<-model4dd(pars = parameters_m4dd3, init = initials_m4,
                        end.time = end.time)
res_model4dd4<-model4dd(pars = parameters_m4dd4, init = initials_m4,
                        end.time = end.time)
res_model4dd5<-model4dd(pars = parameters_m4dd5, init = initials_m4,
                        end.time = end.time)
res_model4dd6<-model4dd(pars = parameters_m4dd6, init = initials_m4,
                        end.time = end.time)
res_model4dd7<-model4dd(pars = parameters_m4dd7, init = initials_m4,
                        end.time = end.time)
res_model4dd8<-model4dd(pars = parameters_m4dd8, init = initials_m4,
                        end.time = end.time)

res_model4fd<-model4fd(pars = parameters_m4fd, init = initials_m4,
                       end.time = end.time)
res_model4fd2<-model4fd(pars = parameters_m4fd2, init = initials_m4,
                        end.time = end.time)
res_model4fd3<-model4fd(pars = parameters_m4fd3, init = initials_m4,
                        end.time = end.time)
res_model4fd4<-model4fd(pars = parameters_m4fd4, init = initials_m4,
                        end.time = end.time)
res_model4fd5<-model4fd(pars = parameters_m4fd5, init = initials_m4,
                        end.time = end.time)
res_model4fd6<-model4fd(pars = parameters_m4fd6, init = initials_m4,
                        end.time = end.time)
res_model4fd7<-model4fd(pars = parameters_m4fd7, init = initials_m4,
                        end.time = end.time)
res_model4fd8<-model4fd(pars = parameters_m4fd8, init = initials_m4,
                        end.time = end.time)

#> m5: LSD
res_model5dd<-model5dd(pars = parameters_m5dd, init = initials_m5,
                        end.time = end.time)
res_model5dd2<-model5dd(pars = parameters_m5dd2, init = initials_m5,
                        end.time = end.time)

res_model5fd<-model5fd(pars = parameters_m5fd, init = initials_m5,
                        end.time = end.time)
res_model5fd2<-model5fd(pars = parameters_m5fd2, init = initials_m5,
                        end.time = end.time)

#> m6: FMD
res_model6dd<-model6dd(pars = parameters_m6dd, init = initials_m6,
                       end.time = end.time)
res_model6dd$results<-res_model6dd$results %>% 
  relocate(M,.after = R) 

res_model6dd2<-model6dd(pars = parameters_m6dd2, init = initials_m6,
                        end.time = end.time)
res_model6dd2$results<-res_model6dd2$results %>% 
  relocate(M,.after = R) 

res_model6dd3<-model6dd(pars = parameters_m6dd3, init = initials_m6,
                        end.time = end.time)
res_model6dd3$results<-res_model6dd3$results %>% 
  relocate(M,.after = R) 

res_model6dd4<-model6dd(pars = parameters_m6dd4, init = initials_m6,
                        end.time = end.time)
res_model6dd4$results<-res_model6dd4$results %>% 
  relocate(M,.after = R) 

res_model6fd<-model6fd(pars = parameters_m6fd, init = initials_m6,
                       end.time = end.time)
res_model6fd$results<-res_model6fd$results %>% 
  relocate(M,.after = R) 

res_model6fd2<-model6fd(pars = parameters_m6fd2, init = initials_m6,
                        end.time = end.time)
res_model6fd2$results<-res_model6fd2$results %>% 
  relocate(M,.after = R) 

res_model6fd3<-model6fd(pars = parameters_m6fd3, init = initials_m6,
                        end.time = end.time)
res_model6fd3$results<-res_model6fd3$results %>% 
  relocate(M,.after = R)

res_model6fd4<-model6fd(pars = parameters_m6fd4, init = initials_m6,
                        end.time = end.time)
res_model6fd4$results<-res_model6fd4$results %>% 
  relocate(M,.after = R) 

#> m7: Brucellosis
res_model7dd<-model7dd(pars = parameters_m7dd, init = initials_m7,
                        end.time = end.time)
res_model7dd$results<-res_model7dd$results %>% 
  relocate(M,.after = R) 

res_model7dd2<-model7dd(pars = parameters_m7dd2, init = initials_m7,
                        end.time = end.time)
res_model7dd2$results<-res_model7dd2$results %>% 
  relocate(M,.after = R) 

res_model7fd<-model7fd(pars = parameters_m7fd, init = initials_m7,
                        end.time = end.time)
res_model7fd$results<-res_model7fd$results %>% 
  relocate(M,.after = R)

res_model7fd2<-model7fd(pars = parameters_m7fd2, init = initials_m7,
                        end.time = end.time)
res_model7fd2$results<-res_model7fd2$results %>% 
  relocate(M,.after = R)

############## 4) MULTIPLE RUNS - ALL MODELS ################################
# RUN the whole thing will take a lot of RAM and TIME, we can run this separately
# use the same parameters as the single run and n_rep simulations 

start <- print(Sys.time())

#> MX RUNS MODEL 1 - NO infection #####
sim_rep_m1 <- replicate(n_rep,(model1(pars = parameters_m1, init = initials_m1,
                                    end.time = end.time)))

#> MX RUNS MODEL 2 SI - Anthrax #####
sim_rep_m2dd<-replicate(n_rep,(model2dd(pars = parameters_m2dd, init = initials_m2,
                                        end.time = end.time)))
sim_rep_m2fd<-replicate(n_rep,(model2fd(pars = parameters_m2fd, init = initials_m2,
                                        end.time = end.time)))

#> MX RUNS MODEL 3 SEI - Bovine tuberculosis #####
sim_rep_m3dd<-replicate(n_rep,(model3dd(pars = parameters_m3dd, init = initials_m3,
                                         end.time = end.time)))
sim_rep_m3dd2<-replicate(n_rep,(model3dd(pars = parameters_m3dd2, init = initials_m3,
                                         end.time = end.time)))
sim_rep_m3dd3<-replicate(n_rep,(model3dd(pars = parameters_m3dd3, init = initials_m3,
                                         end.time = end.time)))
sim_rep_m3dd4<-replicate(n_rep,(model3dd(pars = parameters_m3dd4, init = initials_m3,
                                         end.time = end.time)))

sim_rep_m3fd<-replicate(n_rep,(model3fd(pars = parameters_m3fd, init = initials_m3,
                                         end.time = end.time)))
sim_rep_m3fd2<-replicate(n_rep,(model3fd(pars = parameters_m3fd2, init = initials_m3,
                                         end.time = end.time)))
sim_rep_m3fd3<-replicate(n_rep,(model3fd(pars = parameters_m3fd3, init = initials_m3,
                                         end.time = end.time)))
sim_rep_m3fd4<-replicate(n_rep,(model3fd(pars = parameters_m3fd4, init = initials_m3,
                                         end.time = end.time)))

#> MX RUNS MODEL 4 SIRS - Hemorrhagic septicemia #####
sim_rep_m4dd<-replicate(n_rep,(model4dd(pars = parameters_m4dd, init = initials_m4,
                                        end.time = end.time)))
sim_rep_m4dd2<-replicate(n_rep,(model4dd(pars = parameters_m4dd2, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4dd3<-replicate(n_rep,(model4dd(pars = parameters_m4dd3, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4dd4<-replicate(n_rep,(model4dd(pars = parameters_m4dd4, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4dd5<-replicate(n_rep,(model4dd(pars = parameters_m4dd5, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4dd6<-replicate(n_rep,(model4dd(pars = parameters_m4dd6, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4dd7<-replicate(n_rep,(model4dd(pars = parameters_m4dd7, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4dd8<-replicate(n_rep,(model4dd(pars = parameters_m4dd8, init = initials_m4,
                                         end.time = end.time)))

sim_rep_m4fd<-replicate(n_rep,(model4fd(pars = parameters_m4fd, init = initials_m4,
                                        end.time = end.time)))
sim_rep_m4fd2<-replicate(n_rep,(model4fd(pars = parameters_m4fd2, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4fd3<-replicate(n_rep,(model4fd(pars = parameters_m4fd3, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4fd4<-replicate(n_rep,(model4fd(pars = parameters_m4fd4, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4fd5<-replicate(n_rep,(model4fd(pars = parameters_m4fd5, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4fd6<-replicate(n_rep,(model4fd(pars = parameters_m4fd6, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4fd7<-replicate(n_rep,(model4fd(pars = parameters_m4fd7, init = initials_m4,
                                         end.time = end.time)))
sim_rep_m4fd8<-replicate(n_rep,(model4fd(pars = parameters_m4fd8, init = initials_m4,
                                         end.time = end.time)))

#> MX RUNS MODEL 5 SEIRS - Lumpy skin disease #####
sim_rep_m5dd<-replicate(n_rep,(model5dd(pars = parameters_m5dd, init = initials_m5,
                                         end.time = end.time)))
sim_rep_m5dd2<-replicate(n_rep,(model5dd(pars = parameters_m5dd2, init = initials_m5,
                                         end.time = end.time)))

sim_rep_m5fd<-replicate(n_rep,(model5fd(pars = parameters_m5fd, init = initials_m5,
                                         end.time = end.time)))
sim_rep_m5fd2<-replicate(n_rep,(model5fd(pars = parameters_m5fd2, init = initials_m5,
                                         end.time = end.time)))

#> MX RUNS MODEL 6 SEIRMS/E - Foot and mouth disease #####
sim_rep_m6dd<-replicate(n_rep,(model6dd(pars = parameters_m6dd, init = initials_m6, 
                                         end.time = end.time)))
sim_rep_m6dd2<-replicate(n_rep,(model6dd(pars = parameters_m6dd2, init = initials_m6,
                                         end.time = end.time)))
sim_rep_m6dd3<-replicate(n_rep,(model6dd(pars = parameters_m6dd3, init = initials_m6,
                                         end.time = end.time)))
sim_rep_m6dd4<-replicate(n_rep,(model6dd(pars = parameters_m6dd4, init = initials_m6,
                                         end.time = end.time)))

sim_rep_m6fd<-replicate(n_rep,(model6fd(pars = parameters_m6fd,  
                                        init = initials_m6,end.time = end.time)))
sim_rep_m6fd2<-replicate(n_rep,(model6fd(pars = parameters_m6fd2,
                                         init = initials_m6,end.time = end.time)))
sim_rep_m6fd3<-replicate(n_rep,(model6fd(pars = parameters_m6fd3,                          
                                         init = initials_m6,end.time = end.time)))
sim_rep_m6fd4<-replicate(n_rep,(model6fd(pars = parameters_m6fd4,                          
                                         init = initials_m6,end.time = end.time)))

#> MX RUNS MODEL 7 SEIRMS/E - Brucellosis #####
sim_rep_m7fd<-replicate(n_rep,(model7fd(pars = parameters_m7fd, init = initials_m7,
                                         end.time = end.time)))
sim_rep_m7fd2<-replicate(n_rep,(model7fd(pars = parameters_m7fd2, init = initials_m7,
                                         end.time = end.time)))

sim_rep_m7dd<-replicate(n_rep,(model7dd(pars = parameters_m7dd, init = initials_m7,
                                         end.time = end.time)))
sim_rep_m7dd2<-replicate(n_rep,(model7dd(pars = parameters_m7dd2, init = initials_m7,
                                         end.time = end.time)))

end <- print(Sys.time())

time_to_simulate <- end - start
time_to_simulate

######## 5) PLOTTING Single and Multiple runs ######## 
## model list ########
### single simulation ----
s <- list(res_model1,
          
          res_model2dd, res_model2fd,   #3
          
          res_model3dd,  res_model3fd,  #5
          res_model3dd2, res_model3fd2, #7
          res_model3fd3, res_model3dd3, #9
          res_model3fd4, res_model3dd4, #11
          
          res_model4dd,  res_model4fd,  #13
          res_model4dd2, res_model4fd2, #15
          res_model4fd3, res_model4dd3, #17
          res_model4fd4, res_model4dd4, #19
          res_model4dd5, res_model4fd5, #21
          res_model4dd6, res_model4fd6, #23
          res_model4fd7, res_model4dd7, #25
          res_model4fd8, res_model4dd8, #27
          
          res_model5dd,  res_model5fd,  #29
          res_model5dd2, res_model5fd2, #31
          
          res_model6dd,  res_model6fd,  #33
          res_model6dd2, res_model6fd2, #35
          res_model6fd3, res_model6dd3, #37
          res_model6fd4, res_model6dd4, #39
          
          res_model7dd,  res_model7fd, #41
          res_model7fd2, res_model7dd2) #43

### multiple simulation ----
sim_rep_m <- list(
          sim_rep_m1, 
          
          sim_rep_m2dd, sim_rep_m2fd,  
          
          sim_rep_m3dd,  sim_rep_m3fd,
          sim_rep_m3dd2, sim_rep_m3fd2,
          sim_rep_m3fd3, sim_rep_m3dd3,
          sim_rep_m3fd4, sim_rep_m3dd4,
          
          sim_rep_m4dd,  sim_rep_m4fd,
          sim_rep_m4dd2, sim_rep_m4fd2,
          sim_rep_m4fd3, sim_rep_m4dd3,
          sim_rep_m4fd4, sim_rep_m4dd4,
          sim_rep_m4dd5, sim_rep_m4fd5,
          sim_rep_m4dd6, sim_rep_m4fd6,
          sim_rep_m4fd7, sim_rep_m4dd7,
          sim_rep_m4fd8, sim_rep_m4dd8,
          
          sim_rep_m5dd,  sim_rep_m5fd,
          sim_rep_m5dd2, sim_rep_m5fd2,
          
          sim_rep_m6dd,  sim_rep_m6fd, #33 
          sim_rep_m6dd2, sim_rep_m6fd2,#35
          sim_rep_m6fd3, sim_rep_m6dd3,#37
          sim_rep_m6fd4, sim_rep_m6dd4,#39
          
          sim_rep_m7dd,  sim_rep_m7fd,  #41
          sim_rep_m7fd2, sim_rep_m7dd2) #43

# model's name  ----
nam<-c('no_infection',
       'Anthrax_DD_beta3e-5','Anthrax_FD_beta001',
       
       "bTB_DD_beta00014", "bTB_FD_beta04",
       "bTB_DD_beta27e-5", "bTB_FD_beta0008",  
       "bTB_FD_beta00014", "bTB_DD_beta4e-6",
       "bTB_FD_beta00063", "bTB_DD_beta2e-5",
       
       'HS_DD_beta03death005','HS_FD_beta100death005',
       'HS_DD_beta03death05', 'HS_FD_beta100death05',
       
       'HS_FD_beta03death005','HS_DD_beta0001death005',
       'HS_FD_beta03death05', 'HS_DD_beta0001death05',
       
       'HS_DD_beta05death005','HS_FD_beta165death005',
       'HS_DD_beta05death05', 'HS_FD_beta165death05',
       
       'HS_FD_beta05death005','HS_DD_beta0002death005',
       'HS_FD_beta05death05', 'HS_DD_beta0002death05',
       
       'LSD_DD_beta0008',  'LSD_FD_beta2',
       'LSD_FD_beta0032',  'LSD_DD_beta1e-4',
       
       'FMD_DD_beta0026', 'FMD_FD_beta7',
       'FMD_DD_beta21',   'FMD_FD_beta6552',
       'FMD_FD_beta0115', 'FMD_DD_beta4e-4',
       'FMD_FD_beta21',   'FMD_DD_beta007',
         
       'Brucellosis_DD_beta5e-6', 'Brucellosis_FD_beta00016',
       'Brucellosis_FD_beta5e-3', 'Brucellosis_DD_beta1e-5')

# Single: convert res_model to data.frame
for(i in 1:length(s)){
 s[[i]]<-s[[i]]$results %>%
 mutate(time_y = time/365) %>% #convert day to year for plotting
  melt(id.vars = c('time','time_y'),
         value.name = 'value', variable.name = 'class')
  
  s[[i]]$model <- paste0(nam[[i]]) # adding the model name
}

# FUNCTION adding time and run number for multiple model's simulations ----
pop_sim_prep <- function(x, n_rep, end.time, melt){ # x = simulation of model, e.g. sim_run_m1
  df<-list()
  for (i in 1:n_rep){
    run <- paste("run", seq(n_rep), sep="")
    names(df)[i]
    df[[i]] <- x[,i]$results[,-c(1)]
    df[[i]]$time_d <- seq(from=1,to=end.time+1,by=1)
  }
  
  df <- map2(df,run, ~cbind(.x, run = .y))   # adding number of replications to the column
  df2 <- data.table::rbindlist(df)          # binding row
  
  if (melt == T) {  #option for melting the data in case we need...
    
    df3 <- melt(df2, id.vars = c('time_d','run'))
    return(df3 = data.frame(df3))
  } 
  
  else  {
    return( df2 = data.frame(df2)) 
  }
  
}

# rearrange df and calculate total population change (%) --------
m <- list()
for (i in 1:length(sim_rep_m)) {
  m[[i]] <- pop_sim_prep(x = sim_rep_m[[i]], n_rep=n_rep, end.time= end.time, melt = F)
  m[[i]] <- m[[i]]%>%
    group_by(run) %>%
    mutate(Ndiff = ((last(N)-first(N))/first(N))*100)%>% #calculate change in the total population at year100, and year0
    mutate(time_y = time_d/365) %>% #convert day to year for plotting
    as.data.frame()
  m[[i]]$model <- paste0(nam[[i]])
  
}

# relocate the columns of FMD and Brucellosis model's compartments as it easier for plotting-labeling
# for m6 & m7
m67<-m[32:43] 

head(m67[[1]])

for (i in 1:length(m67)) {
  m67[[i]]<-m67[[i]]%>% 
    relocate(M,.after = R) 
} 

# check columns before plotting: the population class should be in an order like SIRN, SIN, SIERN, if not back to the relocate
for (i in 1:length(m67)){
  names(m67)[i]
  print(head(m67[[i]]))
}

# Modelling end here #########################

# Save files #########################
# We can skip this :: save the data frame  (.rds) for working next time
# for (i in 1:length(s)) { saveRDS(s[[i]], file = paste0("df_m",i,"_",nam[[i]],"_1run.rds")) }

# for (i in 1:length(m)) {saveRDS(m[[i]], file = paste0("df_m",i,"_",nam[[i]],"_100runs.rds")) }

#Anthrax = FD, bTB = FD , HS FD=DD are similar so select DD and present FD in the supp,
#LSD = FD as indirect contact, cow density might not relate with transmission

str(s)
#-----------------------------------------------------

# Start again here:
# 
# PLOTTING FOR Main text ####
# single run ######

# select model for main text
s1<-list(s[[1]], 
         s[[3]],
         s[[4]],
         s[[12]],
         s[[22]],
         s[[28]],
         s[[36]],
         s[[42]])

for (i in 1:length(s1)){
  names(s1)[i]
  print(head(s1[[i]]))
}

my_theme <- theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))

nam2<- c('No infection',
         'Anthrax FD',
         'bTB DD',
         'HS DD: fatality 0.53%',
         'HS DD: fatality 5.8%',
         'LSD DD',
         'FMD FD',
         'Brucellosis FD')
p <- list()  
for (i in 1:length(s1)) {
  # Condition
  model <- s1[[i]]$model
  # > p1 - plot 1 run no infection ########
  if (any(model == "no_infection")){
    p[[i]]<- ggplot(s1[[i]]) + 
      geom_line(aes(x = time_y ,y = value,  color = class))  +
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[[i]],") ", nam2[[i]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c( 'adult','subadult','calf','total' ),
                          values = c('a'='seagreen4',
                                     'sa'='firebrick',
                                     'c'='dodgerblue3',
                                     'N'='#153030'))+ 
    
    my_theme +
    theme(legend.position = "none")   
    # guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))
  } 
  else if (any(model == "Anthrax")) {
    #> p2 - plot 1 run SI Anthrax ######
    p[[i]] <- s1[[i]] %>%  filter(class %in% c("S","I","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[[i]],") ", nam2[[i]]))  +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','I', 'total') , #this one label is manual
                          values = c('S'='seagreen4',
                                     'I'='firebrick',
                                     'N'='#153030'))+
      my_theme +
      theme(legend.position = "none")   
    # guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))
    }
  else if(any(model == "bTB")) {
    #> p3 - plot 1 run SEI bTB ######
    p[[i]] <- s1[[i]] %>% filter(class %in% c("S","E","I", "N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[[i]],") ", nam2[[i]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','E','I','total' ),
                          values = c('S'='seagreen4',
                                     'E'='darkorange2',
                                     'I'='firebrick',
                                     "N"='#153030'))+ 
      my_theme +
      theme(legend.position = "none")   
    # guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))
  } 
  else if (any(model == "HS_DD_beta03death005")) {
    #> p4 plot 1 run SIRS HS ######
    p[[i]] <- s1[[i]] %>% filter(class %in% c("S","I","R","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[[i]],") ", nam2[[i]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','I','R','total' ),
                          values = c('S'='seagreen4',
                                     'I'='firebrick',
                                     "R"='dodgerblue3',
                                     "N"='#153030'))+ 
      my_theme +
      theme(legend.position = "none")   
    # guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))
    } 
  else if (any(model == "HS_05")) {
    #> p4 plot 1 run SIRS HS ######
    p[[i]] <- s1[[i]] %>% filter(class %in% c("S","I","R","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[[i]],") ", nam2[[i]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','I','R','total' ),
                          values = c('S'='seagreen4',
                                     'I'='firebrick',
                                     "R"='dodgerblue3',
                                     "N"='#153030'))+ 
      my_theme+
      theme(legend.position = "none") 
    # guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))
  }
  else if  (any(model =="LSD")) {
    #> p5 plot 1 run SEIR LSD ######
    p[[i]]<- s1[[i]] %>% filter(class %in% c("S","E","I","R","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[[i]],") ", nam2[[i]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','E','I',"R",'total' ),
                          values = c('S'='seagreen4',
                                     'E'='darkorange2',
                                     'I'='firebrick',
                                     "R"='dodgerblue3',
                                     "N"='#153030'))+ 
      my_theme +
      theme(legend.position = "none")   
    # guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))
  }  else if  (any(model == "FMD")) {
    #> p6 plot 1 run SEIRMS/E FMD ######
    p[[i]]<-s1[[i]] |>filter(class %in% c("S","E","I","R","M","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[[i]],") ", nam2[[i]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','E','I',"R",'M','total' ),#'total change (%)'),
                          values = c('S'='seagreen4',
                                     'E'='darkorange2',
                                     'I'='firebrick',
                                     "R"='dodgerblue3',
                                     "M"='mediumorchid4',
                                     "N"='#153030')) +
      my_theme +
      theme(legend.position = "none")   
    # guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))

  } else {
    #> p7 plot 1 run SEIRMS/E Brucellosis ######
    p[[i]]<-s1[[i]]%>% filter(class %in% c("S","E","I","R","M","N")) %>%
      ggplot() + 
      geom_line(aes(x = time_y, y = value, color = class))+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[[i]],") ", nam2[[i]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','E','I',"R",'M','total' ),#'total change (%)'),
                          values = c('S'='seagreen4',
                                     'E'='darkorange2',
                                     'I'='firebrick',
                                     "R"='dodgerblue3',
                                     "M"='mediumorchid4',
                                     "N"='#153030'))+
      my_theme +
      theme(legend.position = "none")   
    # guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))     
  }
}

print(p)

# plot multiple runs ###### 

str(m)

# select 100 simulations  ---------
m1<-list(m[[1]], 
         m[[3]],
         m[[4]],
         m[[12]],
         m[[22]],
         m[[28]],
         m67[[5]],
         m67[[11]])

for (i in 1:length(m1)){
  names(m1)[i]
  print(head(m1[[i]]))
}

pl <- list()
for (i in 1:length(m1)) {
  # Condition
  model <- m1[[i]]$model
  #  > pl1 - plot 100 runs no infection ########
  if (any(m1[[i]]$model == "no_infection")) {
    pl[[i]]<-ggplot(m1[[i]]) + 
      geom_line(aes(x = time_y ,y = a,  group = run, color = 'adult'), linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = sa, group = run, color = 'subadult'),linewidth = 0.1, alpha = 0.12) +
      geom_line(aes(x = time_y, y = c,  group = run, color = 'calf' ), linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = N,  group = run, color = 'total'), linewidth = 0.1, alpha = 0.12) +
      labs(x="Years", y= "Population", 
            title= paste0(LETTERS[9:16][i],") ", nam2[[i]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c( 'adult','subadult','calf','total'),
                          values = c('seagreen4',
                                     'firebrick',
                                     'dodgerblue3',
                                     '#153030'),
                          breaks = c('adult','subadult','calf','total'))+ 
      my_theme+
      
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth =0.7 )))+
      
      stat_summary(m1[[i]], mapping =aes( x = time_y, y = a, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(m1[[i]], mapping =aes( x = time_y, y = sa, group = 1), fun=mean, geom="line", colour= 'firebrick',linewidth = 0.5)+
      stat_summary(m1[[i]], mapping =aes( x = time_y, y = c, group = 1), fun=mean, geom="line", colour= 'dodgerblue3',linewidth = 0.5)+
      stat_summary(m1[[i]], mapping =aes( x = time_y, y = N, group = 1), fun=mean, geom="line", colour= '#153030',linewidth = 0.5)
  }
  else if (any(model == "Anthrax_FD_beta001")) {       
    #> pl2 - plot 100 runs SI Anthrax ######
    pl[[i]] <-ggplot(m1[[i]])+ 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'), linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'), linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ), linewidth = 0.1, alpha = 0.12)+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[9:16][i],") ", nam2[[i]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','I', 'total') , #this one label is manual
                          values = c('seagreen4',
                                     'firebrick',
                                     '#153030'),
                          breaks = c('S','I','total'))+
      my_theme+
      guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom ="line", colour="firebrick",linewidth = 0.5)
  }
  else if(any(model == "bTB_DD_beta00014")) {
    #> pl3 - plot 100 runs SEI bTB ######
    pl[[i]] <-ggplot(m1[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="years", y= "population", 
           title= paste0(LETTERS[9:16][i],") ", nam2[[i]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "population",
                          labels = c('S','E','I','total'),
                          values = c('seagreen4',
                                     'darkorange2',
                                     'firebrick',
                                     '#153030'),
                          breaks = c('S','E','I','total'))+ #blackgreen
      my_theme+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5) +
      stat_summary(m1[[i]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)
     #blackgreen
  }
      
  else if(any(model == "HS_DD_beta03death005")) {
    #> pl4 - plot 100 runs SIRS HS ######
    pl[[i]] <-ggplot(m1[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[9:16][i],") ", nam2[[i]]))  +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','I','R','total'),
                          values = c('seagreen4',
                                     'firebrick',
                                     'dodgerblue3',
                                     '#153030'),
                          breaks = c('S','I','R','total'))+ #blackgreen
      my_theme+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(m1[[i]], mapping = aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)
  }
  else if(any(model == "HS_DD_beta05death05")) {
    #> pl4 - plot 100 runs SIRS HS ######
    pl[[i]] <-ggplot(m1[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[9:16][i],") ", nam2[[i]]))  +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','I','R','total'),
                          values = c('seagreen4',
                                     'firebrick',
                                     'dodgerblue3',
                                     '#153030'),
                          breaks = c('S','I','R','total'))+ #blackgreen
      my_theme+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(m1[[i]], mapping = aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)
  }
  else if(any(model == "LSD_DD_beta0008")){
    #> pl5 - plot 100 runs SEIRS LSD ######
    pl[[i]] <-ggplot(m1[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[9:16][i],") ", nam2[[i]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','E','I','R','total'),
                          values = c('seagreen4',
                                     'darkorange2',
                                     'firebrick',
                                     'dodgerblue3',
                                     '#153030'),
                          breaks = c('S','E','I','R','total'))+ #blackgreen
      my_theme+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)#blackgreen
  }
  else if(any(model == "FMD_FD_beta0115"))  {
    #> pl6 - plot 100 runs SEIRMS/E FMD ######
    pl[[i]] <-ggplot(m1[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[9:16][i],") ", nam2[[i]]))  +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','E','I','R','M','total'),
                          values = c('seagreen4',
                                     'darkorange2',
                                     'firebrick',
                                     'dodgerblue3',
                                     'mediumorchid4',
                                     '#153030'),
                          breaks = c('S','E','I','R','M','total'))+ #blackgreen
      my_theme+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)
  }
  else {
    #> pl7 - plot 100 runs SEIRMS/E Brucellosis ######
    pl[[i]] <-ggplot(m1[[i]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="Years", y= "Population", 
           title= paste0(LETTERS[9:16][i],") ", nam2[[i]]))  +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','E','I','R','M','total'),
                          values = c('seagreen4',
                                     'darkorange2',
                                     'firebrick',
                                     'dodgerblue3',
                                     'mediumorchid4',
                                     '#153030'),
                          breaks = c('S','E','I','R','M','total'))+ #blackgreen
      my_theme+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
      stat_summary(m1[[i]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)
  }
}   
print(pl)

# Fig. 1 : Combine single and 100 simulations plots  #######
# use 'Patchwork'
pp1 <- wrap_plots(p, ncol=1) & plot_annotation(title = '1 simulation') & theme(plot.title = element_text(hjust = 0.01))
pp2 <- wrap_plots(pl, ncol=1) & plot_annotation(title = '100 simulations') & theme(plot.title = element_text(hjust = 0.01))
pp3<-(wrap_elements(pp1)|wrap_elements(pp2))
pp4<-wrap_elements(pp3 + plot_layout(widths = c(0.76,0.98)))
pp4

ggsave("patchwork_allmodels8.png",pp4, width = 27, height = 44, units = 'cm', dpi = 600)

# Fig. 2 : compare LSD & FMD & Brucellosis ######
# separate small panel: LSD, FMD & Brucellosis DD, FD & DD* (rescale) 

# > Select the models from the list
m2<-list(m[[28]],m[[30]],m[[31]], #LSD: DD beta = 0.008, FD beta = 0.032, DD rescale 0.032/N
         m[[32]],m[[36]],m[[37]], #FMD: DD beta = 0.026, FD beta = 0.115, DD rescale 0.115/N
         m[[40]],m[[42]],m[[43]]) #Brucellosis: DD beta = 5e-6, FD beta = 5e-3, DD rescale 5e-3/N

head(m2)         
#Fig 2,3,4
nam_fig2 <- list("LSD DD, beta = 0.008",
                 "LSD FD, beta = 0.032",
                 "LSD DD (rescale), beta = 1e-4", 
                 
                 "FMD DD, beta = 0.026",
                 "FMD FD, beta = 0.115",
                 "FMD DD (rescale), beta = 3e-4", 
                 
                 'Brucellosis DD, beta = 5.5e-6',
                 'Brucellosis FD, beta = 5.5e-3',
                 "Brucellosis DD (rescale), beta = 1e-5") 


lsd_dd<-ggplot(m2[[1]]) + geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
    geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
    labs(x="Years", y= "Population", 
         title= (nam_fig2[[1]])) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c('S','E','I','R','total'),
                        values = c('seagreen4',
                                   'darkorange2',
                                   'firebrick',
                                   'dodgerblue3',
                                   '#153030'),
                        breaks = c('S','E','I','R','total'))+ #blackgreen
    my_theme+
    theme(legend.position = "none")+
    stat_summary(m2[[1]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
    stat_summary(m2[[1]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
    stat_summary(m2[[1]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
    stat_summary(m2[[1]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
    stat_summary(m2[[1]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)

lsd_fd<-ggplot(m2[[2]]) + geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="Years", y= "Population", 
       title= (nam_fig2[[2]])) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "Compartment",
                      labels = c('S','E','I','R','total'),
                      values = c('seagreen4',
                                 'darkorange2',
                                 'firebrick',
                                 'dodgerblue3',
                                 '#153030'),
                      breaks = c('S','E','I','R','total'))+ #blackgreen
  my_theme+ 
  theme(legend.position = "none")+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
  stat_summary(m2[[2]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
  stat_summary(m2[[2]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m2[[2]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
  stat_summary(m2[[2]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m2[[2]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)

lsd_dd_re<-ggplot(m2[[3]]) + geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="Years", y= "Population", 
       title= (nam_fig2[[3]])) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "Compartment",
                      labels = c('S','E','I','R','total'),
                      values = c('seagreen4',
                                 'darkorange2',
                                 'firebrick',
                                 'dodgerblue3',
                                 '#153030'),
                      breaks = c('S','E','I','R','total'))+ #blackgreen
  my_theme+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
  stat_summary(m2[[3]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
  stat_summary(m2[[3]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m2[[3]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
  stat_summary(m2[[3]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m2[[3]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)


fmd_dd <-ggplot(m2[[4]]) + 
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="Years", y= "Population", 
       title= (nam_fig2[[4]])) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "Compartment",
                      labels = c('S','E','I','R','M','total'),
                      values = c('seagreen4',
                                 'darkorange2',
                                 'firebrick',
                                 'dodgerblue3',
                                 'mediumorchid4',
                                 '#153030'),
                      breaks = c('S','E','I','R','M','total'))+ #blackgreen
  my_theme+
  theme(legend.position = "none")+
  stat_summary(m2[[4]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
  stat_summary(m2[[4]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m2[[4]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
  stat_summary(m2[[4]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m2[[4]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
  stat_summary(m2[[4]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)

fmd_fd <-ggplot(m2[[5]]) + 
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="Years", y= "Population", 
       title= (nam_fig2[[5]])) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "Compartment",
                      labels = c('S','E','I','R','M','total'),
                      values = c('seagreen4',
                                 'darkorange2',
                                 'firebrick',
                                 'dodgerblue3',
                                 'mediumorchid4',
                                 '#153030'),
                      breaks = c('S','E','I','R','M','total'))+ #blackgreen
  my_theme+
  theme(legend.position = "none")+
  stat_summary(m2[[5]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
  stat_summary(m2[[5]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m2[[5]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
  stat_summary(m2[[5]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m2[[5]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
  stat_summary(m2[[5]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)

fmd_fd_re <-ggplot(m2[[6]]) + 
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="Years", y= "Population", 
       title= (nam_fig2[[6]])) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "Compartment",
                      labels = c('S','E','I','R','M','total'),
                      values = c('seagreen4',
                                 'darkorange2',
                                 'firebrick',
                                 'dodgerblue3',
                                 'mediumorchid4',
                                 '#153030'),
                      breaks = c('S','E','I','R','M','total'))+ #blackgreen
  my_theme+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
  stat_summary(m2[[6]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
  stat_summary(m2[[6]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m2[[6]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
  stat_summary(m2[[6]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m2[[6]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
  stat_summary(m2[[6]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)

bru_dd <-ggplot(m2[[7]]) + 
      geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
      geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
      geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
      labs(x="Years", y= "Population", 
           title= (nam_fig2[[7]])) +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "Compartment",
                          labels = c('S','E','I','R','M','total'),
                          values = c('seagreen4',
                                     'darkorange2',
                                     'firebrick',
                                     'dodgerblue3',
                                     'mediumorchid4',
                                     '#153030'),
                          breaks = c('S','E','I','R','M','total'))+ #blackgreen
      my_theme+
      theme(legend.position = "none")+
      stat_summary(m2[[7]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
      stat_summary(m2[[7]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
      stat_summary(m2[[7]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
      stat_summary(m2[[7]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
      stat_summary(m2[[7]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
      stat_summary(m2[[7]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)
  
bru_fd <-ggplot(m2[[8]]) + 
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="Years", y= "Population", 
       title= (nam_fig2[[8]])) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "Compartment",
                      labels = c('S','E','I','R','M','total'),
                      values = c('seagreen4',
                                 'darkorange2',
                                 'firebrick',
                                 'dodgerblue3',
                                 'mediumorchid4',
                                 '#153030'),
                      breaks = c('S','E','I','R','M','total'))+ #blackgreen
  my_theme+
  theme(legend.position = "none")+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
  stat_summary(m2[[8]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
  stat_summary(m2[[8]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m2[[8]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
  stat_summary(m2[[8]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m2[[8]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
  stat_summary(m2[[8]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)

bru_fd_re <-ggplot(m2[[9]]) + 
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
  labs(x="Years", y= "Population", 
       title= (nam_fig2[[9]])) +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "Compartment",
                      labels = c('S','E','I','R','M','total'),
                      values = c('seagreen4',
                                 'darkorange2',
                                 'firebrick',
                                 'dodgerblue3',
                                 'mediumorchid4',
                                 '#153030'),
                      breaks = c('S','E','I','R','M','total'))+ #blackgreen
  my_theme+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
  stat_summary(m2[[9]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
  stat_summary(m2[[9]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
  stat_summary(m2[[9]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
  stat_summary(m2[[9]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
  stat_summary(m2[[9]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
  stat_summary(m2[[9]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)


p1<-(lsd_dd|lsd_fd|lsd_dd_re)
p2<-(fmd_dd|fmd_fd|fmd_fd_re)
p3<-(bru_dd|bru_fd|bru_fd_re)

w <- wrap_plots(p1/p2/p3) + plot_annotation(tag_levels = 'A') & 
  theme (plot.tag = element_text(face = 'bold'))
w

#ggsave("patchwork_lsd_fmd_bru2.png",w, width = 35, height = 24, units = 'cm', dpi = 400)

#> summarise basic stats #######  
# single run
for (i in 1:length(s)) {
  print (s[[i]] |> group_by(class,model) |>
          dplyr::summarise(Max=max(value),
                           Min=min(value)))
  }

# multiple runs
# add disease name
nam_dis<-c("No infection", 
           "Anthrax","Anthrax",
           "bTB","bTB","bTB","bTB",
           "bTB","bTB","bTB","bTB",
           "HS","HS","HS","HS", "HS","HS","HS","HS",
           "HS","HS","HS","HS", "HS","HS","HS","HS",
           "LSD","LSD","LSD","LSD",
           "FMD","FMD","FMD","FMD",
           "FMD","FMD","FMD","FMD",
           "Brucellosis","Brucellosis","Brucellosis","Brucellosis")

mx<-list()
for (i in 1:length(m)){
  
  mx[[i]] <- melt(m[[i]], id.vars = c('time_y','time_d','run','model',"Ndiff"),
                    value.name = 'value', variable.name = 'class')
  mx[[i]]$model2<- nam_dis[[i]]
  
}

for (i in 1:length(mx)){
  print(mx[[i]] |> group_by(class,model) |>
           dplyr::summarise(Mean = mean(value),Max=max(value), Min=min(value)))
}

# combining multiple runs df list
mx2 <-data.table::rbindlist(mx)
str(mx2)
#select some columns: Ndiff,run,model
dft <- mx2 %>% 
  dplyr::select(Ndiff,run,model,model2) %>% 
  drop_na() %>%
  distinct() %>%
  data.frame()

View(dft)

#calculating mean,median by models
dft2<-dft|>
  group_by(model)|>
  summarise(Mean = mean(Ndiff))|>
  arrange(desc(Mean))     

print(dft2,n=43) #n = 43, number of models

write.csv(dft, "df_ndiff_allmodels.csv")
write.csv(dft2, "df_ndiff_mean.csv")

## load df % population change (files in GitHub repo) ----------------------------------------------------------------------
dft <- read.csv("./df_ndiff_allmodels.csv")
str(dft)

table(dft$run) # check this
table(dft$model)
# Fig. 3 Boxplot #####
# select the models
want <- c('no_infection','Anthrax_DD_beta3e-5','Anthrax_FD_beta001', 
          'HS_FD_beta03death005', 'HS_DD_beta05death05',
          'bTB_DD_beta00014' , 'bTB_FD_beta00063',
          'LSD_DD_beta0008' , 'LSD_FD_beta0032', 'LSD_DD_beta1e-4',
          'FMD_DD_beta0026', 'FMD_FD_beta0115','FMD_DD_beta4e-4',
          'Brucellosis_DD_beta5e-6', 'Brucellosis_FD_beta5e-3','Brucellosis_DD_beta1e-5')

dft3<-dft |> filter(model %in% want)
dft3
dft3$model3<- recode_factor(dft3$model,
                            'no_infection'   = 'No infection',
                            'LSD_DD_beta0008' = 'LSD DD',
                            'LSD_FD_beta0032' = 'LSD FD',
                            'LSD_DD_beta1e-4' = 'LSD DD*', 
                            
                            'FMD_DD_beta0026' = 'FMD DD',
                            'FMD_FD_beta0115' = 'FMD FD',
                            'FMD_DD_beta4e-4' =  'FMD DD*',
                            
                            'Brucellosis_DD_beta5e-6' = 'Brucellosis DD',
                            'Brucellosis_FD_beta5e-3' = 'Brucellosis FD',
                            'Brucellosis_DD_beta1e-5' = 'Brucellosis DD*',
                            
                            'Anthrax_DD_beta3e-5' = 'Anthrax DD*',
                            'Anthrax_FD_beta001'  =  'Anthrax FD',
                            
                            'HS_DD_beta05death05' ='HS DD',
                            'HS_FD_beta03death005' = 'HS FD', #mortality 5.8%
                            
                            'bTB_DD_beta00014' = 'bTB DD', 
                            'bTB_FD_beta00063' =  'bTB FD')
dft3$model2<-as.factor(dft3$model2)
dft3$model3<-as.factor(dft3$model3)

table(dft3$model2)
table(dft3$model3)

# Manually reorder the disease factor
dft3$model2 <- factor(dft3$model2, levels = 
                        c( "No infection",
                           "Anthrax",
                           "HS",
                           "bTB",
                           "LSD",
                           "FMD",
                           "Brucellosis"))

#Anthrax,HS,bTB, LSD,FMD,Bru
# one box per variety
pbase<-dft3%>%
  ggplot( aes(x=model3, y=Ndiff)) + 
  geom_boxplot() +
  facet_wrap(~model2, scale="free_x")+ 
  coord_cartesian(ylim = c(-120, 500))+
  ggtitle("Gaur population change by infectious diseases")+
  xlab("")+
  ylab("Population change (%)")+
  scale_y_continuous(limits = c(-200, 500), 
                     breaks = seq(from = -200, to = 500, by = 100)) +
  theme(plot.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))
pbase
#ggsave("gaur_ndiff_test_base.png",pbase,width = 40, height = 30, units = 'cm', dpi = 400)

#### 1) ggstat facet  

#This ggostat code may not be compatible with different R versions. ><
#Anthrax,HS,bTB, LSD,FMD,Bru
nam_bp<-c("LSD",
          "FMD",
          "Brucellosis",
          "Anthrax",
          "HS",
          "bTB")
dft_dis<-dft3 |> filter(!model2 %in% "no infection")
table(dft_dis$model2)

no <- dft3 |> filter(model2 %in% c("No infection")) 
table(no$model2)
str(no)
plt<-list()
for (i in 1:length(nam_bp)){
  plt[[i]] <- dft_dis %>% filter(model2 %in% c(nam_bp[[i]])) %>%
    grouped_ggbetweenstats(
      x=model3,
      y=Ndiff,
      grouping.var = model2,
      k=0,
      plot.type = "boxviolin",
      pairwise.comparisons=F,
      bf.message = F,
      results.subtitle = FALSE,
      centrality.point.args = list(size = 3, color = "darkred"),
      centrality.label.args = list(size = 4, 
                                   nudge_x = 0.5, 
                                   nudge_y = 0.5, 
                                   segment.linetype = 4, 
                                   min.segment.length = 0),
      xlab = "",
      ylab = "Population change (%)")+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    ggplot2::scale_y_continuous(limits = c(-200, 500), 
                                breaks = seq(from = -200, to = 500, by = 100)) +
    scale_colour_manual(values = c("lightcoral",
                                   "steelblue1",
                                   '#D9C756')) +
    theme(plot.title = element_text(size = 17),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          strip.text = element_text(size = 15))
}

plt_no <- no %>% 
  ggbetweenstats(
    x=model2,
    y=Ndiff,
    k=0,
    plot.type = "boxviolin",
    pairwise.comparisons=F,
    bf.message = F,
    results.subtitle = FALSE,
    centrality.point.args = list(size = 3, color = "darkred"),
    centrality.label.args = list(size = 4, 
                                 nudge_x = 0.5, 
                                 nudge_y = 0.5, 
                                 segment.linetype = 4, 
                                 min.segment.length = 0),
    title= "No infection",
    xlab = "",
    ylab = "Population change (%)")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))+
  scale_colour_manual(values = c("#44A57CFF")) +
  theme(plot.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15))

plt_no

plt1<-wrap_elements((plt[[1]]|plt[[2]]|plt[[3]]) & plot_annotation(subtitle = "B)") & 
                                                      theme(title = element_text(size = 18)))
plt2<-wrap_elements((plt_no|plt[[4]]|plt[[5]]|plt[[6]]) & 
                      plot_annotation(title = "Gaur population change by infectious diseases",
                                      subtitle = "A)") & theme(title = element_text(size = 18),
                                                             plot.title = element_text(face = 'bold')))     
pbox<-(plt2/plt1)
pbox

#ggsave("gaur_ndiff_test.png",pbox,width = 50, height = 28, units = 'cm', dpi = 600)

# Done for the main text ...

# Figures in the supplementary materials --------------------------------------------
## Fig. Anthrax ####
#'Anthrax_FD_beta01','Anthrax_DD_beta3e-5'
a1 <- list(s[[2]], s[[3]])
a2 <- list(m[[2]], m[[3]])

for (i in 1:length(a1)) {
  print(head(a1[[i]]))
  print(head(a2[[i]]))
}

nam_an<-c("Anthrax FD, beta = 0.01",
          "Anthrax DD*, beta = 3e-5")

an1<-list()

for (i in 1:length(a1)) {
  an1[[i]] <- a1[[i]] %>%  filter(class %in% c("S","I","N")) %>%
    ggplot() +
    geom_line(aes(x = time_y, y = value, color = class))+
    labs(x="Years", y= "Population",
         title= nam_an[[i]]) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c('S','I', 'total') , #this one label is manual
                        values = c('S'='seagreen4',
                                   'I'='firebrick',
                                   'N'='#153030'))+
    my_theme +
    theme(legend.position = "none")
}

print(an1)

an2<-list()

for (i in 1:length(a2)) {
  an2[[i]]<-ggplot(a2[[i]])+ 
    
    geom_line(aes(x = time_y ,y = S, group = run, color = 'S'), linewidth = 0.1, alpha = 0.12) + 
    geom_line(aes(x = time_y, y = I, group = run, color = 'I' ), linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = N, group = run, color = 'total'), linewidth = 0.1, alpha = 0.12)+
    labs(x="Years", y= "Population", title=  nam_an[[i]])  +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c('S','I', 'total') , #this one label is manual
                        values = c('seagreen4',
                                   'firebrick',
                                   '#153030'),
                        breaks = c('S','I','total'))+
    my_theme +
    guides(color = guide_legend(override.aes = list(alpha = 1,linewidth=0.7)))+
    stat_summary(a2[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
    stat_summary(a2[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
    stat_summary(a2[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom ="line", colour="firebrick",linewidth = 0.5)
}

print(an2)

#anthrax combine plots

p_an1<-wrap_elements(wrap_plots(an1)+
                       plot_annotation(title = "1 simulation",
                                       tag_levels = "A"))
p_an2<-wrap_elements(wrap_plots(an2)+
                       plot_annotation(title = "100 simulstions",
                                       tag_levels = "A")+
                       plot_layout(guides = "collect") &
                       theme(legend.position = "bottom"))

p_an3<-p_an1/p_an2+plot_layout(heights = c(0.9,1))
p_an3
#ggsave("p_an32.png",p_an3,width = 20, height = 17, units = 'cm', dpi = 400)

## Fig. bTB ####
# "bTB_DD_beta00014", "bTB_FD_beta04"
#"bTB_DD_beta27e-5", "bTB_FD_beta0008"  
#"bTB_FD_beta00014", "bTB_DD_beta4e-6"
#"bTB_FD_beta00063", "bTB_DD_beta2e-5"

btb1<-(s[4:11]) 
btb2<-(m[4:11])
for (i in 1:length(btb1)) {
  print(head(btb1[[i]]))
  print(head(btb2[[i]]))
}
nam_btb<-c("bTB DD, beta = 0.0014",
           "bTB FD*, beta = 0.4",
           
           "bTB DD beta = 2.7e-5",
           "bTB FD*, beta = 0.008",
           
           "bTB FD, beta = 0.0014",
           "bTB DD*, beta = 4e-6",
           
           "bTB FD, beta = 0.0063",
           "bTB DD*, beta = 2e-5")

b1<-list()

for (i in 1:length(btb1)) {
  b1[[i]] <- btb1[[i]] %>% filter(class %in% c("S","E","I", "N")) %>%
    ggplot() +
    geom_line(aes(x = time_y, y = value, color = class))+
    labs(x="Years", y= "Population", title =  nam_btb[[i]]) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c('S','E','I','total' ),
                        values = c('S'='seagreen4',
                                   'E'='darkorange2',
                                   'I'='firebrick',
                                   "N"='#153030'))+
    my_theme 
  
}
#test
print(b1[[1]])

b2<-list()
for (i in 1:length(btb2)) {
  b2[[i]] <-ggplot(btb2[[i]]) +
    geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) +
    geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
    labs(x="Years", y= "Population", title = nam_btb[[i]])  +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c('S','E','I','total'),
                        values = c('seagreen4',
                                   'darkorange2',
                                   'firebrick',
                                   '#153030'),
                        breaks = c('S','E','I','total'))+ #blackgreen
    my_theme+
    theme(legend.position = "bottom")+
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
    stat_summary(btb2[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
    stat_summary(btb2[[i]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
    stat_summary(btb2[[i]], mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
    stat_summary(btb2[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)
}
#test
b2[[1]]

#Bovine TB combine plots
p_btb1<-wrap_plots(b1, ncol=2) +
  plot_annotation(title = '1 simulation') +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_btb1


p_btb2<-wrap_plots(b2, ncol=2) +
  plot_annotation(title = '100 simulations') +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_btb2

#ggsave("s_btb_single.png",p_btb1, width = 25, height = 30, units = 'cm', dpi = 400)
#ggsave("s_btb_multiple.png",p_btb2, width = 25, height = 30, units = 'cm', dpi = 400)

## Fig. HS ####

#'HS_DD_beta03death005','HS_FD_beta100death005'
#'HS_DD_beta03death05', 'HS_FD_beta100death05'
#'HS_FD_beta03death005','HS_DD_beta0001death005'
#'HS_FD_beta03death05', 'HS_DD_beta0001death05'

#'HS_DD_beta05death005','HS_FD_beta165death005'
#'HS_DD_beta05death05', 'HS_FD_beta165death05'
#'HS_FD_beta05death005','HS_DD_beta0002death005'
#'HS_FD_beta05death05', 'HS_DD_beta0002death05'

hs1<-s[12:27]
hs2<-m[12:27]

for (i in 1:length(hs1)) {
  print(head(hs1[[i]]))
  print(head(hs2[[i]]))
}

nam_hs <- c("HS DD, beta = 0.3, fatality = 0.5%",   
            "HS FD*, beta = 100, fatality = 0.5%",
            
            "HS DD, beta = 0.3, fatality = 5.8%", 
            "HS FD*, beta = 100, fatality = 5.8%",
                
            "HS FD, beta = 0.3, fatality = 0.5%", 
            "HS DD*, beta = 0.001, fatality = 0.5%", 
                
            "HS FD, beta = 0.3, fatality = 5.8%", 
            "HS DD*, beta = 0.001, fatality = 5.8%",
                
            "HS DD, beta = 0.5, fatality = 0.5%", 
            "HS FD*, beta = 165, fatality = 0.5% ",
                
            "HS DD, beta = 0.5, fatality = 5.8%", 
            "HS FD*, beta = 165, fatality = 5.8%",
                
            "HS FD, beta = 0.5, fatality = 0.5%", 
            "HS DD*, beta = 0.002, fatality = 0.5%", 
                
            "HS FD, beta = 0.5, fatality = 5.8%", 
            "HS DD*, beta = 0.002, fatality = 5.8%" )

h1<-list()
for (i in 1:length(hs1)) {
  h1[[i]] <- hs1[[i]] %>% filter(class %in% c("S","I","R","N")) %>%
    ggplot() +
    geom_line(aes(x = time_y, y = value, color = class))+
    labs(x="Years", y= "Population", title= paste0(nam_hs[[i]])) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c('S','I','R','total' ),
                        values = c('S'='seagreen4',
                                   'I'='firebrick',
                                   "R"='dodgerblue3',
                                   "N"='#153030'))+
    my_theme +     
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))
}

h2<-list()
for (i in 1:length(hs2)){
  h2[[i]] <-ggplot(hs2[[i]]) + 
    geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
    geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
    labs(x="Years", y= "Population", title= nam_hs[[i]]) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c('S','I','R','total'),
                        values = c('seagreen4',
                                   'firebrick',
                                   'dodgerblue3',
                                   '#153030'),
                        breaks = c('S','I','R','total'))+ #blackgreen
    my_theme +
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
    stat_summary(hs2[[i]], mapping = aes( x = time_y, y  = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
    stat_summary(hs2[[i]], mapping = aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
    stat_summary(hs2[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
    stat_summary(hs2[[i]], mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)
}
p_hs1<-wrap_elements(wrap_plots(h1, ncol=4)+
                       plot_annotation(title = "1 simulation",
                                       tag_levels = "A")+
                       plot_layout(guides = "collect") &
                       theme(legend.position = "bottom"))

p_hs2<-wrap_elements(wrap_plots(h2, ncol=4)+
                       plot_annotation(title = "100 simulations",
                                       tag_levels = "A")+
                       plot_layout(guides = "collect") &
                       theme(legend.position = "bottom"))
p_hs1
p_hs2

#ggsave("s_hs_single.png",p_hs1,width = 40, height = 25, units = 'cm', dpi = 400)
#ggsave("s_hs_multiple.png",p_hs2,width = 40, height = 25, units = 'cm', dpi = 400)

## Fig. LSD ####
#'LSD_DD_beta0008',  'LSD_FD_beta2',
#'LSD_FD_beta0032',  'LSD_DD_beta1e-4',

lsd1<-s[28:31]
lsd2<-m[28:31]

for (i in 1:length(lsd1)) {
  print(head(lsd1[[i]]))
  print(head(lsd2[[i]]))
}

nam_lsd<-c("LSD DD, beta = 0.008",
           "LSD FD*, beta = 2.4",
           
           "LSD FD, beta = 0.032",
           "LSD DD*, beta = 0.0001")
l1<-list()
for (i in 1:length(lsd1)) {
  l1[[i]]<- lsd1[[i]] %>% filter(class %in% c("S","E","I","R","N")) %>%
    ggplot() + 
    geom_line(aes(x = time_y, y = value, color = class))+
    labs(x="Years", y= "Population", 
         title=  nam_lsd[[i]]) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c('S','E','I',"R",'total' ),
                        values = c('S'='seagreen4',
                                   'E'='darkorange2',
                                   'I'='firebrick',
                                   "R"='dodgerblue3',
                                   "N"='#153030'))+ 
    my_theme +
    theme(legend.position = "none")   
}

l2<-list()
for(i in 1:length(lsd2)) {
  
  l2[[i]] <-ggplot(lsd2[[i]]) + 
    geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
    geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
    labs(x="Years", y= "Population", 
         title=  nam_lsd[[i]]) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c('S','E','I','R','total'),
                        values = c('seagreen4',
                                   'darkorange2',
                                   'firebrick',
                                   'dodgerblue3',
                                   '#153030'),
                        breaks = c('S','E','I','R','total'))+ #blackgreen
    my_theme + 
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
    stat_summary(lsd2[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+
    stat_summary(lsd2[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
    stat_summary(lsd2[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
    stat_summary(lsd2[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
    stat_summary(lsd2[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)#blackgreen
}

#LSD combine plots
p_l1<-wrap_elements(wrap_plots(l1)+
                      plot_annotation(title = "1 simulation",
                                      tag_levels = "A"))

p_l2<-wrap_elements(wrap_plots(l2)+
                      plot_annotation(title = "100 simulstions",
                                      tag_levels = "A")+
                      plot_layout(guides = "collect") &
                      theme(legend.position = "bottom"))

p_l3<-p_l1/p_l2+plot_layout(heights = c(0.9,1))
p_l3
#ggsave("s_lsd.png",p_l3,width = 25, height = 30, units = 'cm', dpi = 400)

## Fig. FMD ####

#'FMD_DD_beta0026', 'FMD_FD_beta7'
#'FMD_DD_beta21',   'FMD_FD_beta6552'
#'FMD_FD_beta0115', 'FMD_DD_beta4e-4'
#'FMD_FD_beta21',   'FMD_DD_beta007'

fmd1 <- s[32:39]
fmd2 <- m[32:39]

for (i in 1:length(fmd1)) {
  print(head(fmd1[[i]]))
  print(head(fmd2[[i]]))
}

nam_fmd<-c("FMD DD, beta = 0.026",
           "FMD FD*, beta = 7.8",
           
           "FMD DD, beta = 21",
           "FMD FD*, beta = 6552",
           
           "FMD FD, beta = 0.115",
           "FMD DD*, beta = 0.0004",
           
           "FMD FD, beta = 21",
           'FMD DD*, beta = 0.073')

f1<-list()
for (i in 1:length(fmd1)) { 
  f1[[i]]<-fmd1[[i]] |>filter(class %in% c("S","E","I","R","M","N")) %>%
    ggplot() + 
    geom_line(aes(x = time_y, y = value, color = class))+
    labs(x="Years", y= "Population", title = nam_fmd[[i]]) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c('S','E','I',"R",'M','total' ),#'total change (%)'),
                        values = c('S'='seagreen4',
                                   'E'='darkorange2',
                                   'I'='firebrick',
                                   "R"='dodgerblue3',
                                   "M"='mediumorchid4',
                                   "N"='#153030')) +
    my_theme +
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))
  
} 

f2<-list()
for (i in 1:length(fmd2)) { 
  f2[[i]] <-ggplot(fmd2[[i]]) + 
    geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
    geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
    labs(x="Years", y= "Population", nam_fmd[[i]]) +
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
    my_theme +
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.7)))+
    stat_summary(fmd2[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
    stat_summary(fmd2[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
    stat_summary(fmd2[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
    stat_summary(fmd2[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
    stat_summary(fmd2[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
    stat_summary(fmd2[[i]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)
}

p_f1<-wrap_elements(wrap_plots(f1,ncol=2)+
                      plot_annotation(title = "1 simulation",
                                      tag_levels = "A")+
                      plot_layout(guides = "collect") &
                      theme(legend.position = "bottom"))
p_f2<-wrap_elements(wrap_plots(f2,ncol=2)+
                      plot_annotation(title = "100 simulstions",
                                      tag_levels = "A")+
                      plot_layout(guides = "collect") &
                      theme(legend.position = "bottom"))
p_f1
p_f2
ggsave("s_fmd_single_rerun.png",p_f1,width = 25, height = 30, units = 'cm', dpi = 400)
#ggsave("s_fmd_multiple.png",p_f2,width = 25, height = 30, units = 'cm', dpi = 400)

## Fig. Brucellosis ####

# 'Brucellosis_DD_beta5e-6', 'Brucellosis_FD_beta00016',
# 'Brucellosis_FD_beta5e-3', 'Brucellosis_DD_beta1e-5'

bru1 <- s[40:43] 
bru2 <- m[40:43]

for (i in 1:length(bru1)) {
  print(head(bru1[[i]]))
  print(head(bru2[[i]]))
}

nam_bru<-c("Brucellosis DD, beta = 5e-6",
           "Brucellosis FD*, beta = 0.0016",
           
           "Brucellosis FD, beta = 0.005",
           "Brucellosis DD*, beta = 1e-5")

br1<-list()

for (i in 1:length(bru1)) {
  br1[[i]]<-bru1[[i]]%>% filter(class %in% c("S","E","I","R","M","N")) %>%
    ggplot() + 
    geom_line(aes(x = time_y, y = value, color = class))+
    labs(x="Years", y= "Population", 
         title= nam_bru[[i]]) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c('S','E','I',"R",'M','total' ),#'total change (%)'),
                        values = c('S'='seagreen4',
                                   'E'='darkorange2',
                                   'I'='firebrick',
                                   "R"='dodgerblue3',
                                   "M"='mediumorchid4',
                                   "N"='#153030'))+
    my_theme +
    theme(legend.position = "none")   
}

br2<-list()
for (i in 1:length(bru2)) {
  br2[[i]] <-ggplot(bru2[[i]]) + 
    geom_line(aes(x = time_y, y = N, group = run, color = 'total'),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),linewidth = 0.1, alpha = 0.12) + 
    geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = M, group = run, color = 'M' ),linewidth = 0.1, alpha = 0.12)+
    labs(x="Years", y= "Population", 
         title=  nam_bru[[i]]) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
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
    stat_summary(bru2[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean, geom="line", colour="#153030",linewidth = 0.5)+ #blackgreen
    stat_summary(bru2[[i]], mapping = aes(x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
    stat_summary(bru2[[i]], mapping = aes(x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",linewidth = 0.5)+
    stat_summary(bru2[[i]], mapping = aes(x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",linewidth = 0.5)+
    stat_summary(bru2[[i]], mapping = aes(x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",linewidth = 0.5)+
    stat_summary(bru2[[i]], mapping = aes(x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",linewidth = 0.5)
}


#combine plots
p_br1<-wrap_elements(wrap_plots(br1,ncol=2)+
                       plot_annotation(title = "1 simulation",
                                       tag_levels = "A"))
#ggsave("s_bru_rerun.png",p_br1,width = 25, height = 30, units = 'cm', dpi = 400)

p_br2<-wrap_elements(wrap_plots(br2,ncol=2)+
                       plot_annotation(title = "100 simulstions",
                                       tag_levels = "A")+
                       plot_layout(guides = "collect") &
                       theme(legend.position = "bottom"))

p_br3<-p_br1/p_br2+plot_layout(heights = c(0.9,1.1))
p_br3
#ggsave("s_bru.png",p_br3,width = 25, height = 30, units = 'cm', dpi = 400)

## Fig. Boxplot: % of the total population change for all models ####
a<-dft %>%
  ggplot(aes(x= reorder(model,-Ndiff), y = Ndiff)) +
  geom_boxplot(position = position_dodge(width = .75),size=0.5)+
  #geom_jitter(color="black", size=0.6, alpha=0.4) +
  theme_ipsum() +
  
  theme( 
    legend.position="none",
    plot.title = element_text(size = 13),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title=element_text(size=11),
    legend.text = element_text(size = 11),
    axis.text=element_text(size=11))+
  ggtitle("Gaur population change by infectious disease models") +
  xlab("") +
  ylab("Population change (%)")+
  coord_flip()+ 
  #scale_x_discrete(labels = paste0("M",seq_along(dft_all2$model)))+
  #scale_fill_manual(values= alpha(b25))+ 
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))

a
# ggsave("gaur_ndiff_label.png",a,width = 25, height = 27, units = 'cm', dpi = 600)

## Fig. PCA biplot ####

#biplot/multivariable
#PCA


#import PCA excel file
df_pca <- read_xlsx("./pca_table.xlsx", sheet = "norm") %>% 
  as.data.frame()

df_pca2<-df_pca %>% dplyr::select(c('beta','incubation',"infectious","fatality"))

rownames(df_pca2)<-df_pca$model

mt<-data.matrix(df_pca2)
head(mt)

gpca <- PCA(mt, graph = FALSE,scale.unit = T)

#eigenvalues
get_eigenvalue(gpca)
head(gpca$var$contrib)
head(gpca$var)
cat <- length(unique(df_pca$disease))
pal <- wesanderson::wes_palette("Cavalcanti1", 
                                length(unique(df_pca$disease)), 
                                type = "continuous")[1:length(unique(df_pca$disease))]

#pal2<- (RColorBrewer::brewer.pal(cat,"Spectral"))

pal3 <- wes_palette("Zissou1", 6, type = "continuous")

RColorBrewer::display.brewer.pal(cat,'Spectral')

#correlation plot
library(corrplot)
corrplot(gpca$var$cos2, is.corr=F,col=pal) #component contribution

# Total contribution on PC1 and PC2
# fviz_contrib(gpca, choice = "ind", axes = 1:4)

fviz_pca_var(gpca, col.var = "black")
# Contributions of variables to PC1
fviz_contrib(gpca, choice = "var", axes = 1, top = 10)

# Contributions of variables to PC2
fviz_contrib(gpca, choice = "var", axes = 2, top = 10)
# Disease
fviz_pca_biplot(gpca, label ="var", col.ind="cos2",
                habillage=factor(df_pca$disease)) + 
  theme_minimal()

fviz_pca_biplot(gpca, 
                geom="point", 
                pointsize = 3,
                alpha = 0.7,
                col.var = "black",
                habillage= factor(df_pca$disease),
                legend.title = "Population",
                title = "Biplot: Disease parameters") 

#PCATools
head(df_pca2)
df_t2<-data.table::transpose(df_pca2)

rownames(df_t2) <- colnames(df_pca2)
colnames(df_t2) <- rownames(df_pca2)
head(df_t2)
colnames(df_t2)

metadata<-df_pca
rownames(metadata)
rownames(metadata) <- df_pca$model
head(metadata)

pc<-pca(df_t2,metadata=metadata,
        scale=T,
        center=T)
pc

getComponents(pc)

pairsplot(pc,
          components = getComponents(pc, c(1:4)),
          triangle = F, trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 0.8,
          gridlines.major = FALSE, 
          gridlines.minor = FALSE,
          colby = 'Nchange_p',
          title = 'Pairs plot', 
          plotaxes = T,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))

limits<-c(-200,-100,0,100,200)

bi_pc<-biplot(pc,
       showLoadings = TRUE,
       lab = pc$metadata$model,
       colby = 'Nchange_p',
       hline = 0, vline = 0,
       shape = 'disease',
       sizeLoadingsNames = 5,
       boxedLoadingsNames = F,
       colLoadingsArrows = 'grey30',
       legendPosition = 'right',
       legendLabSize = 11, legendIconSize = 5,
       xlim = c(-10,10),
       ylim = c(-5,5),
       title = "PCA biplot",
       subtitle="Diseases parameters contribute to the % of the population change",
       titleLabSize = 16,
       subtitleLabSize = 15,
       max.overlaps = 15 #ggrepel
)+
  scale_color_gradientn(colours = rev(pal3))+
  guides(color = guide_colorbar(limits = limits))+
  labs(color = "population(%)")
bi_pc
#ggsave("s_biplot_label.png",bi_pc,width = 23, height = 20, units = 'cm', dpi = 600)


#--- done --- :) #
