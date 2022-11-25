############## NOTES ##################
# Reference R Code:  https://doi.org/10.1098/rsif.2021.0638  (Transmission models indicate Ebola virus persistence in non-human primate populations is unlikely)
# https://github.com/dtsh2/ebola_model

### SEIR of infectious disease model####

############# LOAD PACKAGES ##########
#package<-c("reshape","EpiDynamics", "plyr","emdbook", "reshape2","stringr","ggplot2")
#install.packages(package)

rm(list = ls())

library(reshape)
library(EpiDynamics)
library(tidyverse)
library(plyr)     
library(reshape2) 
library(stringr)
library(emdbook)  
library(ggplot2); theme_set(theme_bw())
library(SciViews)

## MODELS ################################
############## MODEL 1 SIR WITH MORTALITY IN I ######
model1 =
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 9)
        change <- matrix(0, nrow = 9, ncol = 4)
        N <- S + E+ I + R
        tau <- 1
        
        rate[1] <- beta * S * I/N    
        change[1, ] <- c(-1, 1, 0, 0) 
        
        rate[2] <- phi * E
        change[2, ] <- c( 0, -1, 1, 0)
        
        rate[3] <- mu_b * N
        change[3, ] <- c(1, 0, 0, 0)
        
        rate[4] <- mu_d* I
        change[4, ] <- c(0, 0, -1, 0)
        
        rate[5] <- (1-rho) * gamma * I
        change[5, ] <- c(0, 0, -1, 1)
        
        rate[6] <- mu_d * S
        change[6, ] <- c(-1, 0, 0, 0)
        
        rate[7] <- mu_d * E
        change[7, ] <- c(0, -1, 0, 0)
        
        rate[8] <- mu_d * R
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

############## --> TEST MODEL 1 WITH BASIC PLOT ######
#Initial population

initials <- c(S = 1000, E = 0, I = 1, R = 0)

end.time <- 20*365

parameters <- c(beta = 0.33,
                phi = 0.3,
                rho = 0.05,
                gamma = 0.3,
                mu_b = 0.5/365, #birth
                mu_d = 0.1/365, #death #1/arv lifespan (days)
                N = sum(initials), tau =1)

res <- model1(pars = parameters, init = initials,
              end.time = end.time)
PlotMods(res)

min(subset(res$results,I==0)$time)

#sum of populations
plot(rowSums(res$results[,2:5]))


#2 ages classes
# population death rate
r=0.1
#time (years)
t=1  

#probability of death
p = 1-(exp(-r*t))

mu_d=p/365

mu_d
############## MODEL (#6) SEIR WITH MORTALITY in Adult & Juvenile ####
#No Epsilon
#No R->S

model6=
  function (pars, init, end.time)  {
    init2 <- init
    Equations <- function(pars, init, end.time) {
      with(as.list(c(pars, init)), {
        rate <- rep(0, 21)
        change <- matrix(0, nrow = 21, ncol = 8)
        
        Nj <- Sj + Ej+ Ij + Rj 
        Na <- Sa + Ea + Ia + Ra
        
        tau <- 1
        
        rate[1] <- beta_j * Sj * (Ij+Ia)/(Na+Nj)        
        change[1, ] <- c(-1, +1, 0, 0, 0, 0, 0, 0) 
        
        rate[2] <- phi_j * Ej
        change[2, ] <- c(0, -1, +1, 0, 0, 0, 0, 0)
        
        rate[3] <- mu_dj * Ij
        change[3, ] <- c(0, 0, -1, 0, 0, 0, 0, 0)
        
        rate[4] <- (1-rho_j) * gamma_j * Ij
        change[4, ] <- c(0, 0, -1, +1, 0, 0, 0, 0)
        
        rate[5] <- mu_dj * Sj
        change[5, ] <- c(-1, 0, 0, 0, 0, 0, 0, 0)
        
        rate[6] <- mu_dj * Ej
        change[6, ] <- c(0, -1, 0, 0, 0, 0, 0, 0)
        
        rate[7] <- mu_dj * Rj
        change[7, ] <- c(0, 0, 0, -1, 0, 0, 0, 0)
        
        rate[8] <- rho_j * gamma_j * Ij
        change[8, ] <- c(0, 0, -1, 0, 0, 0, 0, 0)
        
        rate[9] <- delta * Sj
        change[9, ] <- c(0, 0, 0, 0, +1, 0, 0, 0)
        
        rate[10] <- delta * Ej
        change[10, ] <- c(0, 0, 0, 0, 0, +1, 0, 0)
        
        rate[11] <- delta * Ij
        change[11, ] <- c(0, 0, 0, 0, 0, 0, +1, 0)
        
        rate[12] <- delta * Rj
        change[12, ] <- c(0, 0, 0, 0, 0, 0, 0, +1)
        
        rate[13] <- mu_b * Na
        change[13, ] <- c(+1, 0, 0, 0, 0, 0, 0, 0)
        
        #Adult
        rate[14] <- beta_a * Sa * (Ij+Ia)/(Na+Nj)        
        change[14, ] <- c(0, 0, 0, 0, -1, +1, 0, 0) 
        
        rate[15] <- phi_a * Ea
        change[15, ] <- c(0, 0, 0, 0, 0, -1, +1, 0)
        
        rate[16] <- mu_da * Ia
        change[16, ] <- c(0, 0, 0, 0, 0, 0, -1, 0)
        
        rate[17] <- (1-rho_a) * gamma_a * Ia
        change[17, ] <- c(0, 0, 0, 0, 0, 0, -1, +1)
        
        rate[18] <- mu_da * Sa
        change[18, ] <- c(0, 0, 0, 0, -1, 0, 0, 0)
        
        rate[19] <- mu_da * Ea
        change[19, ] <- c(0, 0, 0, 0, 0, -1, 0, 0)
        
        rate[20] <- mu_da * Ra
        change[20, ] <- c(0, 0, 0, 0, 0, 0, 0, -1)
        
        rate[21] <- rho_a * gamma_a * Ia
        change[21, ] <- c(0, 0, 0, 0, 0, 0, -1, 0)
        
        init <- c(Sj = Sj, Ej =Ej, Ij = Ij, Rj = Rj, Sa = Sa, Ea = Ea, Ia = Ia, Ra = Ra)
        for (i in 1:21) {
          num <- rpois(1, rate[i] * tau)
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
          init <- init + change[i, ] * num.min
        }
        return(init)
      })
    }
    
    Sj <- Ej <- Ij <- Rj <- Sa <- Ea <- Ia <- Ra <- double()
    t <- 0
    time <- seq(0, end.time, by = pars["tau"])
    for (t in time) {
      tmp <- Equations(pars, init, end.time)
      Sj <- c(Sj, init["Sj"])
      Ej <- c(Ej, init["Ej"])
      Ij <- c(Ij, init["Ij"])
      Rj <- c(Rj, init["Rj"])
      
      Sa <- c(Sa, init["Sa"])
      Ea <- c(Ea, init["Ea"])
      Ia <- c(Ia, init["Ia"])
      Ra <- c(Ra, init["Ra"])
      init <- tmp
    }
    return(list(pars = pars, init = init2, time = time, results = data.frame(time, 
                                                                             Sj, Ej, Ij, Rj,
                                                                             Sa, Ea, Ia, Ra)))
  }

#parameters

beta_a = 0.3
beta_j = 0.3
phi_a = 0.3
phi_j = 0.3/2
rho_a = 0.05
rho_j = 0.5
gamma_a = 0.3
gamma_j = 0.3/2

####Birth rate & Death rate ####
mu_b = ln(69/47)/(((15*365))+1/(22*365))*2
mu_da = 1/(22*365) #1/maximum lifespan (days) ; lifespan = 22 years
mu_dj = 1/(22*365)*2 

mu_b = 1/(365*2) # average lifespan = 8 years
mu_da = 1/(10*365) #
mu_dj = 1/(10*365)*2

mu_b = ln(69/47)/(((15*365))+1/(22*365))*2
mu_da = 1/(10*365) #
mu_dj = 1/(10*365)*2

###############################

delta = 1/(365*3)

Sj = 500
Ej= 0
Ij = 0
Rj = 0

Sa = 500
Ea = 0
Ia = 1
Ra = 0
############## --> TEST MODEL 6 WITH BASIC PLOT ######
initials <- c(Sj = 500, Ej= 0, Ij = 0, Rj = 0, Sa = 500 , Ea = 0, Ia = 1, Ra = 0)
end.time <- 20*365

parameters <- c(beta_a = beta_a, 
                beta_j = beta_j,
                phi_a = phi_a,
                phi_j = phi_j,
                rho_j = rho_j,
                rho_a = rho_a,
                gamma_a = gamma_a,
                gamma_j = gamma_j,
                mu_b = mu_b, 
                mu_da = mu_da, 
                mu_dj=mu_dj,
                Nj = sum(Sj,Ej,Ij,Rj),
                Na = sum(Sa,Ea,Ia,Ra), 
                tau = 1,
                #epsilon = 2e-5, 
                delta = delta)

res <- model6(pars = parameters, init = initials,
              end.time = end.time)
PlotMods(res)
df_meta_sir<-res$results

min(subset(res$results,Ia==0)$time)

min(subset(res$results,Ij==0)$time)

max(subset(res$results,Ia!=0)$time)

max(subset(res$results,Ij!=0)$time)

#sum of populations
plot(rowSums(res$results[,2:9]))

head(res$result)

