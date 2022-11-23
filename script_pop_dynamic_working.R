# Infectious disease modelling code
# - 18 Nov 2022-
# # Reference R Code:  https://doi.org/10.1098/rsif.2021.0638  (Transmission models indicate Ebola virus persistence in non-human primate populations is unlikely)
# https://github.com/dtsh2/ebola_model

############## LOAD PACKAGES ##########

library(reshape)
library(EpiDynamics)
library(plyr)     
library(reshape2) 
library(stringr)
library(emdbook)  
library(ggplot2); theme_set(theme_bw())

## MODELS ################################

# Population dynamic model, no infection, 3 age classes
# c = calf
# sa = Subadult
# a = adult

## set initial values ## 
#initials <- c(a = 1000, sa = 100, c = 100)

############## MODEL 1 population dynamic model, no infection, 3 age-classes ######
# unit == day (per day)
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
    return(list(pars = pars, init = init2, time = time, results = data.frame(time, 
                                                                             c, sa, a)))
  }


############## --> TEST MODEL 1 Population dynamic, no infection ######

#gaur population (Khao Phaeng Ma Non-Hunting Area, Thailand)
N = 300 

#estimate the age structure proportion
#calf:subadult:adult ratio
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

initials <- c(c = c, sa = sa, a = a )

end.time <- 100*365 #predict for ... years

parameters <- c(mu_b = 0.34/365, 
                mu_c = 0.27/365, 
                mu_sa = 0.15/365,
                mu_a = 0.165/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials), 
                tau = 1)


res_g <- model1(pars = parameters, init = initials,
                end.time = end.time)

PlotMods(res_g)

res_g$total<-rowSums(res_g$results[,2:4])

plot(rowSums(res_g$results[,2:4]), 
     main = "gaur total population",
     xlab="time",ylab="animal")

res_gaur_df<-data.frame(res_g$total)
str(res_gaur_df)


# banteng population
N = 290
#calf:subadult:adult ratio
rat = 1.1+1.3+1
c = N/rat*1.1
sa = N/rat*1.3
a = N/rat*1
c+sa+a

initials <- c(c = c, sa = sa, a = a )
end.time <- 20*365 #predict for ... years

parameters <- c(mu_b = 0.4/365, 
                mu_c = 0.26/365, 
                mu_sa = 0.26/365,
                mu_a = 0.15/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials), 
                tau = 1)


res_bt <- model1(pars = parameters, init = initials,
                 end.time = end.time)
PlotMods(res_bt)

#sum of populations
plot(rowSums(res_bt$results[,2:4]), 
     main = "banteng total population",
     xlab="time",ylab="animal")

#buffalo population (HKK, Thailand)
N = 70
rat = 1+6+5
N/rat
c = (N/rat)*1
sa = (N/rat)*6
a = (N/rat)*5


initials <- c(c = c, sa = sa, a = a )
end.time <- 100*365 #predict for ... years

parameters <- c(mu_b = 0.37/365, 
                mu_c = 0.27/365, 
                mu_sa = 0.15/365,
                mu_a = 0.15/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials), 
                tau = 1)

res_buf<- model1(pars = parameters, init = initials,
                 end.time = end.time)
PlotMods(res_buf)

#sum of populations
res_buf$total<-rowSums(res_buf$results[,2:4])
plot(res_buf$total, main = "buffalo total population", 
     xlab="time",ylab="animal")

View(res_buf$results)

res_buf_df<-data.frame(res_buf$total)
View(res_buf_df)
#serow population
#assume, mortality as an average for mammal
N = 120 

rat = 1+1+1
N/rat
c = round((N/rat)*1,0)
sa = round((N/rat)*1,0)
a = round((N/rat)*1,0)
N==a+sa+c

parameters <- c(mu_b = 0.7/365, 
                mu_c = 0.5/365, 
                mu_sa = 0.15/365,
                mu_a = 0.28/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials), 
                tau = 1)

res_se<- model1(pars = parameters, init = initials,
                end.time = end.time)
PlotMods(res_se)

#sum of populations
res_se$total<-rowSums(res_se$results[,2:4])
plot(res_se$total, main = "serow total population", 
     xlab="time",ylab="animal")

#goral population
#assume
N = 292

rat = 1+1+1
N/rat
c = round((N/rat)*1,0)
sa = round((N/rat)*1,0)
a = round((N/rat)*1,0)
N==a+sa+c

parameters <- c(mu_b = 0.5/365, 
                mu_c = 0.45/365, 
                mu_sa = 0.27/365,
                mu_a = 0.17/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials), 
                tau = 1)


res_gor<- model1(pars = parameters, init = initials,
                 end.time = end.time)
PlotMods(res_gor)

#sum of populations
res_gor$total<-rowSums(res_gor$results[,2:4])
plot(res_gor$total, main = "goral total population", 
     xlab="time",ylab="animal")

