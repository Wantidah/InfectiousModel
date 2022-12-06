############## LOAD PACKAGES ##########
library(EpiDynamics)
library(dplyr)   
library(tidyverse)
library(reshape2) 
library(stringr)
library(ggplot2); theme_set(theme_bw())

############## 2) SI model (Anthrax)  #####
model2=
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
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
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

# gaur population
N=300

#calf:subadult:adult ratio
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

initials <- c(Sc = c, Ic = 0, Ssa = sa, Isa = 0, Sa = (a-1), Ia = 1 )

end.time <- 100*365 #predict for ... years

#SI parameters
parameters <- c(
  beta_c = 5e-5,
  gamma_c = 1,
  rho_c = 1,
  beta_sa = 5e-5,
  gamma_sa = 1,
  rho_sa = 1,
  beta_a = 5e-5,
  gamma_a = 1,
  rho_a = 1,
  epsilon = 2e-5,
  N = sum(initials),
  tau=1,
  mu_b = 0.34/365, 
  mu_c = 0.27/365, 
  mu_sa = 0.15/365,
  mu_a = 0.165/365,
  delta_c = 1/365,
  delta_sa = 1/(3*365)
)

# TEST
# single run
res_model2 <-model2(pars = parameters, init = initials,
       end.time = end.time)

str(res_model2)

#minimum I,N extinction
min(subset(res_model2$results,I==0)$time)
min(subset(res_model2$results,N==0)$time)

#plot epi
#PlotMods(res_model2)

#convert to data.frame, change days -> years
res_model2<-res_model2$results %>%
  mutate(time_y = time/365) %>% #convert day to year for plotting
  as.data.frame()

# plot SI Anthrax ######
p<-ggplot() + 
  geom_line(data = res_model2,aes(x = time_y ,y = S, color = 'S')) + 
  geom_line(data = res_model2,aes(x = time_y, y = I, color = 'I' ))+
  geom_line(data = res_model2, aes(x = time_y, y = N,color = 'total'))+
  
  labs(x="years", y= "population",
       title= 'Gaur population with anthrax infection, 1 simulation, 100 years') +
  
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  
  scale_color_manual( name = "population",
                      labels = c('S','I','total' ),
                      values = c('S'='seagreen4',
                                 'I'='firebrick',
                                 "total"='#153030'))+ 
  theme_bw() +
  theme( plot.title = element_text(size = 18),
         axis.title.x = element_text(size = 15),
         axis.title.y = element_text(size = 15),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=13))+
  guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))

ggsave("gaur_anthrax_1sim_100y_all.png",p, width = 25, height = 15, units = 'cm', dpi = 600)
