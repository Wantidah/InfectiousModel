############## LOAD PACKAGES ##########
library(EpiDynamics)
library(dplyr)   
library(tidyverse)
library(reshape2) 
library(stringr)
library(ggplot2); theme_set(theme_bw())

############## 3) SEI model (Bovine tuberculosis) ##### 
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
          num.min <- min(num, init[which(change[i, ] < 
                                           0)])
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
    
     return(list(pars = pars, 
                init = init2, 
                time = time, 
                results = results))
}

#gaur population 
N=300
#calf:subadult:adult ratio
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

initials <- c(Sc = c, Ec = 0, Ic = 0, Ssa = sa, Esa = 0, Isa = 0, Sa = a, Ea = 0, Ia = 1)

end.time <- 100*365 #predict for ... years

#SEI parameter
parameters <- c( 
  beta_c = 0.043/30,
  beta_sa = 0.043/30,
  beta_a = 0.043/30,
  phi_c = 0.21/30,
  phi_sa = 0.21/30,
  phi_a = 0.21/30,
  gamma_c = 0,
  gamma_sa = 0,
  gamma_a = 0,
  rho_c = 0,
  rho_sa = 0,
  rho_a = 0.1, 
  epsilon = 2e-5,
  N = sum(initials),
  tau=1,
  
  mu_b = 0.34/365, 
  mu_bI = (0.34/365)*(1-0.27), #Ia birth rate reduce by = 27%   
  mu_c = 0.27/365, 
  mu_sa = 0.15/365,
  mu_a = 0.165/365,
  delta_c = 1/365,
  delta_sa = 1/(3*365)
)

# TEST
# single run
res_model3 <- model3(pars = parameters, init = initials,
                     end.time = end.time)

str(res_model3)

#minimum I,N extinction
min(subset(res_model3$results,I==0)$time)
min(subset(res_model3$results,N==0)$time)

#plot epi
#PlotMods(res_model3)

#convert to data.frame, change days -> years
res_model3<-res_model3$results %>%
  mutate(time_y = time/365) %>% #convert day to year for plotting
  as.data.frame()

# plot SIR HS ######
p<-ggplot() + 
  geom_line(data = res_model3,aes(x = time_y ,y = S, color = 'S')) + 
  geom_line(data = res_model7,aes(x = time_y, y = E,  color = 'E' ))+
  geom_line(data = res_model3,aes(x = time_y, y = I, color = 'I' ))+
  geom_line(data = res_model3, aes(x = time_y, y = N,color = 'total'))+
  
  labs(x="years", y= "population",
       title= 'Gaur population with bTB infection, 1 simulation, 100 years') +
  
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  
  scale_color_manual( name = "population",
                      labels = c('S','E','I','total' ),
                      values = c('S'='seagreen4',
                                 'E'='darkorange2',
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

ggsave("gaur_bTB_1sim_100y_all.png",p, width = 25, height = 15, units = 'cm', dpi = 600)


