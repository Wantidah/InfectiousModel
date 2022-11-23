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
    return(list(pars = pars, init = init2, time = time, results = data.frame(time, 
                                                                             Sc, Ic, Ssa, Isa, Sa,Ia)))
  }

#estimate the age structure proportion
N=300
#calf:subadult:adult ratio
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

#SI parameters
#same transmission rate (beta), case fatality rate (rho)
#gaur
parameters <- c(
  beta_c = 5e-5,
  gamma_c = (1/(1/24))/365,
  rho_c = 1,
  beta_sa = 5e-5,
  gamma_sa = (1/(1/24))/365,
  rho_sa = 1,
  beta_a = 5e-5,
  gamma_a = (1/(1/24))/365,
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

#single run
initials <- c(Sc = c, Ic = 0, Ssa = sa, Isa = 0, Sa = a, Ia = 1 )
end.time <- 100*365 #predict for ... years


## MULTIPLE RUNS ################################

#### set initial values - change for meta populaions
initials <- c(Sc = c, Ic = 0, Ssa = sa, Isa = 0, Sa = a, Ia = 1 )
end.time <- 20*365 #predict for ... years
n_rep <- 10
############## MODEL 2) SI Anthrax MULTIPLE RUNS  #####


sim_run_anthrax_rep<-replicate(n_rep,(model2(pars = parameters, init = initials,
                                             end.time = end.time)))


####### plot SI ######
res_si_gaur <- model2(pars = parameters, init = initials,
                      end.time = end.time)
PlotMods(res_si_gaur)


#sum of populations
res_si_gaur$total<-rowSums(res_si_gaur$results[,2:6])

plot(rowSums(res_si_gaur$results[,2:6]), main = "Anthrax: gaur total population", 
     xlab="time",ylab="animal")


res_si_gaur_df<-data.frame(res_si_gaur$total)
View(res_si_gaur_df)
#buffalo
N = 70
rat = 1+6+5
N/rat
c = (N/rat)*1
sa = (N/rat)*6
a = (N/rat)*5

N==a+sa+c

initials <- c(Sc = c, Ic = 0, Ssa = sa, Isa = 0, Sa = a, Ia = 1 )

end.time <- 100*365 

parameters <- c(
  beta_c = 0.0001,
  gamma_c = (1/(1/24))/365,
  rho_c = 1,
  beta_sa = 0.0001,
  gamma_sa = (1/(1/24))/365,
  rho_sa = 1,
  beta_a = 0.0001,
  gamma_a = (1/(1/24))/365,
  rho_a = 1,
  epsilon = 2e-5,
  N = sum(initials),
  tau=1,
  
  mu_b = 0.37/365, 
  mu_c = 0.27/365, 
  mu_sa = 0.15/365,
  mu_a = 0.15/365,
  delta_c = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials), 
  tau = 1)

#plot
res_si <- model2(pars = parameters, init = initials,
                 end.time = end.time)
PlotMods(res_si)

View(res_si$results)

min(subset(res_si$results,Isa==0)$time)

#sum of populations
res_si$total<-rowSums(res_si$results[,2:6])
plot(rowSums(res_si$results[,2:6]), main = "Anthrax: buffalo total population", 
     xlab="time",ylab="animal")


#combine plot
#get the total population compared between non-infectious and infectious
non<-res_gaur_df
inf<-res_si_gaur_df

tot_df <-cbind(non, inf)
tot_df$time <-seq.int(nrow(tot_df))

str(tot_df)

ggplot() + 
  geom_line(data = tot_df,aes(x = time ,y = res_si_gaur.total, color = 'res_si_gaur.total')) + 
  geom_line(data = tot_df,aes(x = time, y = res_g.total, color = 'res_g.total' ))+
  labs(x="days", y="total population (N)",
       title='Gaur total population in 100 years') +
  scale_color_manual(name = "N",
                     labels = c('non-infection','anthrax'),
                     values = c('#009988','#cc6677'))+ #0c7bdc #0077bb
  theme( plot.title = element_text(size = 18),
         axis.title.x = element_text(size = 15),
         axis.title.y = element_text(size = 15),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=13))

?scale_colour_manual
