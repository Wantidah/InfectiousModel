############## 4) SIRS model  (Hemorrhagic septicemia)  #####
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
    return(list(pars = pars, 
                init = init2, 
                time = time, 
                results = data.frame(time, 
                                     Sc,  Ic, Rc, Ssa, Isa, Rsa, Sa, Ia, Ra)))
  }


#estimate the age structure proportion
N=300
#calf:subadult:adult ratio
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

initials <- c(Sc = c,  Ic = 0, Rc = 0, Ssa = sa, Isa = 0, Rsa = 0, Sa = a, Ia = 1, Ra = 0)

end.time <- 100*365 #predict for ... years

#SEI parameter
#same transmission rate (beta), case fatality rate (rho)
#gaur
parameters <- c( 
  beta_c = 0.33/365,
  beta_sa = 0.33/365,
  beta_a = 0.33/365,
  gamma_c  =1/3/365,
  gamma_sa =1/3/365,
  gamma_a =1/3/365,
  rho_c = 0.9,
  rho_sa = 0.43,
  rho_a = 0.43,
  omega_c = 1/180/365,
  omega_sa = 1/180/365, 
  omega_a = 1/180/365,
  epsilon = 2e-5,
  mu_b = 0.34/365, 
  mu_c = 0.27/365, 
  mu_sa = 0.15/365,
  mu_a = 0.165/365,
  delta_c = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials),
  tau=1
)


####### plot SIR ######
res_sir_gaur <- model4(pars = parameters, init = initials,
                       end.time = end.time)

min(subset(res_sir_gaur$results,N==0)$time)

PlotMods(res_sir_gaur)
str(res_sir_gaur)

#sum of populations
res_sir_gaur$total<-rowSums(res_sir_gaur$results[,2:9])

plot(rowSums(res_sir_gaur$results[,2:9]), main = "HS: gaur total population", 
     xlab="time",ylab="animal")

res_sir_gaur_df<-data.frame(res_sir_gaur$total)

View(res_sir_gaur_df)

#combine plot
#get the total population compared between non-infectious and infectious
non<-res_gaur_df
inf<-res_sir_gaur_df

tot_df <-cbind(non, inf)
tot_df$time <-seq.int(nrow(tot_df))

str(tot_df)
tail(tot_df)

end.time 

if (end.time < 36500) {
  ggplot() + 
    geom_line(data = tot_df,aes(x = time ,y = res_sir_gaur.total, color = 'res_si_gaur.total')) + 
    geom_line(data = tot_df,aes(x = time, y = res_g.total, color = 'res_g.total' ))+
    labs(x="days", y="total population (N)",
         title='Gaur total population in 20 years') +
    scale_color_manual(name = "N",
                       labels = c('non-infection','HS'),
                       values = c('#009988','#cc6677'))+ #0c7bdc #0077bb
    theme( plot.title = element_text(size = 18),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 15),
           legend.title=element_text(size=11),
           legend.text = element_text(size = 11),
           axis.text=element_text(size=13)) 
  
} else {
  ggplot() + 
    geom_line(data = tot_df,aes(x = time ,y = res_sir_gaur.total, color = 'res_si_gaur.total')) + 
    geom_line(data = tot_df,aes(x = time, y = res_g.total, color = 'res_g.total' ))+
    labs(x="days", y="total population (N)",
         title='Gaur total population in 100 years') +
    scale_color_manual(name = "N",
                       labels = c('non-infection','HS'),
                       values = c('#009988','#cc6677'))+ #0c7bdc #0077bb
    theme( plot.title = element_text(size = 18),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 15),
           legend.title=element_text(size=11),
           legend.text = element_text(size = 11),
           axis.text=element_text(size=13))
}
