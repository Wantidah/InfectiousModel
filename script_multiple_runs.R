###### Multiple runs ###########

#### set initial values - change for populations
N = 300 

rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

end.time <- 10 * 365
n_rep <- 

############## MODEL 1 - POPULATION DEMO - NO INFECTION #####
# unit == day (per day)
initials <- c(c = c, sa = sa, a = a )

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

parameters <- c(mu_b = 0.34/365, 
                mu_c = 0.27/365, 
                mu_sa = 0.15/365,
                mu_a = 0.165/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials), 
                tau = 1)


sim_rep_m1<-replicate(n_rep,(model1(pars = parameters, init = initials,
                                             end.time = end.time)))
str(sim_rep_m1[,1]$results)

############## MODEL 2 Anthrax MULTIPLE RUNS  #####
initials <- c(Sc = c, Ic = 0, Ssa = sa, Isa = 0, Sa = a, Ia = 1 )

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

parameters <- c(
  beta_c = 5e-5, beta_sa = 5e-5, beta_a = 5e-5,
  gamma_c = (1/(1/24))/365, gamma_sa = (1/(1/24))/365, gamma_a = (1/(1/24))/365,
  rho_c = 1, rho_sa = 1,  rho_a = 1,
  epsilon = 2e-5,
  mu_b = 0.34/365, 
  mu_c = 0.27/365, mu_sa = 0.15/365, mu_a = 0.165/365,
  delta_c = 1/365, delta_sa = 1/(3*365),
  N = sum(initials),
  tau=1)

sim_rep_m2<-replicate(n_rep,(model2(pars = parameters, init = initials,
                                             end.time = end.time)))
str(sim_rep_m2[,1]$results)
############## MODEL 3 Bovine tuberculosis  MULTIPLE RUNS  #####
initials <- c(Sc = c, Ec = 0, Ic = 0, Ssa = sa, Esa = 0, Isa = 0, Sa = a, Ea = 0, Ia = 1)

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
    
    return(list(pars = pars, init = init2, time = time, results = results))
  }

parameters <- c( 
  beta_c = 0.043/30, beta_sa = 0.043/30, beta_a = 0.043/30,
  phi_c = 0.21/30, phi_sa = 0.21/30, phi_a = 0.21/30,
  gamma_c = 0, gamma_sa = 0, gamma_a = 0,
  rho_c = 0, rho_sa = 0, rho_a = 0.1, 
  epsilon = 2e-5,
  N = sum(initials),
  tau=1,
  mu_b = 0.34/365,  mu_bI = (0.34/365)*(1-0.27), #Ia birth rate reduce by = 27%   
  mu_c = 0.27/365,  mu_sa = 0.15/365, mu_a = 0.165/365,
  delta_c = 1/365, delta_sa = 1/(3*365))

sim_rep_m3<-replicate(n_rep,(model3(pars = parameters, init = initials,
                                             end.time = end.time)))
str(sim_rep_m3[,1]$results)

############## MODEL 4  Hemorrhagic septicemia   MULTIPLE RUNS  #####
initials <- c(Sc = c,  Ic = 0, Rc = 0, Ssa = sa, Isa = 0, Rsa = 0, Sa = a, Ia = 1, Ra = 0)

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

parameters <- c( 
  beta_c = 0.33/365, beta_sa = 0.33/365, beta_a = 0.33/365,
  gamma_c  =1/3/365, gamma_sa =1/3/365, gamma_a =1/3/365,
  rho_c = 0.9, rho_sa = 0.43, rho_a = 0.43,
  omega_c = 1/180/365, omega_sa = 1/180/365, omega_a = 1/180/365,
  epsilon = 2e-5,
  mu_b = 0.34/365, mu_c = 0.27/365, mu_sa = 0.15/365, mu_a = 0.165/365,
  delta_c = 1/365, delta_sa = 1/(3*365),
  N = sum(initials),
  tau=1)

sim_rep_m4<-replicate(n_rep,(model4(pars = parameters, init = initials,
                                    end.time = end.time)))
str(sim_rep_m4[,1]$results)

############## MODEL 5  Lumpy Skin Disease   MULTIPLE RUNS  #####
initials <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = a, Ea = 0, Ia = 1, Ra = 0)

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


#SEIR parameter
parameters <- c( 
  beta_c = 0.043/30, beta_sa = 0.043/30, beta_a = 0.043/30, 
  phi_c = 1/7/365,phi_sa = 1/7/365,phi_a = 1/7/365,
  gamma_c = 0, gamma_sa = 0,gamma_a = 0,
  rho_c = 0.03,rho_sa = 0.03, rho_a = 0.03, 
  omega_c = (1/365)/365, omega_sa =  (1/365)/365, omega_a = (1/365)/365,
  epsilon = 2e-5,
  mu_b = 0.34/365,  mu_bI = (0.34/365)*(1-0.27), #Ia birth rate reduce by = 27%   
  mu_c = 0.27/365,  mu_sa = 0.15/365, mu_a = 0.165/365,
  delta_c = 1/365, delta_sa = 1/(3*365),
  N = sum(initials),
  tau=1)

sim_rep_m5<-replicate(n_rep,(model5(pars = parameters, init = initials,
                                    end.time = end.time)))
str(sim_rep_m5[,1]$results)
############## MODEL 6 Foot and mouth disease  MULTIPLE RUNS  #####
initials <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = a, Ea = 0, Ia = 1, Ra = 0)

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

parameters <- c( 
  beta_c = 0.52/365, beta_sa = 0.52/365, beta_a = 0.52/365,
  phi_c = 1/8/365, phi_sa = 1/6/365, phi_a = 1/6/365,
  gamma_c = 1/5/365, gamma_sa = 1/5/365, gamma_a = 1/5/365,
  rho_c = 0.1, rho_sa = 0.05, rho_a = 0.03, 
  alpha = 0.5,
  omega_c = (1/120)/365, omega_sa =  (1/120)/365, omega_a = (1/565)/365, omega_m = (1/144)/365,
  epsilon = 2e-5,
  mu_b = 0.34/365, mu_bI = (0.34/365)*(0.9), #Ia birth rate reduce by = 10%  (assume)
  mu_c = 0.27/365, mu_sa = 0.15/365, mu_a = 0.165/365,
  delta_c = 1/365, delta_sa = 1/(3*365),
  N = sum(initials),
  tau=1)

sim_rep_m6<-replicate(n_rep,(model6(pars = parameters, init = initials,
                                    end.time = end.time)))
str(sim_rep_m6[,1]$results)
############## MODEL 7 Bovine brucellosis MULTIPLE RUNS  #####
initials <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = a, Ea = 0, Ia = 1, Ra = 0)

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
                        Sc, Ec, Ic, Rc, Ssa, Esa, Isa, Rsa, Sa, Ea, Ia, Ra, M, Sm)%>% 
      dplyr::mutate(N = rowSums(across(-c(time), na.rm=TRUE)))%>% 
      dplyr::mutate(S = rowSums(across(c(Sa,Ssa,Sc,Sm)), na.rm=TRUE))%>% 
      dplyr::mutate(E = rowSums(across(c(Ea,Esa,Ec)), na.rm=TRUE))%>% 
      dplyr::mutate(I = rowSums(across(c(Ia,Isa,Ic)), na.rm=TRUE))%>% 
      dplyr::mutate(R = rowSums(across(c(Ra,Rsa,Rc)), na.rm=TRUE))
    
    return (list(pars = pars, init = init2, time = time, results = results))
    
  }

parameters <- c( 
  beta_c = 2/365, beta_sa = 2/365,beta_a = 2/365,
  phi_c = 1/14/365, phi_sa = 1/14/365, phi_a = 1/14/365,
  gamma_c = 1/2*365, gamma_sa = 1/2*365, gamma_a = 1/2*365,
  rho_c = 0.1, rho_sa = 0.05, rho_a = 0.03, 
  alpha = 0.9,
  omega_c = (1/180)/365, omega_sa =  (1/180)/365, omega_a = (1/180)/365, omega_m = (1/180)/365,
  epsilon = 2e-5,
  mu_b = 0.34/365, mu_bI = (0.34/365)*(0.5), #Ia birth rate reduce by = 10%  (assume)
  mu_c = 0.27/365, mu_sa = 0.15/365, mu_a = 0.165/365,
  delta_c = 1/365, delta_sa = 1/(3*365),
  N = sum(initials),
  tau=1
)

sim_rep_m7<-replicate(n_rep,(model7(pars = parameters, init = initials,
                                    end.time = end.time)))
str(sim_rep_m7[,1]$results)
############# 

###### FUNCTION single_pop_sim_prep ##############
#change from results$I to results$Ia
#However, this Ia need to be changed to total I and N (total population) later
single_pop_sim_prep <- function(x, n_rep, end.time){ 
  # x = simulation of model, e.g. sim_run_m1
  mat = matrix(NA, nrow=n_rep, ncol = end.time+1)
  for (i in 1:n_rep){
    mat[i,]<-x[,i]$results$I #***need to add I, N before running***#
  }
  colnames(mat) = paste("time", seq(from=1,to=end.time+1,by=1), sep="")
  rownames(mat) = paste("run", seq(n_rep), sep="")
  dat = as.data.frame(mat)
  dat$run = rownames(dat)
  mdat = melt(dat, id.vars="run")
  mdat$time = as.numeric(gsub("time", "", mdat$variable))
  mdat
}

## SINGLE POPULATION PRINT PREPS #######
########### list simulations multiple runs

res_mx<-list(sim_rep_m1,
             sim_rep_m2,
             sim_rep_m3,
             sim_rep_m4,
             sim_rep_m5,
             sim_rep_m6)

for (i in 1:length(res_mx)){
  assign(paste0("df_mx", i), single_pop_sim_prep(x=res_mx[[i]],n_rep = n_rep, end.time = end.time))
}

res_mx_p<-rbind(df_mx1,
                df_mx2,
                df_mx3,
                df_mx4,
                df_mx5,
                df_mx6)

res_mx_p$model<-c(rep('1',dim(df_mx1)[1]),
                  rep('2',dim(df_mx2)[1]),
                  rep('3',dim(df_mx3)[1]),
                  rep('4',dim(df_mx4)[1]),
                  rep('5',dim(df_mx5)[1]),
                  rep('6',dim(df_mx6)[1]))

## SINGLE POPULATION PLOTS #######
for (i in unique(res_mx_p$model)){
  subdata <- subset(res_mx_p, model == i)
  pdf(paste("plot_ts_mx", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot(subdata, aes(x=time, y=value, group=run)) +
          theme_bw() +
          theme(panel.grid=element_blank()) +
          geom_line(size=0.2, alpha=0.15)+
          ylab('Numbers') + xlab('time')+
          stat_summary(aes(group = 1), fun.y=mean, geom="line", colour="black",size = 1.1))
  dev.off()
}

## plot extinction times

df.agg <- aggregate(time ~ run + value + model, res_mx_p, min)
df.ag<-(df.agg[df.agg$value==0,c('model','time')])

neworder <- c("1","2","3","4","5","6")
library(plyr)  ## or dplyr (transform -> mutate)
df.ag <- arrange(transform(df.ag,
                           model=factor(model,levels=neworder)),model)
labs <- c('1' = "SI",
          '2' = "SEI",
          '3' = "SIRS",
          '4' = "SEIR",
          '5' = "SEIRM[FMD]",
          '6' = "SEIRM[Bru]")
p<-ggplot(df.ag, aes(x=time))+
  geom_histogram(color="black", fill="grey")+
  facet_wrap(model~., ncol = 6,labeller = labeller(model = labs))+
  scale_x_continuous(breaks = c(0, 800, 1600), labels = c("0", "800", "1600"))
pdf("extinctions.pdf", width = 8, height = 3)
p
dev.off()

p_yr<-ggplot(df.ag, aes(x=time))+
  geom_histogram(color="black", fill="grey")+
  facet_wrap(model~., ncol = 6,labeller = labeller(model = labs))+
  scale_x_continuous(limits = c(0,1000),breaks = c(0, 400, 800), labels = c("0", "400", "800")) +
  ylim(0,60)
pdf("extinctions_1yr.pdf", width = 6, height = 3)
p_yr
dev.off()
