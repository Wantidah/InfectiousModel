#TRY
#model 7 the population extinction from brucellosis
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

#> MODEL 7 SEIRMS/E - Brucellosis  #####
# gaur population 
N = 300 

#estimate the age structure proportion
#calf:subadult:adult ratio
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 
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

## SET UP FUNCTIONS FOR MX METRICS ################################

############## FUNCTION  1 COUNT EXTINCTIONS ########

my_min_ext_I<-function(x=get_time$results, y=get_time$results$I,...){
  res_no<-vector()
  if(!is.na(min(subset(x,y==0)$time))){
    res_no = min(subset(x,y==0)$time)
  } else{
    res_no = end.time
  }
  res_no
}

my_min_ext_N<-function(x=get_time$results, y=get_time$results$N,...){
  res_no<-vector()
  if(!is.na(min(subset(x,y==0)$time))){
    res_no = min(subset(x,y==0)$time)
  } else{
    res_no = end.time
  }
  res_no
}

############## FUNCTION  2 COUNT EXTINCTIONS - NO INF #########
#count I == 0
my_imp_ext_I<-function(x=time,y=I,...){
  res_no<-vector()
  for (i in 2:length(x)){
    if(is.na(y[i]))
    {res_no[i-1] = NA} 
    else
      if(y[i-1]>0 & y[i]==0 & !is.na(y[i])){ 
        res_no[i-1] = 1
      } else{
        res_no[i-1] = 0
      }
  }
  sum(res_no, na.rm = T)
}

#count N == 0
my_imp_ext_N<-function(x=time,y=N,...){
  res_no<-vector()
  for (i in 2:length(x)){
    if(is.na(y[i]))
    {res_no[i-1] = NA} 
    else
      if(y[i-1]>0 & y[i]==0 & !is.na(y[i])){ 
        res_no[i-1] = 1
      } else{
        res_no[i-1] = 0
      }
  }
  sum(res_no, na.rm = T)
}

############## FUNCTION  3 COUNT PERSISTENCE TIMES ########

my_imp_ext_na_I<-function(x=time,y=I,...){
  res_no<-vector()
  for (i in 2:length(x)){
    if(is.na(y[i]))
    {res_no[i-1] = NA}
    else
      if(y[i-1]>0 & y[i]==0 & !is.na(y[i])){
        res_no[i-1] = 1
      } else{
        res_no[i-1] = 0
      }
  }
  length(res_no[!is.na(res_no)])
}

my_imp_ext_na_N<-function(x=time,y=N,...){
  res_no<-vector()
  for (i in 2:length(x)){
    if(is.na(y[i]))
    {res_no[i-1] = NA}
    else
      if(y[i-1]>0 & y[i]==0 & !is.na(y[i])){
        res_no[i-1] = 1
      } else{
        res_no[i-1] = 0
      }
  }
  length(res_no[!is.na(res_no)])
}


## SET UP FUNCTIONS FOR MODEL RUNS LOOP THROUGH PARS ################################

# brucellosis
beta <- seq(from = 0.002/365, to = 2/365, by = 0.00055)
beta
rho <- seq(from = 0.01, to = 0.1, by = 0.01)
rho
#epsilon <- seq(from = 2e-5 , to = 0.001, by = 0.0001)
#epsilon

dI <- as.data.frame(matrix(NA, length(beta),length(rho)))
res_min_I<-array(unlist(dI), dim=c(length(beta), length(rho)))
res_num_ext_I<-array(unlist(dI), dim=c(length(beta), length(rho)))
res_time_inf_I<-array(unlist(dI), dim=c(length(beta), length(rho)))

dN <- as.data.frame(matrix(NA, length(beta),length(rho)))
res_min_N<-array(unlist(dN), dim=c(length(beta), length(rho)))
res_num_ext_N<-array(unlist(dN), dim=c(length(beta), length(rho)))
res_time_inf_N<-array(unlist(dN), dim=c(length(beta), length(rho)))

############## FUNCTION  6 MODEL7 - brucellosis PAR LOOPS & MIN TIME TO EXTINCTION #######

my_fun_model7_N<-function() {
  for (k in 1:length(beta)){
    for (i in 1:length(rho)) {
      parameters_m7 <- c(beta_c = beta[k], beta_sa = beta[k], beta_a = beta[k],
                      phi_c = 1/14, phi_sa = 1/14, phi_a = 1/14,
                      gamma_c = 1/(2*365), gamma_sa = 1/(2*365), gamma_a = 1/(2*365),
                      rho_c = rho[i], rho_sa = rho[i], rho_a = rho[i], 
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
      initials_m7 <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, 
                       Sa = (a-1), Ea = 0, Ia = 1, Ra = 0)
      
      get_time <- model7(pars = parameters_m7, init = initials_m7,
                         end.time = end.time)
      res_min_N[k,i]<- my_min_ext_N(x=get_time$results,y=get_time$results$N)
      res_num_ext_N[k,i]<- my_imp_ext_N(x=get_time$time,y=get_time$results$N)
      res_time_inf_N[k,i]<- my_imp_ext_na_N(x=get_time$time,y=get_time$results$N)
      
    }}
  dN <- list(res_min_N,res_num_ext_N,res_time_inf_N)
  dN
}


my_fun_model7_I<-function() {
  for (k in 1:length(beta)){
    for (i in 1:length(rho)) {
      parameters_m7 <- c(beta_c = beta[k], beta_sa = beta[k], beta_a = beta[k],
                         phi_c = 1/14, phi_sa = 1/14, phi_a = 1/14,
                         gamma_c = 1/(2*365), gamma_sa = 1/(2*365), gamma_a = 1/(2*365),
                         rho_c = rho[i], rho_sa = rho[i], rho_a = rho[i], 
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
      initials_m7 <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, 
                       Sa = (a-1), Ea = 0, Ia = 1, Ra = 0)
      
      get_time <- model7(pars = parameters_m7, init = initials_m7,
                         end.time = end.time)
      # res_min[k,i]<- min(subset(get_time$results,I==0)$time)
      res_min_I[k,i]<- my_min_ext_I(x=get_time$results,y=get_time$results$I)
      res_num_ext_I[k,i]<- my_imp_ext_I(x=get_time$time,y=get_time$results$I)
      res_time_inf_I[k,i]<- my_imp_ext_na_I(x=get_time$time,y=get_time$results$I)
        
    }}
  dI <- list(res_min_I,res_num_ext_I,res_time_inf_I)
  dI
}
## PLOT DATA PREPARATION #################################

out_put_fun_min<-function(x, par1, par2, ...){ # x = out, par1, par2 = beta_a, rho_a, ...
  output <- array(unlist(x), dim = c(length(par1),length(par2),n_rep*3))
  output_min<-output[,,seq(from=1,to=n_rep*3,by=3)]
  output_min
}

out_put_fun_num<-function(x, par1, par2, ...){ # x = out, par1, par2 = beta_a, rho_a, ...
  output <- array(unlist(x), dim = c(length(par1),length(par2),n_rep*3))
  output_num<-output[,,seq(from=2,to=n_rep*3,by=3)]
  output_num
}

out_put_fun_time<-function(x, par1, par2, ...){ # x = out, par1, par2 = beta_a, rho_a, ...
  output <- array(unlist(x), dim = c(length(par1),length(par2),n_rep*3))
  output_time<-output[,,seq(from=3,to=n_rep*3,by=3)]/end.time
  output_time
}

############## FUNCTION 15 PREP OUTPUT FOR PLOTTING TIME TO EXTINCTION or NUMBER OF EXTINCTIONS #############

my_plot_min<-function(x, par1, par2, par1_n, par2_n){ # x = data run, e.g. big_run, par1, par2 are parameters varied, e.g. beta_a, rho_a, par1_n, par2_n are names, e.g. 'beta', 'rho', op is output, either time or extinctions
  plot_res<-apply(x, c(1,2), mean, na.rm = T)
  row.names(plot_res)<-par1; colnames(plot_res)<-par2
  df <- melt(plot_res)
  colnames(df)<-c(par1_n,par2_n,'duration')
  df
}

############## FUNCTION 16 PREP OUTPUT FOR PLOTTING PROPORTION OF OUTBREAKS PERSISTING ########

my_plot_num<-function(x, par1, par2, par1_n, par2_n){ # x = data run, e.g. big_run, par1, par2 are parameters varied, e.g. beta_a, rho_a
  x[!is.finite(x)] <- 0
  plot_res<-apply(x, c(1,2), sum, na.rm = T)/n_rep
  row.names(plot_res)<-par1; colnames(plot_res)<-par2
  df <- melt(plot_res)
  colnames(df)<-c(par1_n,par2_n,'outbreaks')
  df
}

############## FUNCTION 17 PREP OUTPUT FOR PLOTTING TIME TO EXTINCTION WITH NO INFINITE FROM PERSISTENCE ######

my_plot_time<-function(x, par1, par2, par1_n, par2_n){ # x = data run, e.g. big_run, par1, par2 are parameters varied, e.g. beta_a, rho_a
  x[!is.finite(x)] <- end.time
  plot_res<-apply(x, c(1,2), mean, na.rm = T)
  row.names(plot_res)<-par1; colnames(plot_res)<-par2
  df <- melt(plot_res)
  colnames(df)<-c(par1_n,par2_n,'time')
  df
}

############## FUNCTION 18 PREP OUTPUT FOR PLOTTING MX SIMULATIONS - SINGLE POP ####
single_pop_sim_prep_N <- function(x, n_rep, end.time){ # x = simulation of model, e.g. sim_run_m1
  mat = matrix(NA, nrow=n_rep, ncol = end.time+1)
  for (i in 1:n_rep){
    mat[i,]<-x[,i]$results$N
  }
  colnames(mat) = paste("time", seq(from=1,to=end.time+1,by=1), sep="")
  rownames(mat) = paste("run", seq(n_rep), sep="")
  dat = as.data.frame(mat)
  dat$run = rownames(dat)
  mdat = melt(dat, id.vars="run")
  mdat$time = as.numeric(gsub("time", "", mdat$variable))
  mdat
}

single_pop_sim_prep_I <- function(x, n_rep, end.time){ # x = simulation of model, e.g. sim_run_m1
  mat = matrix(NA, nrow=n_rep, ncol = end.time+1)
  for (i in 1:n_rep){
    mat[i,]<-x[,i]$results$I
  }
  colnames(mat) = paste("time", seq(from=1,to=end.time+1,by=1), sep="")
  rownames(mat) = paste("run", seq(n_rep), sep="")
  dat = as.data.frame(mat)
  dat$run = rownames(dat)
  mdat = melt(dat, id.vars="run")
  mdat$time = as.numeric(gsub("time", "", mdat$variable))
  mdat
}


n_rep = 10
end.time = 20 * 365

############## RUN MODEL 7 #####

big_run_model7_I<-replicate(n_rep,my_fun_model7_I())
big_run_model7_N<-replicate(n_rep,my_fun_model7_N())

## MODEL OUTPUTS PREPARATION #######################

############## MODELS 1-6 OUTPUT PREPARATION #############################
mod_res<-list(big_run_model7_I,
              big_run_model7_N)

for (i in 1:length(mod_res)){
  assign(paste0("Res_min_", i), out_put_fun_min(x=mod_res[[i]],par1 = beta,par2 = rho))
}

for (i in 1:length(mod_res)){
  assign(paste0("Res_num_", i), out_put_fun_num(x=mod_res[[i]],par1 = beta,par2 = rho))
}

for (i in 1:length(mod_res)){
  assign(paste0("Res_time_", i), out_put_fun_time(x=mod_res[[i]],par1 = beta,par2 = rho))
}


plot_res_min<-list(Res_min_1,
                   Res_min_2)

plot_res_min <- lapply(plot_res_min,function(x) replace(x,is.infinite(x),end.time))

plot_res_num<-list(Res_num_1,Res_num_2)
plot_res_time<-list(Res_time_1,Res_time_2)

plot_res_min <- lapply(plot_res_min,function(x) replace(x,is.infinite(x),end.time))

for (i in 1:length(plot_res_min)){
  assign(paste0("df_min_", i), my_plot_min(x=plot_res_min[[i]],par1 = beta, par2 = rho, par1_n = 'beta', par2_n = 'rho'))
}

for (i in 1:length(plot_res_num)){
  assign(paste0("df_num_", i), my_plot_num(x=plot_res_num[[i]],par1 = beta,par2 = rho, par1_n = 'beta', par2_n = 'rho'))
}

for (i in 1:length(plot_res_time)){
  assign(paste0("df_time_", i), my_plot_time(x=plot_res_time[[i]],par1 = beta,par2 = rho, par1_n = 'beta', par2_n = 'rho'))
}

n = length(plot_res_time)

for(i in 1:n){
  t1 <- do.call(cbind, mget(paste0("df_min_", 1:n) ) )
}
for(i in 1:n){
  t2 <- do.call(cbind, mget(paste0("df_num_", 1:n) ) )
}
for(i in 1:n){
  t3 <- do.call(cbind, mget(paste0("df_time_", 1:n) ) )
}

res<-cbind(t1,t2,t3)

colnames(res) <- colnames(res) %>% str_replace(".*.beta", "beta")
colnames(res) <- colnames(res) %>% str_replace(".*.rho", "rho")

res_all = melt(res, id.vars=c("beta",'rho'))

class(res_all)

## PLOT MODEL OUTPUTS TILES #############################

## PLOT MODEL OUTPUTS TILES #############################

r1d<- res_all %>% filter(stringr::str_detect(variable, "1.duration"))
str(r1d)
table(r1d$variable)
plot(r1d$value)
hist(r1d$value)
p<-(ggplot(r1d, aes(x = beta, y = rho, fill = value))+
      geom_tile()+
      scale_fill_gradientn(colors=colorRampPalette(c("whitesmoke","royalblue","seagreen","orange","red","brown"))(500),
                           name="Average\noutbreak\nduration", na.value = "grey") +
      labs(x = expression(beta),y=expression(rho)) +
      theme_bw())
p
ggsave("extinction_duration-bru_I.png",p,width = 10, height = 10, units = 'cm', dpi  = 600)

r1o<- res_all %>% filter(stringr::str_detect(variable, "1.outbreak"))
hist(r1o$value)
p<-(ggplot(r1o, aes(x = beta, y = rho, fill = value))+
      geom_tile()+
      scale_fill_gradientn(colors=colorRampPalette(c("brown","red","orange","seagreen","royalblue","whitesmoke"))(500),
                           name="Average\nnumber\nextinctions") +
      labs(x = expression(beta),y=expression(rho)) +
      theme_bw())
p
ggsave("extinction_outbreaks-bru_I.png",p,width = 10, height = 10, units = 'cm', dpi  = 600)

r1t<- res_all %>% filter(stringr::str_detect(variable, "1.time"))
str(r1t)
table(r1t$variable)
plot(r1t$value)
hist(r1t$value)

p<-(ggplot(r1t, aes(x = beta, y = rho, fill = value))+
      geom_tile()+
      scale_fill_gradientn(colors=colorRampPalette(c("brown","red","orange","seagreen","royalblue","whitesmoke"))(500),
                           name="Average\nnumber\nextinctions") +
      labs(x = expression(beta),y=expression(rho)) +
      theme_bw())

p
ggsave("extinction_outbreaks-bru_I.png",p,width = 10, height = 10, units = 'cm', dpi  = 600)



r2d<- res_all %>% filter(stringr::str_detect(variable, "2.duration"))
str(r2d)
table(r2d$variable)
plot(r2d$value) # N = 0 
hist(r2d$value)
p<-(ggplot(r2d, aes(x = beta, y = rho, fill = value))+
      geom_tile()+
      scale_fill_gradientn(colors=colorRampPalette(c("whitesmoke","royalblue","seagreen","orange","red","brown"))(500),
                           name="Average\noutbreak\nduration", na.value = "grey") +
      labs(x = expression(beta),y=expression(rho)) +
      theme_bw())
p



r2o<- res_all %>% filter(stringr::str_detect(variable, "2.outbreaks"))
str(r2o)
table(r2o$variable)
plot(r2o$value)
p<-(ggplot(r2o, aes(x = beta, y = rho, fill = value))+
      geom_tile()+
      scale_fill_gradientn(colors=colorRampPalette(c("brown","red","orange","seagreen","royalblue","whitesmoke"))(500),
                           name="Average\nnumber\nextinctions") +
      labs(x = expression(beta),y=expression(rho)) +
      theme_bw())

p

r2t<- res_all %>% filter(stringr::str_detect(variable, "2.time"))
str(r2t)
table(r2t$variable)
plot(r2t$value)
hist(r2t$value)
p<-(ggplot(r2t, aes(x = beta, y = rho, fill = value))+
      geom_tile()+
      scale_fill_gradientn(colors=colorRampPalette(c("brown","red","orange","seagreen","royalblue","whitesmoke"))(500),
                           name="Average\nnumber\nextinctions") +
      labs(x = expression(beta),y=expression(rho)) +
      theme_bw())

p

## SINGLE POPULATION PRINT PREPS #######

res_mx<-list(sim_run_m1_1,
             sim_run_m1_2,
             sim_run_m1_3,
             sim_run_m2_1,
             sim_run_m2_2,
             sim_run_m3,
             sim_run_m4,
             sim_run_m5,
             sim_run_m7_1,
             sim_run_m7_2)

for (i in 1:length(res_mx)){
  assign(paste0("df_mx", i), single_pop_sim_prep_N(x=res_mx[[i]],n_rep = n_rep, end.time = end.time))
}

for (i in 1:length(res_mx)){
  assign(paste0("df_mx", i), single_pop_sim_prep_I(x=res_mx[[i]],n_rep = n_rep, end.time = end.time))
}


res_mx_p<-rbind(df_mx1,
                df_mx2,
                df_mx3,
                df_mx4,
                df_mx5,
                df_mx6,
                df_mx7,
                df_mx8,
                df_mx9,
                df_mx10)

res_mx_p$model<-c(rep('1',dim(df_mx1)[1]),
                  rep('2',dim(df_mx2)[1]),
                  rep('3',dim(df_mx3)[1]),
                  rep('4',dim(df_mx4)[1]),
                  rep('5',dim(df_mx5)[1]),
                  rep('6',dim(df_mx6)[1]),
                  rep('7',dim(df_mx7)[1]),
                  rep('8',dim(df_mx8)[1]),
                  rep('9',dim(df_mx9)[1]),
                  rep('10',dim(df_mx10)[1]))

## SINGLE POPULATION PLOTS #######

for (i in unique(res_mx_p$model)){
  subdata <- subset(res_mx_p, model == i)
  pdf(paste("plot_ts_mx", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot(s2, aes(x=time, y=N, group=run)) +
          theme_bw() +
          theme(panel.grid=element_blank()) +
          geom_line(size=0.2, alpha=0.15)+
          ylab('Numbers') + xlab('time')+
          stat_summary(aes(group = 1), fun=mean, geom="line", colour="black",size = 1.1))
  dev.off()
}

############### PLOT ALL TS SINGLE RUNS ##############

res_p_ts<-list(sim_run_m1_1,
               sim_run_m1_2,
               sim_run_m1_3,
               sim_run_m2_1,
               sim_run_m2_2,
               sim_run_m3,
               sim_run_m4,
               sim_run_m5,
               sim_run_m6,
               sim_run_m6v,
               sim_run_m7_1,
               sim_run_m7_2)

# Make plots. Single runs
for (i in 1:length(res_p_ts)) {
  pdf(paste("plotts", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot() + 
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$S), color = "black", size =1.2) +
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$I), color = "red", size =1.2) +
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$R), color = "seagreen", size =1.2, linetype = "dotted") +
          xlab('Time') +
          ylab('Infection State Numbers'))
  dev.off()
}

############### PLOT ALL TS SINGLE RUNS I ONLY ##############

# Make plots. Single runs
for (i in 1:length(res_p_ts)) {
  pdf(paste("plotts_i", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot() + 
          geom_line(data = res_p_ts[[i]][[4]], aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$I), color = "red", size =1.2) +
          xlab('Time') +
          ylab('Infection Numbers'))
  dev.off()
}

# Make plots. Single runs - SCALED
for (i in 1:length(res_p_ts)) {
  pdf(paste("plotts_i_scale", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot() + 
          geom_line(data = res_p_ts[[i]][[4]], 
                    aes(x = res_p_ts[[i]][[4]]$time, y = res_p_ts[[i]][[4]]$I), color = "red", size =1.2) +
          xlab('Time') +
          ylim(0, 500) + 
          ylab('Infection Numbers'))
  dev.off()
}

## plot extinction times

df.agg <- aggregate(time ~ run + value + model, res_mx_p, min)

df.ag <- (df.agg[df.agg$value==0,c('model','time')])

neworder <- c("1","2","3","4","5","9","10","6","7","8")
library(plyr)  ## or dplyr (transform -> mutate)
df.ag <- arrange(transform(df.ag,
                           model=factor(model,levels=neworder)),model)
labs <- c('1' = "SIR",
          '2' = "SI*R",
          '3' = "SI*R+",
          '4' = "S[I]R",
          '5' = "S[I]*R",
          '6' = "SIR*",
          '7' = 'S[I]*R',
          '8' = 'S[I]R*',
          '9' = 'SEIR',
          '10' = 'SEI*R')

p<-ggplot(df.ag, aes(x=time))+
  geom_histogram(color="black", fill="grey")+
  facet_wrap(model~., ncol = 5,labeller = labeller(model = labs))+
  scale_x_continuous(breaks = c(0, 800, 1600), labels = c("0", "800", "1600"))
pdf("extinctions.pdf", width = 5, height = 3)
p
dev.off()

p_yr<-ggplot(df.ag, aes(x=time))+
  geom_histogram(color="black", fill="grey")+
  facet_wrap(model~., ncol = 5,labeller = labeller(model = labs))+
  scale_x_continuous(limits = c(0,1000),breaks = c(0, 400, 800), labels = c("0", "400", "800")) +
  ylim(0,60)
pdf("extinctions_1yr.pdf", width = 5, height = 3)
p_yr
dev.off()
