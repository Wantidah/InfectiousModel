###### Multiple runs ###########

#### set initial values - change for populations
N = 300 

rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

end.time <- 20 * 365
n_rep <- 10

############## MODEL 1 - POPULATION DEMO - NO INFECTION #####
parameters <- c(mu_b = 0.34/365, 
                mu_c = 0.27/365, 
                mu_sa = 0.15/365,
                mu_a = 0.165/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials), 
                tau = 1)
initials <- c(c = c, sa = sa, a = a )

sim_model1_rep<-replicate(n_rep,(model1(pars = parameters, init = initials,
                                             end.time = end.time)))
sim_model1_rep

############## MODEL 2 Anthrax MULTIPLE RUNS  #####
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

initials <- c(Sc = c, Ic = 0, Ssa = sa, Isa = 0, Sa = a, Ia = 1 )

sim_run_anthrax_rep<-replicate(n_rep,(model2(pars = parameters, init = initials,
                                             end.time = end.time)))

############## MODEL 3 Bovine tuberculosis  MULTIPLE RUNS  #####
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
  mu_b = 0.34/365, 
  mu_bI = (0.34/365)*(1-0.27), #Ia birth rate reduce by = 27%   
  mu_c = 0.27/365, 
  mu_sa = 0.15/365,
  mu_a = 0.165/365,
  delta_c = 1/365,
  delta_sa = 1/(3*365),  
  epsilon = 2e-5,
  N = sum(initials),
  tau=1
)

initials <- c(Sc = c, Ec = 0, Ic = 0, Ssa = sa, Esa = 0, Isa = 0, Sa = a, Ea = 0, Ia = 1)

sim_run_tb_rep<-replicate(n_rep,(model3(pars = parameters, init = initials,
                                             end.time = end.time)))

############## MODEL 4  Hemorrhagic septicemia   MULTIPLE RUNS  #####
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

initials <- c(Sc = c,  Ic = 0, Rc = 0, Ssa = sa, Isa = 0, Rsa = 0, Sa = a, Ia = 1, Ra = 0)
sim_run_hs_rep<-replicate(n_rep,(model4(pars = parameters, init = initials,
                                        end.time = end.time)))

############## MODEL 5  Lumpy Skin Disease   MULTIPLE RUNS  #####
parameters <- c( 
  beta_c = 0.043/30,
  beta_sa = 0.043/30,
  beta_a = 0.043/30,
  phi_c = 1/7/365,
  phi_sa = 1/7/365,
  phi_a = 1/7/365,
  gamma_c = 0,
  gamma_sa = 0,
  gamma_a = 0,
  rho_c = 0.03,
  rho_sa = 0.03,
  rho_a = 0.03, 
  omega_c = (1/365)/365,
  omega_sa =  (1/365)/365, 
  omega_a = (1/365)/365,
  epsilon = 2e-5,
  mu_b = 0.34/365, 
  mu_bI = (0.34/365)*(1-0.27), #Ia birth rate reduce by = 27%   
  mu_c = 0.27/365, 
  mu_sa = 0.15/365,
  mu_a = 0.165/365,
  delta_c = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials),
  tau=1
)

initials <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = a, Ea = 0, Ia = 1, Ra = 0)

sim_run_lsd_rep<-replicate(n_rep,(model5(pars = parameters, init = initials,
                                        end.time = end.time)))

############## MODEL 6 Foot and mouth disease  MULTIPLE RUNS  #####
parameters <- c( 
  beta_c = 0.52/365,
  beta_sa = 0.52/365,
  beta_a = 0.52/365,
  phi_c = 1/8/365,
  phi_sa = 1/6/365,
  phi_a = 1/6/365,
  gamma_c = 1/5/365,
  gamma_sa = 1/5/365,
  gamma_a = 1/5/365,
  rho_c = 0.1,
  rho_sa = 0.05,
  rho_a = 0.03, 
  alpha = 0.5,
  omega_c = (1/120)/365,
  omega_sa =  (1/120)/365, 
  omega_a = (1/565)/365,
  omega_m = (1/144)/365,
  epsilon = 2e-5,
  mu_b = 0.34/365, 
  mu_bI = (0.34/365)*(0.9), #Ia birth rate reduce by = 10%  (assume)
  mu_c = 0.27/365, 
  mu_sa = 0.15/365,
  mu_a = 0.165/365,
  delta_c = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials),
  tau=1
)

initials <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = a, Ea = 0, Ia = 1, Ra = 0)

sim_run_fmd_rep<-replicate(n_rep,(model6(pars = parameters, init = initials,
                                         end.time = end.time)))

############## MODEL 7 Bovine brucellosis MULTIPLE RUNS  #####
parameters <- c( 
  beta_c = 2/365,
  beta_sa = 2/365,
  beta_a = 2/365,
  phi_c = 1/14/365,
  phi_sa = 1/14/365,
  phi_a = 1/14/365,
  gamma_c = 1/2*365,
  gamma_sa = 1/2*365,
  gamma_a = 1/2*365,
  rho_c = 0.1,
  rho_sa = 0.05,
  rho_a = 0.03, 
  alpha = 0.9,
  omega_c = (1/180)/365,
  omega_sa =  (1/180)/365, 
  omega_a = (1/180)/365,
  omega_m = (1/180)/365,
  epsilon = 2e-5,
  mu_b = 0.34/365, 
  mu_bI = (0.34/365)*(0.5), #Ia birth rate reduce by = 10%  (assume)
  mu_c = 0.27/365, 
  mu_sa = 0.15/365,
  mu_a = 0.165/365,
  delta_c = 1/365,
  delta_sa = 1/(3*365),
  N = sum(initials),
  tau=1
)
initials <- c(Sc = c, Ec = 0, Ic = 0, Rc = 0, M = 0, Sm = 0, Ssa = sa, Esa = 0, Isa = 0, Rsa = 0, Sa = a, Ea = 0, Ia = 1, Ra = 0)

sim_run_bru_rep<-replicate(n_rep,(model7(pars = parameters, init = initials,
                                         end.time = end.time)))
###### FUNCTION single_pop_sim_prep
#change from results$I to results$Ia
#However, this Ia need to be changed into N (total population) later
single_pop_sim_prep <- function(x, n_rep, end.time){ # x = simulation of model, e.g. sim_run_m1
  mat = matrix(NA, nrow=n_rep, ncol = end.time+1)
  for (i in 1:n_rep){
    mat[i,]<-x[,i]$results$Ia #Dave's one is results$I
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
###########list simulations multiple runs

res_mx<-list(
             sim_run_anthrax_rep,
             sim_run_tb_rep,
             sim_run_hs_rep,
             sim_run_lsd_rep,
             sim_run_fmd_rep,
             sim_run_bru_rep)

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
