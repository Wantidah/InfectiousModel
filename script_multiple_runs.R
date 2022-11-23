###### Multiple runs ###########

#### set initial values - change for populations
N = 300 

rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

initials <- c(c = c, sa = sa, a = a )

end.time <- 20 * 365

n_rep <- 100

############## MODEL 1 - POPULATION DEMO NO INFECTION #####

parameters <- c(mu_b = 0.34/365, 
                mu_c = 0.27/365, 
                mu_sa = 0.15/365,
                mu_a = 0.165/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials), 
                tau = 1)
sim_run_non_rep<-replicate(n_rep,(model1(pars = parameters, init = initials,
                                             end.time = end.time)))

############## MODEL 2 SI Anthrax MULTIPLE RUNS  #####
sim_run_anthrax_rep<-replicate(n_rep,(model2(pars = parameters, init = initials,
                                             end.time = end.time)))