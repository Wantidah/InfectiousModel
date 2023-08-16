# creating population dynamics models: containing 5 species parameter values

############## LOAD PACKAGES ##########
library(EpiDynamics)
library(dplyr)  
library(plyr)
library(tidyverse)
library(reshape2) 
library(stringr)
library(ggplot2); theme_set(theme_bw())

set.seed(111)
# Population dynamic model, no infection, 3 age classes
# c = calf
# sa = Subadult
# a = adult
# unit == day (per day)

############## MODEL 1 population dynamic model, no infection, 3 age-classes 
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
          num.min <- min(num, init[which(change[i, ] <0)])
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
    results<-data.frame(time, a,  sa,  c)%>% 
      dplyr::mutate(N = rowSums(across(-c(time)), na.rm=TRUE))
    
    return(list(pars = pars, init = init2, time = time, results = results))

}
# setting parameter values     
#predicting time  (unit = days)
# this can be adjusted, in the manuscript used end.time 100*365, n_rep = 100
end.time <- 2*365 
n_rep <- 2
# gaur population #########
N = 300 
# calf:subadult:adult ratio
rat = 1.3+1.3+1.5
N/rat
c = round((N/rat)*1.3, 0)
sa = round((N/rat)*1.3,0) 
a = round((N/rat)*1.5,0) 

initials1 <- c(c = c, sa = sa, a = a )
pm1<- c(mu_b = 0.34/365, 
                mu_c = 0.27/365, 
                mu_sa = 0.15/365,
                mu_a = 0.165/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials1), 
                tau = 1)
# banteng population #########
N = 290
#calf:subadult:adult ratio
rat = 1.1+1.3+1
c = round((N/rat)*1.1,0)
sa = round((N/rat)*1.3, 0)
a = round((N/rat*1), 0)
initials2 <- c(c = c, sa = sa, a = a )
pm2 <- c(mu_b = 0.4/365, 
                    mu_c = 0.26/365, 
                    mu_sa = 0.26/365,
                    mu_a = 0.15/365,
                    delta_c = 1/365,
                    delta_sa = 1/(3*365),
                    N = sum(initials2), 
                    tau = 1)
#buffalo population #########
N = 70
#calf:subadult:adult ratio
rat = 1+6+5
N/rat
c = round((N/rat)*1,0)
sa = round((N/rat)*6,0)
a = round((N/rat)*5,0)
initials3 <- c(c = c, sa = sa, a = a )
pm3 <- c(mu_b = 0.37/365, 
                    mu_c = 0.27/365, 
                    mu_sa = 0.15/365,
                    mu_a = 0.15/365,
                    delta_c = 1/365,
                    delta_sa = 1/(3*365),
                    N = sum(initials3), 
                    tau = 1)
           
#serow population  #########
N = 120 
#calf:subadult:adult ratio
rat = 1+1+1
N/rat
c = round((N/rat)*1,0)
sa = round((N/rat)*1,0)
a = round((N/rat)*1,0)
initials4 <- c(c = c, sa = sa, a = a )
pm4<- c(mu_b = 0.7/365, 
                   mu_c = 0.5/365, 
                   mu_sa = 0.15/365,
                   mu_a = 0.28/365,
                   delta_c = 1/365,
                   delta_sa = 1/(3*365),
                   N = sum(initials4), 
                   tau = 1)
           
# goral population  #########
# assume
N = 292
rat = 1+1+1
N/rat
c = round((N/rat)*1,0)
sa = round((N/rat)*1,0)
a = round((N/rat)*1,0)
initials5 <- c(c = c, sa = sa, a = a )

pm5<- c(mu_b = 0.5/365, 
                   mu_c = 0.45/365, 
                   mu_sa = 0.27/365,
                   mu_a = 0.17/365,
                   delta_c = 1/365,
                   delta_sa = 1/(3*365),
                   N = sum(initials5), 
                   tau = 1)
           
# list of paramters
init<-list(initials1,initials2,initials3,initials4,initials5)
pm <- list (pm1,pm2,pm3,pm4,pm5)

dfs<-list()
dfm<-list()

nam<-c('Gaur', 'Banteng', 'Buffalo','Serow','Goral')

#  single, multiple run -------------
for (i in 1:length(pm)) {
  dfs[[i]]<-model1(pars = pm[[i]], init = init[[i]],
                  end.time = end.time)
  #PlotMods(df[[i]]) #Epi PlotMods
  dfm[[i]] <- replicate(n_rep,(model1( pm[[i]], init = init[[i]],
                                              end.time = end.time)))
}
#plot sum of populations

for (i in 1:length(pm)) {
plot(dfs[[i]]$results$N, 
     main = paste0(nam[[i]]," ","total population"), 
     xlab="time",ylab="numbers") 
}

#convert to data.frame, change days -> years
dfs2<-list()
for(i in 1:length(dfs)){
dfs2[[i]]<-df[[i]]$results %>%
  mutate(time_y = time/365) %>% #convert day to year for plotting
  melt(id.vars = c('time','time_y'), 
       value.name = 'value', variable.name = 'class')
 
dfs2[[i]]$species <- paste0(nam[[i]]) # adding the model name 
}

# summarise min-max by class
for(i in 1:length(dfs2)){
  print(dfs2[[i]] |> group_by(class,species) |>
          dplyr::summarise(Max=max(value),
                           Min=min(value))) 
}

# FUNCTION adding time and run number for multiple model's simulations ----
pop_sim_prep <- function(x, n_rep, end.time, melt){ # x = simulation of model, e.g. sim_run_m1
  df<-list()
  for (i in 1:n_rep){
    run <- paste("run", seq(n_rep), sep="")
    names(df)[i]
    df[[i]] <- x[,i]$results[,-c(1)]
    df[[i]]$time_d <- seq(from=1,to=end.time+1,by=1)
  }
  
  df <- map2(df,run, ~cbind(.x, run = .y))   # adding number of replications to the column
  df2 <- data.table::rbindlist(df)          # binding row
  
  if (melt == T) {  #option for melting the data in case we need...
    
    df3 <- melt(df2, id.vars = c('time_d','run'))
    return(df3 = data.frame(df3))
  } 
  
  else  {
    return( df2 = data.frame(df2)) 
  }
  
}

# rearrange df and calculate total population change (%) --------
sm <- list()
for (i in 1:length(pm)) {
  sm[[i]] <- pop_sim_prep(x = dfm[[i]], n_rep=n_rep, end.time= end.time, melt = F)
  sm[[i]] <- sm[[i]]%>%
    group_by(run) %>%
    mutate(Ndiff = ((last(N)-first(N))/first(N))*100)%>% #calculate change in the total population at year100, and year0
    mutate(time_y = time_d/365) %>% #convert day to year for plotting
    as.data.frame()
  sm[[i]]$model <- paste0(nam[[i]])
  
}
# plot population dynamics ######

# 1 run
pp<-list()

for(i in 1:length(dfs2)){
pp[[i]]<-ggplot(dfs2[[i]]) + 
  geom_line(aes(x = time_y ,y = value,  color = class), linewidth=1)  +
  labs(x="Years", y= "Compartment", title= nam[[i]]) +
  #scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "Population",
                    labels = c( 'adult','subadult','calf','total' ),
                      values = c('a'='seagreen4',
                                 'sa'='firebrick',
                                 'c'='dodgerblue3',
                                 'N'='#153030'))+ 
  
  my_theme+
  theme(legend.position = "none")   
 #ggsave(paste0(nam[[i]],"_pop_1sim_test.png"),pp[[i]], width = 25, height = 15, units = 'cm', dpi = 600)
}
print(pp)

pp2<-list()

for (i in 1:length(sm)) {
  pp2[[i]]<-ggplot(sm[[i]]) + 
    geom_line(aes(x = time_y ,y = a,  group = run, color = 'adult'), linewidth = 0.1, alpha = 0.12) + 
    geom_line(aes(x = time_y, y = sa, group = run, color = 'subadult'),linewidth = 0.1, alpha = 0.12) +
    geom_line(aes(x = time_y, y = c,  group = run, color = 'calf' ), linewidth = 0.1, alpha = 0.12)+
    geom_line(aes(x = time_y, y = N,  group = run, color = 'total'), linewidth = 0.1, alpha = 0.12) +
    labs(x="years", y= "population", title= paste0(LETTERS[1:5][i],") ", nam[[i]])) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "Compartment",
                        labels = c( 'adult','subadult','calf','total'),
                        values = c('seagreen4',
                                   'firebrick',
                                   'dodgerblue3',
                                   '#153030'),
                        breaks = c('adult','subadult','calf','total'))+ 
    my_theme +
    
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth =0.7 )))+
    
    stat_summary(sm[[i]], mapping =aes( x = time_y, y = a, group = 1), fun=mean, geom="line", colour='seagreen4',linewidth = 0.5)+
    stat_summary(sm[[i]], mapping =aes( x = time_y, y = sa, group = 1), fun=mean, geom="line", colour= 'firebrick',linewidth = 0.5)+
    stat_summary(sm[[i]], mapping =aes( x = time_y, y = c, group = 1), fun=mean, geom="line", colour= 'dodgerblue3',linewidth = 0.5)+
    stat_summary(sm[[i]], mapping =aes( x = time_y, y = N, group = 1), fun=mean, geom="line", colour= '#153030',linewidth = 0.5)
}

print(pp2)

ps1 <- wrap_plots(pp, ncol=1) & plot_annotation(title = '1 simulation') & theme(plot.title = element_text(hjust = 0.01))
ps2 <- wrap_plots(pp2, ncol=1) & plot_annotation(title = '100 simulations') & theme(plot.title = element_text(hjust = 0.01))
ps3<-(wrap_elements(pp1)|wrap_elements(pp2))
ps4<-wrap_elements(pp3 + plot_layout(widths = c(0.77,0.96)))
ps4
ggsave("s_pop_5sp.png",ps4, width = 25, height = 30, units = 'cm', dpi = 600)

#-- DONE :) -- #