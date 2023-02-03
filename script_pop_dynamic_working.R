# creating population dynamics models: containing 5 species parameter values

############## LOAD PACKAGES ##########
library(EpiDynamics)
library(dplyr)  
library(plyr)
library(tidyverse)
library(reshape2) 
library(stringr)
library(emdbook)  
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
      dplyr::mutate(N = rowSums(across(-c(time), na.rm=TRUE)))
    
    return(list(pars = pars, init = init2, time = time, results = results))

}
# setting parameter values     
#predicting time  (unit = days)
end.time <- 100*365 

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
df<-list()
nam<-c('Gaur', 'Banteng', 'Buffalo','Serow','Goral')

# TEST single run -------------
for (i in 1:length(pm)) {
  df[[i]]<-model1(pars = pm[[i]], init = init[[i]],
                  end.time = end.time)
  #PlotMods(df[[i]]) #Epi PlotMods
}
#plot sum of populations
p<-list()
for (i in 1:length(pm)) {
p[[i]]<-plot(df[[i]]$results$N, 
     main = paste0(nam[[i]]," ","total population"), 
     xlab="time",ylab="numbers") }

print(p)

#convert to data.frame, change days -> years
df2<-list()
for(i in 1:length(df)){
df2[[i]]<-df[[i]]$results %>%
  mutate(time_y = time/365) %>% #convert day to year for plotting
  melt(id.vars = c('time','time_y'), 
       value.name = 'value', variable.name = 'class')
 df2[[i]]$species <- paste0(nam[[i]]) # adding the model name 
  str(df2[[i]])}

# summarise min-max by class
for(i in 1:length(df2)){
  print(df2[[i]] |> group_by(class,species) |>
          dplyr::summarise(Max=max(value),
                           Min=min(value))) 
  }

# plot population dynamics ######
pp<-list()
for(i in 1:length(df2)){
pp[[i]]<-ggplot(df2[[i]]) + 
  geom_line(aes(x = time_y ,y = value,  color = class), linewidth=1)  +
  labs(x="years", y= "population",
       title= nam[[i]]) +
  #scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                    labels = c( 'adult','subadult','calf','total' ),
                      values = c('a'='seagreen4',
                                 'sa'='firebrick',
                                 'c'='dodgerblue3',
                                 'N'='#153030'))+ 
  
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth=1)))
 #ggsave(paste0(nam[[i]],"_pop_1sim_test.png"),pp[[i]], width = 25, height = 15, units = 'cm', dpi = 600)
}
print(pp)

#scale: ylim(0, 500)
ps<-list()
for(i in 1:length(df2)){
  ps[[i]]<-ggplot(df2[[i]]) + 
    geom_line(aes(x = time_y ,y = value,  color = class), linewidth=1)  +
    labs(x="years", y= "population",
         title= nam[[i]]) +
    ylim(0, 500)+
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "population",
                        labels = c( 'adult','subadult','calf','total' ),
                        values = c('a'='seagreen4',
                                   'sa'='firebrick',
                                   'c'='dodgerblue3',
                                   'N'='#153030'))+ 
    
    theme_bw() +
    theme( plot.title = element_text(size = 13),
           axis.title.x = element_text(size = 12),
           axis.title.y = element_text(size = 12),
           legend.title=element_text(size=11),
           legend.text = element_text(size = 11),
           axis.text=element_text(size=11))+
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth=1)))
  #ggsave(paste0(nam[[i]],"_pop_1sim_sacle_test.png"),ps[[i]], width = 25, height = 15, units = 'cm', dpi = 600)
}
print(ps)

#-- DONE :) -- #