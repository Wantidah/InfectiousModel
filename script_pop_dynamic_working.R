# Infectious disease modelling code
# - 18 Nov 2022-
# # Reference R Code:  https://doi.org/10.1098/rsif.2021.0638  (Transmission models indicate Ebola virus persistence in non-human primate populations is unlikely)
# https://github.com/dtsh2/ebola_model

rm(list=ls())
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
set.seed(111)
############## Run the model function ############### 
############## MODEL 1 population dynamic model, no infection, 3 age-classes ######
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
  
# TEST MODEL 1 Population dynamic, no infection ######
# set initial values
end.time <- 100*365 #predict for ... years
    
#gaur population #########
  N = 300 
  #calf:subadult:adult ratio
  rat = 1.3+1.3+1.5
  N/rat
  c = round((N/rat)*1.3, 0)
  sa = round((N/rat)*1.3,0) 
  a = round((N/rat)*1.5,0) 
#same initials population for every species
initials_m1 <- c(c = c, sa = sa, a = a )

pm1<- c(mu_b = 0.34/365, 
                mu_c = 0.27/365, 
                mu_sa = 0.15/365,
                mu_a = 0.165/365,
                delta_c = 1/365,
                delta_sa = 1/(3*365),
                N = sum(initials_m1), 
                tau = 1)
# banteng population #########
N = 290
#calf:subadult:adult ratio
rat = 1.1+1.3+1
c = N/rat*1.1
sa = N/rat*1.3
a = N/rat*1
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
c = (N/rat)*1
sa = (N/rat)*6
a = (N/rat)*5
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

#goral population  #########
#assume
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

# test plotting
for (i in length(pm)) {
 
  df[[i]]<-model1(pars = pm[[i]], init = init[[i]],
                  end.time = end.time)
  }

  
  #plot sum of populations
    plot(df[[i]]$N, main = paste0(nam[[i]]," ","total population"), 
       xlab="time",ylab="numbers")

df

# TEST single run -------------
res_model1 <- model1(pars = pm1, init = initials_m1,
                     end.time = end.time)
str(res_model1)

#plot epi
png("gaur_pop_1sim_epiplot.png",width = 25, height = 15, units = 'cm', res = 600)
PlotMods(res_model1)
dev.off()

#minimum N extinction
min(subset(res_model1$results,N==0)$time)

#convert to data.frame, change days -> years
r1<-res_model1$results %>%
  mutate(time_y = time/365) %>% #convert day to year for plotting
  melt(id.vars = c('time','time_y'), 
       value.name = 'value', variable.name = 'class')

# the df may need to be in an order we want. Just for easier when plotting the line graph 
# in a class column it should be: adult, subadult, calf, total
str(r1)

#check min, max, mean for all class, and we can set the  limit of y scale
r1 |> group_by(class) |>
  summarise(Mean = mean(value),
            Max=max(value),
            Min=min(value))

# plot population dynamic ######
p1<-ggplot(r1) + 
  geom_line(aes(x = time_y ,y = value,  color = class))  +
 
  labs(x="years", y= "population",
       title= 'No infection') +
  
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
  guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))

print(p1)

ggsave("gaur_pop_1sim.png",p1, width = 25, height = 15, units = 'cm', dpi = 600)

#scale: ylim(0, 1000)
p1s<-ggplot(r1) + 
  geom_line(aes(x = time_y ,y = value,  color = class))  +
  
  labs(x="years", y= "population",
       title= 'A) no infection') +
  
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  ylim(0, 1000)+
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
  guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))

print(p1s)

ggsave("gaur_pop_1sim_scale.png",p1s, width = 25, height = 15, units = 'cm', dpi = 600)


