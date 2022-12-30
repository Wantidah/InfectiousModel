# try to make dataframe with all the population class
# ongoning
# the for loop is works and can make a dataframe we want, but kust have to find a nicer way of plotting graph

############## LOAD PACKAGES ##########
library(EpiDynamics)
library(plyr)     
library(reshape2) 
library(stringr)
library(emdbook)  
library(ggplot2); theme_set(theme_bw())
library(dplyr)   
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggstatsplot)

n_rep = 5
end.time = 365

m6<-replicate(n_rep,(model6(pars = parameters_m6, init = initials_m6,
                                    end.time = end.time)))

# TEST: single_pop_sim_prep Dave's adjust version
single_pop_sim_prep <- function(x, n_rep, end.time, melt){ # x = simulation of model, e.g. sim_run_m1
  df<-list()
  #loop for storing new df
  for (i in 1:n_rep){
    run <- paste("run", seq(n_rep), sep="")
    names(df)[i]
    df[[i]]<- x[,i]$results[,-c(1)]
    df[[i]]$time_d <- seq(from=1,to=end.time+1,by=1)
  }
    df<-map2(df,run, ~cbind(.x, run = .y))   # adding n_rep to the column
    df2<- data.table::rbindlist(df)          # binding row
    
    if (melt == T) {  #option for melting the data in case we need...
     
    #df3 <- gather(df2, key = class, value = value, -c(time_d,run))
      df3 <- melt(df2, id.vars = c('time_d','run'))
    return(df3 = data.frame(df3))
     } 
   
     else  {
      return( df2 = data.frame(df2)) 
     }

   }


sim_rep_m<-list(sim_rep_m1,
                sim_rep_m2,
                sim_rep_m3,
                sim_rep_m4,
                sim_rep_m5,
                sim_rep_m6,
                sim_rep_m7)

m<-list()

nam<-c('no_infection',
       'Anthrax',
       'bTB',
       'HS',
       'LSD',
       'FMD',
       'Brucellosis')

#group and calculate total population change (%) loop
for (i in 1:length(sim_rep_m)) {
  m[[i]]<- single_pop_sim_prep(x = sim_rep_m[[i]], n_rep=n_rep, end.time= end.time, melt = F) 
 
  m[[i]]<- m[[i]]%>%
    group_by(run) %>%
    mutate(Ndiff = ((last(N)-first(N))/first(N))*100)%>% #calculate change in the total population at year100, and year0
    mutate(time_y = time_d/365) %>% #convert day to year for plotting
    as.data.frame()
  
  m[[i]]$model <- paste0(nam[[i]]) # adding name
  
  }

res_melt<-list()

for (i in 1:length(m)){

  res_melt[[i]] <- melt(m[[i]], id.vars = c('time_y','time_d','run','model'),
           value.name = 'value', variable.name = 'class')
}

str(res_melt)

nam2<-c('m1_noinf',
        'm2_anthrax',
        'm3_bTB',
        'm4_HS',
        'm5_LSD',
        'm6_FMD',
        'm7_brucellosis')

for (i in 1:length(res_melt)) {
  saveRDS(res_melt[[i]], file = paste0("df_",nam2[[i]],"_100runs.rds")) }

# PLOT MODEL OUTPUTS: population dynamic line graphs #############################
# multiple run plots
head(lm[[1]])
 
pl1<-ggplot() + 
      geom_line(data = lm[[1]],aes(x = time_y ,y = a,  color = 'adult'),size = 0.1, alpha = 0.12) + 
      geom_line(data = lm[[1]],aes(x = time_y, y = sa, color = 'subadult' ),size = 0.1, alpha = 0.12) +
      geom_line(data = lm[[1]],aes(x = time_y, y = c,  color = 'calf' ),size = 0.1, alpha = 0.12)+
      geom_line(data = lm[[1]], aes(x = time_y, y = N, color = 'total'),size = 0.1, alpha = 0.12) +
   labs(x="years", y= "population",
        title= 'H) No infection') +
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
  guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
  stat_summary(lm[[1]], mapping =aes( x = time_y, y = a, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
  stat_summary(lm[[1]], mapping = aes( x = time_y, y = sa, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
  stat_summary(lm[[1]], mapping = aes( x = time_y, y = c, group = 1), fun=mean, geom="line", colour="darkorange2",size = 0.5)+
  stat_summary(lm[[1]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",size = 0.5) #blackgreen

  print(pl1)
    
    ggsave("gaur_no_infection_100sim.png",pl1, width = 25, height = 15, units = 'cm', dpi = 600)

  p2<-ggplot() + 
    geom_line(data = s[[i]],aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
    geom_line(data = s[[i]],aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
    geom_line(data = s[[i]], aes(x = time_y, y = N,group = run, color = 'total'), size = 0.1, alpha = 0.12)+
    
    labs(x="years", y= "population",
         title= 'Gaur population with anthrax infection, 1 simulation') +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "population",
                        labels = c('S','I','total' ),
                        values = c('S'='seagreen4',
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
    
    print(p2)
    #ggsave("gaur_anthrax_100sim.png",p,width = 25, height = 15, units = 'cm', dpi = 600)
}

## m3 - bTB plot ######  
else if (s[[i]]$model=="bTB") {
  p3<-ggplot() + 
    geom_line(data = s[[i]],aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
    geom_line(data = s[[i]],aes(x = time_y, y = E, group = run, color = 'E' ),size = 0.1, alpha = 0.12)+
    geom_line(data = s[[i]],aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
    geom_line(data = s[[i]], aes(x = time_y, y = N,group = run, color = 'total'), size = 0.1, alpha = 0.12)+
    
    labs(x="years", y= "population",
         title = paste0('Gaur population with bTB infection, 100 simulations')) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "population",
                        labels = c('S','E','I','total' ),
                        values = c('S'='seagreen4',
                                   'E'='darkorange2',
                                   'I'='firebrick',
                                   "total"='#153030'))+ #blackgreen
    theme_bw() +
    theme( plot.title = element_text(size = 18),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 15),
           legend.title=element_text(size=11),
           legend.text = element_text(size = 11),
           axis.text=element_text(size=13))+
    guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
    
    stat_summary(s[[i]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
    stat_summary(s[[i]], mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",size = 0.5)+
    stat_summary(s[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
    stat_summary(s[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",size = 0.5) #blackgreen
  
  print(p3)
  #ggsave("gaur_bTB_100sim.png",p,width = 25, height = 15, units = 'cm', dpi = 600)
}

## m4 - HS plot ###### 
else if (s[[i]]$model=="HS") {
  p4<-ggplot() + 
    geom_line(data = s[[i]],aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
    geom_line(data = s[[i]],aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
    geom_line(data = s[[i]],aes(x = time_y, y = R, group = run, color = 'R' ),size = 0.1, alpha = 0.12)+
    geom_line(data = s[[i]], aes(x = time_y, y = N,group = run, color = 'total'), size = 0.1, alpha = 0.12)+
    
    labs(x="years", y= "population",
         title = paste0('Gaur population with HS infection, 100 simulations')) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "population",
                        labels = c('S','I','R','total' ),
                        values = c('S'='seagreen4',
                                   'I'='firebrick',
                                   "R"='dodgerblue3',
                                   "total"='#153030'))+ #blackgreen
    theme_bw() +
    theme( plot.title = element_text(size = 18),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 15),
           legend.title=element_text(size=11),
           legend.text = element_text(size = 11),
           axis.text=element_text(size=13))+
    guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
    
    stat_summary(m[[i]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
    stat_summary(m[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
    stat_summary(m[[i]], mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",size = 0.5)+
    stat_summary(m[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",size = 0.5) #blackgreen
  
  print(p4)
  ggsave("gaur_HS_100sim.png",p, width = 25, height = 15, units = 'cm', dpi = 600)
}

## m5 - LSD plot ######  
else if (m[[i]]$model=="LSD") {
  p5<-ggplot() + 
    geom_line(data = m[[i]],aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
    geom_line(data = m[[i]],aes(x = time_y, y = E, group = run, color = 'E' ),size = 0.1, alpha = 0.12)+
    geom_line(data = m[[i]],aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
    geom_line(data = m[[i]],aes(x = time_y, y = R, group = run, color = 'R' ),size = 0.1, alpha = 0.12)+
    geom_line(data = m[[i]], aes(x = time_y, y = N,group = run, color = 'total'), size = 0.1, alpha = 0.12)+
    
    labs(x="years", y= "population",
         title = paste0('Gaur population with LSD infection, 100 simulations')) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "population",
                        labels = c('S','E','I','R','total' ),
                        values = c('S'='seagreen4',
                                   'E'='darkorange2',
                                   'I'='firebrick',
                                   "R"='dodgerblue3',
                                   "total"='#153030'))+ #blackgreen
    theme_bw() +
    theme( plot.title = element_text(size = 18),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 15),
           legend.title=element_text(size=11),
           legend.text = element_text(size = 11),
           axis.text=element_text(size=13))+
    guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
    
    stat_summary(m[[i]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
    stat_summary(m[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
    stat_summary(m[[i]], mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",size = 0.5)+
    stat_summary(m[[i]], mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",size = 0.5)+
    stat_summary(m[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",size = 0.5) #blackgreen
  
  ggsave("gaur_LSD_100sim.png",p,width = 25, height = 15, units = 'cm', dpi = 600)
}

## m6 - FMD plot ######
else if (m[[i]]$model=="FMD") {
  p<-ggplot() + 
    geom_line(data = m[[i]], aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
    geom_line(data = m[[i]], aes(x = time_y, y = E, group = run, color = 'E' ),size = 0.1, alpha = 0.12)+
    geom_line(data = m[[i]], aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
    geom_line(data = m[[i]], aes(x = time_y, y = R, group = run, color = 'R' ),size = 0.1, alpha = 0.12)+
    geom_line(data = m[[i]], aes(x = time_y, y = M, group = run, color = 'M'),size = 0.1, alpha = 0.12)+
    geom_line(data = m[[i]], aes(x = time_y, y = N,group = run, color = 'total'), size = 0.1, alpha = 0.12)+
    
    labs(x="years", y= "population",
         title = paste0('Gaur population with FMD infection, 100 simulations')) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "population",
                        labels = c('S','E','I',"R",'M','total' ),#'total change (%)'),
                        values = c('S'='seagreen4',
                                   'E'='darkorange2',
                                   'I'='firebrick',
                                   "R"='dodgerblue3',
                                   "M"='mediumorchid4',
                                   "total"='#153030'))+ 
    theme_bw() +
    theme( plot.title = element_text(size = 18),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 15),
           legend.title=element_text(size=11),
           legend.text = element_text(size = 11),
           axis.text=element_text(size=13))+
    guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
    
    stat_summary(m[[i]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
    stat_summary(m[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
    stat_summary(m[[i]], mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",size = 0.5)+
    stat_summary(m[[i]], mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",size = 0.5)+
    stat_summary(m[[i]], mapping = aes( x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",size = 0.5)+
    stat_summary(m[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",size = 0.5) #blackgreen
  
  ggsave(paste0("gaur_fmd_1sim.png"),p, width = 25, height = 15, units = 'cm', dpi = 600)
  
}

## m7 - Brucellosis plot ######
else  {
  p<-ggplot() + 
    geom_line(data = m[[i]], aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
    geom_line(data = m[[i]], aes(x = time_y, y = E, group = run, color = 'E' ),size = 0.1, alpha = 0.12)+
    geom_line(data = m[[i]], aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
    geom_line(data = m[[i]], aes(x = time_y, y = R, group = run, color = 'R' ),size = 0.1, alpha = 0.12)+
    geom_line(data = m[[i]], aes(x = time_y, y = M, group = run, color = 'M'),size = 0.1, alpha = 0.12)+
    geom_line(data = m[[i]], aes(x = time_y, y = N,group = run, color = 'total'), size = 0.1, alpha = 0.12)+
    
    labs(x="years", y= "population",
         title = paste0('Gaur population with brucellosis infection, 100 simulations')) +
    scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    scale_color_manual( name = "population",
                        labels = c('S','E','I',"R",'M','total' ),#'total change (%)'),
                        values = c('S'='seagreen4',
                                   'E'='darkorange2',
                                   'I'='firebrick',
                                   "R"='dodgerblue3',
                                   "M"='mediumorchid4',
                                   "total"='#153030'))+ 
    theme_bw() +
    theme( plot.title = element_text(size = 18),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 15),
           legend.title=element_text(size=11),
           legend.text = element_text(size = 11),
           axis.text=element_text(size=13))+
    guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
    
  ggsave(paste0("gaur_brucellosis_1sim.png"),p, width = 25, height = 15, units = 'cm', dpi = 600)
}



# multiple run plots ###### 
# ########this one can be skipped
# save the data frame  (.rds) for working next time
# Load the .rds file back and put it as a list ########
l <- list.files(path = getwd(), pattern = "_100runs.rds")
l
lm = lapply(l, readRDS)

str(lm)

#summerize basic stat
res_melt<-list()
for (i in 1:length(lm)){
  
  res_melt[[i]] <- melt(lm[[i]], id.vars = c('time_y','time_d','run','model',"Ndiff"),
                        value.name = 'value', variable.name = 'class')
}

for (i in 1:length(res_melt)){
  print(res_melt[[i]] |> group_by(class) |>
          summarise(Median = median(value), 
                    Mean = mean(value),
                    Max  =max(value),
                    Min=min(value)))
}

lm6<-lm[[7]]%>% 
  relocate(M,.after = R) %>%
  relocate(N,.after = M) 
lm<-list(lm1,lm2,lm3,lm4,lm5,lm6,lm7)

for (i in 1:length(lm)){
  names(lm)[i]
  print(head(lm[[i]]))
}

res_melt[[i]] |> group_by(class) |>
  summarise(Median = median(value), 
            Mean = mean(value),
            Max=max(value),
            Min=min(value))

saveRDS(lm1, file = "df_m1_noinf_100runs.rds")

head(m1)

pl1<-ggplot(data = m[[1]]) + 
  geom_line(aes(x = time_y ,y = a,  group = run, color = 'adult'), size = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = sa, group = run, color = 'subadult'),size = 0.1, alpha = 0.12) +
  geom_line(aes(x = time_y, y = c,  group = run, color = 'calf' ),size = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = N,  group = run, color = 'total'),size = 0.1, alpha = 0.12) +
  
  labs(x="years", y= "population",
       title= 'H) No infection') +
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
  guides(color = guide_legend(override.aes = list(alpha = 1, size=1)))
  stat_summary(m[[1]], mapping =aes( x = time_y, y = a, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
  stat_summary(m[[1]], mapping =aes( x = time_y, y = sa, group = 1), fun=mean, geom="line", colour= 'firebrick',size = 0.5)+
  stat_summary(m[[1]], mapping =aes( x = time_y, y = c, group = 1), fun=mean, geom="line", colour= 'dodgerblue3',size = 0.5)+
  stat_summary(m[[1]], mapping =aes( x = time_y, y = N, group = 1), fun=mean, geom="line", colour= '#153030',size = 0.5)

print(pl1)

ggsave("gaur_noinf_100sim.png",pl1,width = 25, height = 15, units = 'cm', dpi = 600)

pl2<-ggplot(data = lm2)+ 
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'), size = 0.1, alpha = 0.12)+
  
  labs(x="years", y= "population",
       title= 'I) Anthrax infection') +
  
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                      labels = c('S','I', 'total') , #this one label is manual
                      values = c('S'='seagreen4',
                                 'I'='firebrick',
                                 'N'='#153030'))+
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
  stat_summary(lm2, mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
  stat_summary(lm2, mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
  stat_summary(lm2, mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",size = 0.5) #blackgreen

print(pl2)



p3 <-ggplot(m3) + 
  geom_line(aes(x = time_y, y = value, color = class))+
  
  labs(x="years", y= "population",
       title= 'C) bTB infection') +
  # ylim(0,400)+
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                      labels = c('S','E','I','total' ),
                      values = c('S'='seagreen4',
                                 'E'='darkorange2',
                                 'I'='firebrick',
                                 "N"='#153030'))+ 
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))

m5 <- m[[5]] #%>% filter(class %in% c("S","E","I", "R","N"))
pl5<- ggplot(m5) + 
  geom_line(aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
  geom_line(aes(x = time_y, y = E, group = run, color = 'E' ),size = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = R, group = run, color = 'R' ),size = 0.1, alpha = 0.12)+
  geom_line(aes(x = time_y, y = N, group = run, color = 'total'), size = 0.1, alpha = 0.12)+
  
  labs(x="years", y= "population",
       title= 'L) LSD infection') +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  scale_color_manual( name = "population",
                      labels = c('S','E','I','R','total' ),#'total change (%)'),
                      values = c('S'='seagreen4',
                                 'E'='darkorange2',
                                 'I'='firebrick',
                                 "R"='dodgerblue3',
                                 "total"='#153030'))+ #blackgreen
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
  
  stat_summary(m5, mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
  stat_summary(m5, mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",size = 0.5)+
  stat_summary(m5, mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
  stat_summary(m5, mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",size = 0.5)+
  stat_summary(m5, mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",size = 0.5)
print(pl5)

head(mx[[5]])
pl5<- ggplot(mx[[5]]) + 
  geom_line(aes(x = time_y ,y = value, group = run, color = class),size = 0.1, alpha = 0.12) +
  labs(x="years", y= "population",
       title= 'E) LSD infection') +
  scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
  #  ylim(0, 1000) +
  theme_bw() +
  theme( plot.title = element_text(size = 13),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title=element_text(size=11),
         legend.text = element_text(size = 11),
         axis.text=element_text(size=11))+
  guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
  
  stat_summary(mx[[5]], mapping =aes( x = time_y, y = value, group = 1), fun=mean, geom="line", colour = class,size = 0.5)+
  scale_color_manual( name = "population",
                      labels = c('S','E','I',"R",'total' ),
                      values = c('S'='seagreen4',
                                 'E'='darkorange2',
                                 'I'='firebrick',
                                 "R"='dodgerblue3',
                                 "N"='#153030'))+ 
  print(pl5)

scale_color_manual