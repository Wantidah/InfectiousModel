# try to make dataframe with all the population class
# ongoning
# the for loop is works and can make a dataframe we want, but kust have to find a nicer way of plotting graph

n_rep = 20
end.time = 3650

m6<-replicate(n_rep,(model6(pars = parameters_m6, init = initials_m6,
                                    end.time = end.time)))

#TEST: create data.frame from model with population class == works :)
df<-list()

for (i in 1:n_rep){
  run <- paste("run", seq(n_rep), sep="")
  names(df)[i]
  df[[i]]<- m6[,i]$results[,-c(1)]
  df[[i]]$time <- seq(from=1,to=end.time+1,by=1)
}
  df2<-map2(df,run, ~cbind(.x, run = .y))
  df3<- data.table::rbindlist(df2) #bind row

str(df3)
table(df3$N)

#TEST 2: create function calculating N change from time t to t+1
n_change<-function(x=time,y=N,...){
  res_no<-vector()
  
}
my_imp_ext_na_I<-function(x=time,y=I,...){
  
  res_no<-vector()
  
  for (i in 2:length(x)) {
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

t<-function(x=x, y=y,...){
  res_no<-vector()
  for (i in 1:length(x)){
    if (x == 1) {
      res_no == 0
    }
    else { res_no = ((y[i+1]-y[i])/y[i+1])*100 }
    
  }
  res_no
}
1:length(w$time)
w<-wow2

w2<-w %>%
  group_by(run) %>%
  mutate(N_change = ((N - lag(N))/N)*100,
         year = time/365)


table(is.na(w2$N_change))
table((w2$N_change==0))
mean(w2$N_change)
?round()
View(w2)

w$N_change<- t(x=wow2$time,y=wow2$N)
str(w)

# single_pop_sim_prep Dave's adjust version
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
     
      df3 <- gather(df2, key = class, value = value, -c(time,run))
     
    return(df3 = data.frame(df3))
     } 
   
     else  {
      return( df2 = data.frame(df2)) 
     }

   }

t<-df_pop %>% gather(key = class, value = value, -c(time_d,time_y,run))
str(t)
head(t)

#dave's
single_pop_sim_prep <- function(x, n_rep, end.time){ # x = simulation of model, e.g. sim_run_m1
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

sim_rep_m<-list(sim_rep_m1,
                sim_rep_m2,
                sim_rep_m3,
                sim_rep_m4,
                sim_rep_m5,
                sim_rep_m6,
                sim_rep_m7)

m<-list()

nam<-c('pop_dynamic',
       'Anthrax',
       'bTB',
       'HS',
       'LSD',
       'FMD',
       'Brucellosis')
nam<-c('Anthrax')
#group and calculate total population change (%) loop--------
for (i in 1:length(sim_rep_m)) {
  m[[i]]<- single_pop_sim_prep(x = sim_rep_m[[i]], n_rep=n_rep, end.time= end.time, melt = F)
  m[[i]]<- m[[i]]%>%
    group_by(run) %>%
    mutate(Ndiff = ((last(N)-first(N))/first(N))*100)%>% #calculate change in the total population at year100, and year0
    mutate(time_y = time_d/365) %>% #convert day to year for plotting
    as.data.frame()
  
    m[[i]]$model <- paste0(nam[[i]])
   
}

#in case we want to save the data frame  (.rds)
for (i in 1:length(m)) {
saveRDS(m[[i]], file = paste0("df_",nam[[i]],".rds")) }
#########this one can skip##############
#in case load the .rds file
#add into the list()
m<-list(df_pop_dynamic,
        df_Anthrax,
        df_bTB,
        df_HS,
        df_LSD,
        df_fmd,
        df_Brucellosis)

#######################################

# PLOT MODEL OUTPUTS: population dynamic line graphs #############################

for (i in 1:length(m)){
  #subdata <- subset(m, model == k)
  
  ## m1 - no infection plot ###### 
  if (subset(m[[i]], model=="pop_dynamic")==T){
  #if (m[[i]]$model=="pop_dynamic"){
     p<-ggplot() + 
      geom_line(data = m[[i]],aes(x = time_y ,y = a,  group = run, color = 'adult'),size = 0.1, alpha = 0.12) + 
      geom_line(data = m[[i]],aes(x = time_y, y = sa, group = run, color = 'subadult' ),size = 0.1, alpha = 0.12) +
      geom_line(data = m[[i]],aes(x = time_y, y = c,  group = run, color = 'calf' ),size = 0.1, alpha = 0.12)+
      geom_line(data = m[[i]], aes(x = time_y, y = N, group = run, color = 'total'),size = 0.1, alpha = 0.12) +
      
      labs(x="years", y= "population",
           title= 'Gaur population - no infection, 100 simulations') +
      scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
      scale_color_manual( name = "population",
                          labels = c('adult','subadult','calf','total' ),
                          values = c('adult'='seagreen4',
                                     'subadult'='firebrick',
                                     "calf"='dodgerblue3',
                                     "total"='#153030'))+ 
      theme_bw() +
      theme( plot.title = element_text(size = 18),
             axis.title.x = element_text(size = 15),
             axis.title.y = element_text(size = 15),
             legend.title=element_text(size=11),
             legend.text = element_text(size = 11),
             axis.text=element_text(size=13))+
      guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
      
      stat_summary(m[[i]], mapping =aes( x = time_y, y = a, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
      stat_summary(m[[i]], mapping = aes( x = time_y, y = sa, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
      stat_summary(m[[i]], mapping = aes( x = time_y, y = c, group = 1), fun=mean, geom="line", colour="#153030",size = 0.5)+
      stat_summary(m[[i]], mapping = aes( x = time_y, y = N, group = 1), fun=mean, geom="line", colour="dodgerblue3",size = 0.5) #blackgreen
    
      ggsave("gaur_no_infection_100sim.png",p, width = 25, height = 15, units = 'cm', dpi = 600)
  }
  
  ## m2 - Anthrax plot  ###### 
  else if (m[[i]]$model=="Anthrax"){
    p<-ggplot() + 
      geom_line(data = m[[i]],aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
      geom_line(data = m[[i]],aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
      geom_line(data = m[[i]], aes(x = time_y, y = N,group = run, color = 'total'), size = 0.1, alpha = 0.12)+
      
      labs(x="years", y= "population",
           title= 'Gaur population with anthrax infection, 1 simulation, 100 years') +
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
      guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
      
      stat_summary(m[[i]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
      stat_summary(m[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
      stat_summary(m[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",size = 0.5) #blackgreen
    
      ggsave("gaur_anthrax_100sim.png",p,width = 25, height = 15, units = 'cm', dpi = 600)
      }
  
  ## m3 - bTB plot ######  
  else if (m[[i]]$model=="bTB") {
    p<-ggplot() + 
      geom_line(data = m[[i]],aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
      geom_line(data = m[[i]],aes(x = time_y, y = E, group = run, color = 'E' ),size = 0.1, alpha = 0.12)+
      geom_line(data = m[[i]],aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
      geom_line(data = m[[i]], aes(x = time_y, y = N,group = run, color = 'total'), size = 0.1, alpha = 0.12)+
      
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
             
             stat_summary(m[[i]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
             stat_summary(m[[i]], mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",size = 0.5)+
             stat_summary(m[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
             stat_summary(m[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",size = 0.5) #blackgreen
           
      ggsave("gaur_bTB_100sim.png",p,width = 25, height = 15, units = 'cm', dpi = 600)
  }
  
  ## m4 - HS plot ###### 
  else if (m[[i]]$model=="HS") {
    p<-ggplot() + 
      geom_line(data = m[[i]],aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
      geom_line(data = m[[i]],aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
      geom_line(data = m[[i]],aes(x = time_y, y = R, group = run, color = 'R' ),size = 0.1, alpha = 0.12)+
      geom_line(data = m[[i]], aes(x = time_y, y = N,group = run, color = 'total'), size = 0.1, alpha = 0.12)+
      
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
    
    ggsave("gaur_HS_100sim.png",p, width = 25, height = 15, units = 'cm', dpi = 600)
  }
  
 ## m5 - LSD plot ######  
 else if (m[[i]]$model=="LSD") {
  p<-ggplot() + 
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
  
  ggsave(paste0("gaur_fmd_100sim.png"),p, width = 25, height = 15, units = 'cm', dpi = 600)
  
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
      
      stat_summary(m[[i]], mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
      stat_summary(m[[i]], mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
      stat_summary(m[[i]], mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",size = 0.5)+
      stat_summary(m[[i]], mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",size = 0.5)+
      stat_summary(m[[i]], mapping = aes( x = time_y, y = M, group = 1), fun=mean, geom="line", colour="mediumorchid4",size = 0.5)+
      stat_summary(m[[i]], mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",size = 0.5) #blackgreen
    
    ggsave(paste0("gaur_brucellosis_100sim.png"),p, width = 25, height = 15, units = 'cm', dpi = 600)
  }
  
 }
#load .rds file
m<-list(df_pop_dynamic,
        df_Anthrax,
        df_bTB,
        df_HS,
        df_LSD,
        df_fmd,
        df_Brucellosis)
names(m[[1]])


mt<-list()
for(i in 1:length(m))  {
  m[[i]]<-m[[i]] %>% 
    dplyr::select(Ndiff,run,model)%>%
    drop_na()%>%
    distinct()
  mt<-data.table::rbindlist(mt)
}

head(dft2)
tail(dft2)
dft2<-data.table::rbindlist(mt)
str(dft2)
table(is.na(dft2))
dft2<-drop_na(dft2)

write_csv(dft2,"df_all.csv")

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggstatsplot)

# Boxplot
box<-dft2 %>%
  ggplot( aes(x=model, y=Ndiff, fill=model)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Gaur population change (%) with and without infectious diseases") +
  xlab("Model")

print(box)
ggsave("gaur_ndiff_100sim.png",box,width = 25, height = 15, units = 'cm', dpi = 600)

v<-dft2 %>%
  ggplot( aes(x=model, y=Ndiff, fill=model)) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Gaur populatin change (%)") +
  xlab("Model")

print(v)

ggsave("gaur_ndiff_vplot_100sim.png",v,width = 25, height = 15, units = 'cm', dpi = 600)

plt <- ggbetweenstats(
  data = dft2,
  x = model,
  y = Ndiff
)

## plot extinction times
head(mt)
str(m)
df.agg <- aggregate(time ~ run + value + model, res_mx_p, min)
df.ag<-(df.agg[df.agg$value==0,c('model','time')])

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
