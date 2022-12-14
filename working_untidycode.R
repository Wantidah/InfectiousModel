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

nam<-c('pop_dynamic',
       'Anthrax',
       'bTB',
       'HS',
       'LSD',
       'FMD',
       'Brucellosis')

sim_rep_m<-list(sim_rep_m7)

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

# ########this one can be skipped
# save the data frame  (.rds) for working next time
for (i in 1:length(m)) {
saveRDS(m[[i]], file = paste0("df_",nam[[i]],".rds")) }

# Load the .rds file back and put it as a list ########
l <- list.files(path = getwd(), pattern = "*.rds")
m = lapply(l, readRDS)
m

#######################################
# OR if we want to have melt data for plotting the graph (probably is easier to melt the data earlier)
# loop for melt 

m2<-list()

for (i in 1:length(m)){

m2[[i]] <- melt(m[[i]], id.vars = c('time_d','time_y','run','model','Ndiff'),
           value.name = 'value', variable.name = 'class')
}
str(m2)

# PLOT MODEL OUTPUTS: population dynamic line graphs #############################
for (i in 1:length(m2)){
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

m3 <-data.table::rbindlist(m2)
str(m3)
write_csv(m3,"df_melt.csv")
#Prepare the dataframe for boxplotting ########
#the % of population change in 100 for all the models
dft<-m3 %>% 
  dplyr::select(Ndiff,run,model)%>% 
  drop_na()%>%
  distinct()




dft$model <- recode_factor(dft$model, pop_dynamic = "no_infection" )
table(dft$model)
str(dft)

write_csv(dft,"df_all.csv")

#load df back
dft <- read_csv("df_all.csv")

#calculating mean,median
s <- dft |>
  group_by(model)|>
  summarise(Median = median(Ndiff),
            Mean = mean(Ndiff))

s

# Boxplot
box<-dft %>%
  ggplot(aes(x= reorder(model,-Ndiff), y = Ndiff,fill=model)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.6, alpha=0.4) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)) +
  ggtitle("Gaur population change by disease models in 100 years") +
  xlab("") +
  ylab("Population change (%)")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                                       breaks = seq(from = -200, to = 500, by = 100))


print(box)
ggsave("gaur_ndiff_box1_100sim.png",box,width = 20, height = 15, units = 'cm', dpi = 600)

#violoin plot
v<-dft %>%
  ggplot(aes(x=reorder(model,-Ndiff),y = Ndiff,fill=model)) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Gaur population change by disease models in 100 years") +
  xlab("")+
  ylab("Population change (%)")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))


print(v)

ggsave("gaur_ndiff_vplot_100sim.png",v,width = 25, height = 15, units = 'cm', dpi = 600)

#reorder by population change max-min
dft$model <- factor(dft$model, 
                     levels = c("no_infection", 
                                "HS",
                                "LSD",
                                "Anthrax",
                                "FMD",
                                "bTB",
                                "Brucellosis"))
#?ggbetweenstats()

plt<-dft%>%ggbetweenstats(
  x=model,
  y=Ndiff,
  k=0,
  plot.type = "boxviolin",
  pairwise.comparisons=F,
  bf.message = F,
  results.subtitle = FALSE,
  centrality.point.args = list(size = 2, color = "darkred"),
  title= "Gaur population change by disease models in 100 years",
  xlab = "",
  ylab = "Population change (%)",
  package = "ggsci",
  palette = "default_jco")+
  ggplot2::scale_y_continuous(limits = c(-200, 500), 
                              breaks = seq(from = -200, to = 500, by = 100))
  
print(plt)

ggsave("gaur_ndiff_boxviolin_max-min2_100sim.png",plt,width = 20, height = 15, units = 'cm', dpi = 600)

## plot extinction times


df.agg <- aggregate(time_y ~ run + value + model, m3, min)
df.ag<-(df.agg[df.agg$value==0,c('model','time_y')])

neworder <- c("1","2","3","4","5","6","7","8")

df.ag <- arrange(transform(df.ag,
                           model=factor(model,levels=neworder)),model)
labs <- c('1' = "no_infectious",
          '2' = "SI (Anthrax)",
          '3' = "",
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
