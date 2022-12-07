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

df_m3<-single_pop_sim_prep(x = sim_rep_m3, n_rep=n_rep, end.time= end.time, melt = T)
df_m3<- single_pop_sim_prep(x = sim_rep_m3, n_rep=n_rep, end.time= end.time, melt = F)
df_m3_d<-single_pop_sim_prep(x = sim_rep_m3, n_rep=n_rep, end.time= end.time)

df_m2<-single_pop_sim_prep(x = sim_rep_m2, n_rep=n_rep, end.time= end.time, melt = F)
df_m4<-single_pop_sim_prep(x=sim_rep_m4, n_rep=n_rep, end.time= end.time, melt = F)

df_m6<-single_pop_sim_prep(x=sim_rep_m6, n_rep=n_rep, end.time= end.time, melt = F)

df_m6<- df_m6%>%
  group_by(run) %>%
  mutate(N_change = ((N - lag(N))/lag(N))*100) %>% #calculate change percentages in the total population
  mutate(time_y = time_d/365) %>% #convert day to year for plotting
  mutate(Ndiff = ((last(N)-first(N))/last(N))*100) %>% #calculate % population change between Nt-N0; NOTE: Ndiff will give the same number for each run group
  as.data.frame()


View(df_m6)
str(df_m6)

#mutate(across(everything(), ~ if_else(row_number() %in% 3:5, ./2, .)))

#plot

png("gaur_fmd_100y_all.png",width = 25, height = 15, units = 'cm', res = 600)
ggplot() + 
  geom_line(data = df_m6,aes(x = time_y ,y = S, group = run, color = 'S'),size = 0.1, alpha = 0.12) + 
  geom_line(data = df_m6,aes(x = time_y, y = E, group = run, color = 'E' ),size = 0.1, alpha = 0.12)+
  geom_line(data = df_m6,aes(x = time_y, y = I, group = run, color = 'I' ),size = 0.1, alpha = 0.12)+
  geom_line(data = df_m6,aes(x = time_y, y = R, group = run, color = 'R' ),size = 0.1, alpha = 0.12)+
  geom_line(data = df_m6, aes(x = time_y, y = M,group = run, color = 'M'), size = 0.1, alpha = 0.12)+
  geom_line(data = df_m6, aes(x = time_y, y = N,group = run, color = 'total'), size = 0.1, alpha = 0.12)+
  #geom_line(data = df_m6, aes(x = time_y, y = N_change, group = run, color = 'total change (%)' ),size = 0.1, alpha = 0.15)+
  
  labs(x="years", y= "population",
       title= 'Gaur population with FMD infection, 100 simulations') +
    
   scale_x_continuous(breaks=seq(0, (365*100), by = 10))+
    
   scale_color_manual( name = "population",
                       labels = c('S',
                                  'E',
                                  'I',"R",
                                  'M',
                                  'total' ),#'total change (%)'),
                       values = c('S'='seagreen4',
                                  'E'='darkorange2',
                                  'I'='firebrick',
                                  "R"='dodgerblue3',
                                  "M"='lavender',
                                  "total"='#153030'))+ #blackgreen
                                  #'total change (%)'='#0D9EAD' #teal
                                  
      theme_bw() +
      theme( plot.title = element_text(size = 18),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 15),
           legend.title=element_text(size=11),
           legend.text = element_text(size = 11),
           axis.text=element_text(size=13))+
      guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
  
  stat_summary(df_m6, mapping =aes( x = time_y, y = S, group = 1), fun=mean, geom="line", colour='seagreen4',size = 0.5)+
  stat_summary(df_m6, mapping = aes( x = time_y, y = E, group = 1), fun=mean, geom="line", colour="darkorange2",size = 0.5)+
  stat_summary(df_m6, mapping = aes( x = time_y, y = I, group = 1), fun=mean, geom="line", colour="firebrick",size = 0.5)+
  stat_summary(df_m6, mapping = aes( x = time_y, y = R, group = 1), fun=mean, geom="line", colour="dodgerblue3",size = 0.5)+
  stat_summary(df_m6, mapping = aes( x = time_y, y = M, group = 1), fun=mean, geom="line", colour="lavender",size = 0.5)+
  stat_summary(df_m6, mapping = aes(x = time_y, y = N, group = 1), fun=mean,geom="line", colour="#153030",size = 0.5) #blackgreen
  #stat_summary(df_m6, mapping = aes(x = time_y, y = N_change, group = 1), fun=mean, geom="line", colour="#0D9EAD",size = 0.5) #teal

  dev.off()
