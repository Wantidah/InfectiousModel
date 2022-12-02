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

t<-function(x=time,y=N,...){
  res_no<-vector()
  for (i in 1:length(x)){
    if (i == 1) {
      res_no = 0
    }
    else { res_no = (y[i+1]-y[i])/100 }
  }
}

# Dave's adjust version
single_pop_sim_prep <- function(x, n_rep, end.time, melt){ # x = simulation of model, e.g. sim_run_m1
  df<-list()
  #loop for storing new df
  for (i in 1:n_rep){
    run <- paste("run", seq(n_rep), sep="")
    names(df)[i]
    df[[i]]<- x[,i]$results[,-c(1)]
    df[[i]]$time <- seq(from=1,to=end.time+1,by=1)
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

wow<-single_pop_sim_prep(x = m6, n_rep=n_rep, end.time= end.time, melt = T)
wow2<- single_pop_sim_prep(x = m6, n_rep=n_rep, end.time= end.time, melt = F)
#total N
dave<-single_pop_sim_prep(x = m6, n_rep=n_rep, end.time= end.time)

table(dave$time)
table(dave$run)
table(wow2$time)
table(wow2$run)
str(dave)
str(wow2)

#plot
#dave's
#for (i in unique(res_mx_p$model)){
#  subdata <- subset(res_mx_p, model == i)
#  pdf(paste("plot_ts_mx", i, ".pdf", sep = ""), width = 4, height = 3)
  print(ggplot(wow2, aes(x=time, y=N, group=run)) +
          theme_bw() +
          theme(panel.grid=element_blank()) +
          geom_line(size=0.2, alpha=0.15)+
          ylab('Numbers') + xlab('time')+
          stat_summary(aes(group = 1), fun=mean, geom="line", colour="black",size = 1.1))

#}