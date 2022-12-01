# try to make dataframe with all the population class
# ongoning
# the for loop is works and can make a dataframe we want, but kust have to find a nicer way of plotting graph

n_rep = 2
end.time = 730

m2<-replicate(n_rep,(model2(pars = parameters, init = initials,
                                    end.time = end.time)))
m2

df<-list()

for (i in 1:n_rep){
  run <- paste("run", seq(n_rep), sep="")
  names(df)[i]
  df[[i]]<- m2[,i]$results[,-c(1)]
  df[[i]]$time <- seq(from=1,to=end.time+1,by=1)
  
  df2<-map2(df, run, ~cbind(.x, run = .y)) %>%
  
}



str(s2)
single_pop_sim_prep <- function(x, n_rep, end.time){ # x = simulation of model, e.g. sim_run_m1
  df<-list()
  
  for (i in 1:n_rep){
    run <- paste("run", seq(n_rep), sep="")
    names(df)[i]
    df[[i]]<- m2[,i]$results[,-c(1)]
    df[[i]]$time <- seq(from=1,to=end.time+1,by=1)
    df2<-map2(df, run, ~cbind(.x, run = .y))
    
  }
  
}

s2 = melt(s, measure.vars = s[,-c("time")] )
s2<- s %>% gather(key = class, value = value, -c(time,run))
str(s2)
View(s2)
dave<-single_pop_sim_prep(x=m2,n_rep = n_rep, end.time = end.time)
str(dave)

print(ggplot(s, aes(x=time, y=value, group = run))) +
        theme_bw() +
        theme(panel.grid=element_blank()) +
        geom_line(aes(col=class), size=0.2, alpha=0.15)+
        ylab('Numbers') + xlab('time')
        stat_summary(aes(group = 1), fun=mean, geom="line", colour="black",size = 1.1))


print(ggplot(s, aes(x=time, y=S, group= run))) +
        theme_bw() +
        theme(panel.grid=element_blank()) +
        geom_line(size=0.2, alpha=0.15)+
        ylab('I (Numbers)') + xlab('time')+
        #ggtitle("bTB infectious of gaur population, 100 runs")+
        theme(plot.title = element_text(size=8))+
        stat_summary(aes(group = 1), fun=mean, geom="line", colour="black",size = 1.1)
