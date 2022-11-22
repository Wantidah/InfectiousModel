library(SimInf)
## Create an 'SEIR' model with 1600 nodes and initialize
## it to run over 4*365 days. Add one infected individual
## to the first node.
u0 <-u0_SIR()
u0$I[1] <- 1
tspan <- seq(from = 1, to = 365*3, by = 1)
model1 <- SIR(u0 = u0,
              tspan = tspan,
              events = events_SIR(),
              beta = 0.2,
              gamma = 0.1)

u0 <-u0_SEIR()
u0$I[1] <- 1
tspan <- seq(from = 1, to = 365*3, by = 1)
model2 <- SEIR(u0 = u0,
               tspan = tspan,
               events = events_SEIR(),
               beta = 0.2,
               epsilon = 0.11,
               gamma = 0.1)

## Run the model to generate a single stochastic trajectory.
result1 <- run(model1, threads = 1)
summary(result1)
result2 <- run(model2, threads = 1)
summary(result2)

sir_col<-c(rgb(0.1,0.1,0.1,0.2),
           'red',
           'darkgreen')
seir_col<-c(rgb(0.1,0.1,0.1,0.2),
            'blue',
            'red',
            'darkgreen')

pdf("sir_1600_structure_iqr.pdf", width = 4, height = 4)
plot(result1, col = c('black','red','darkgreen'),par(mar=c(3,4,1,1)))
dev.off()

pdf("sir_1600_structure_100.pdf", width = 4, height = 4)
plot(result1, node = 1:100, range = FALSE, col = sir_col,par(mar=c(3,4,1,1)))
dev.off()

pdf("sir_1600_structure_boxplot.pdf", width = 4, height = 4)
boxplot(result1, col = c('grey','red','darkgreen'),ylab = 'N',par(mar=c(4,5,1,1)))
dev.off()

pdf("seir_1600_structure_iqr.pdf", width = 4, height = 4)
plot(result2, col = c('black','blue','red','darkgreen'),par(mar=c(3,4,1,1)))
dev.off()

pdf("seir_1600_structure_100.pdf", width = 4, height = 4)
plot(result2, node = 1:100, range = FALSE, col = seir_col,par(mar=c(3,4,1,1)))
dev.off()

pdf("seir_1600_structure_boxplot.pdf", width = 4, height = 4)
boxplot(result2, col = c('grey','blue','red','darkgreen'),ylab = 'N',par(mar=c(4,5,1,1)))
dev.off()

pdf("spatial_structure_coupling.pdf", width = 7, height = 5)
par(mfrow=c(1,2))
plot(NULL, xlim = c(0, 365*3), ylim = c(0, 0.001),
     ylab = "Prevalance", xlab = "Time", main = 'SIR',font.main = 1, bty = "n")

couple<-seq(from = 0.1,to = 0.9, by = 0.1)
for (i in 1:length(couple)){
  gdata(model1, "coupling") <- couple[i]
  replicate(100, {
    result <- run(model = model1, threads = 1)
    p <- prevalence(model = result, formula = I ~ S + I + R, type = "pop")
    lines(p, col = rgb(1-couple[i],couple[i],couple[i],0.3), lty = 2)
  })
}
col_p<-vector()
for (i in 1:length(couple)){
  col_p[i] = rgb(1-couple[i],couple[i],couple[i],0.3)
}

plot(NULL, xlim = c(0, 365*3), ylim = c(0, 0.001),
     ylab = "", xlab = "Time", main = 'SEIR',font.main = 1, bty = "n")
legend('topright', legend=couple[1:9],col=col_p, lty=2, cex=0.8, bty = "n")
couple<-seq(from = 0.1,to = 0.9, by = 0.1)
for (i in 1:length(couple)){
  gdata(model2, "coupling") <- couple[i]
  replicate(100, {
    result <- run(model = model2, threads = 1)
    p <- prevalence(model = result, formula = I ~ S + E + I + R, type = "pop")
    lines(p, col = rgb(1-couple[i],couple[i],couple[i],0.3), lty = 2)
  })
}
dev.off()
