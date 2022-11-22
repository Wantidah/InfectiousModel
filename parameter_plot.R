#plot diseases parameters

library(readxl)
p <- read_excel("Desktop/parameters.xlsx")
View(p) 

p$diseases<-as.factor(p$diseases)
p$sp<-as.factor(p$sp)
str(p)

p$diseases

png(filename = paste0("/Users/whorpien/Library/CloudStorage/OneDrive-MasseyUniversity/InfectiousDisease_RCode/figures/all_diseases.png"),
    res=600, units = "cm", width= 25, height =20)

plot_histogram(p,title = "diseases parameters")

dev.off()

fmd<- p %>% filter(diseases=="FMD")
hs<- p %>% filter(diseases=="HS")
ba<- p %>% filter(diseases=="Brucella abortus")
bm<- p %>% filter(diseases=="Brucella melitensis")
btb<- p %>% filter(diseases=="BTB")
ath<- p %>% filter(diseases=="Antrax")
lsd<- p %>% filter(diseases=="LSD")

png(filename = paste0("/Users/whorpien/Library/CloudStorage/OneDrive-MasseyUniversity/InfectiousDisease_RCode/figures/anthrax.png"),
    res=600, units = "cm", width= 25, height =20)
plot_histogram(ath,title = "Anthrax")
dev.off()

png(filename = paste0("/Users/whorpien/Library/CloudStorage/OneDrive-MasseyUniversity/InfectiousDisease_RCode/figures/fmd.png"),
    res=600, units = "cm", width= 25, height =20)
plot_histogram(fmd,title = "FMD")
dev.off()

png(filename = paste0("/Users/whorpien/Library/CloudStorage/OneDrive-MasseyUniversity/InfectiousDisease_RCode/figures/hs.png"),
    res=600, units = "cm", width= 25, height =20)
plot_histogram(hs,title = "HS")
dev.off()

png(filename = paste0("/Users/whorpien/Library/CloudStorage/OneDrive-MasseyUniversity/InfectiousDisease_RCode/figures/b_abortus.png"),
    res=600, units = "cm", width= 25, height =20)
plot_histogram(ba,title = "B. abortus")
dev.off()

png(filename = paste0("/Users/whorpien/Library/CloudStorage/OneDrive-MasseyUniversity/InfectiousDisease_RCode/figures/b_melitensis.png"),
    res=600, units = "cm", width= 25, height =20)
plot_histogram(bm,title = "B. melitensis")
dev.off()

png(filename = paste0("/Users/whorpien/Library/CloudStorage/OneDrive-MasseyUniversity/InfectiousDisease_RCode/figures/btb.png"),
    res=600, units = "cm", width= 25, height =20)
plot_histogram(btb,title = "Bovine TB")
dev.off()

png(filename = paste0("/Users/whorpien/Library/CloudStorage/OneDrive-MasseyUniversity/InfectiousDisease_RCode/figures/lsd.png"),
    res=600, units = "cm", width= 25, height =20)
plot_histogram(lsd,title = "LSD")
dev.off()

#list working

d<-c("Antrax", "Brucella abortus", "Brucella melitensis", "BTB", "FMD", "HS", "LSD")
d

ds<-list(ath,ba,bm,btb,fmd,hs,lsd)
ds

g<-list()

for(i in ds) {

  g[[i]]<-plot_histogram(ds[i],title = d[[i]])
png(filename = paste0("/Users/whorpien/Library/CloudStorage/OneDrive-MasseyUniversity/InfectiousDisease_RCode/figures/",d[[i]],".png"),
res=600, units = "cm", width= 20, height =20)
print(g[[i]])
dev.off()
}

clin<-ggplot(fmd,aes(infclin_d))#infectious_cli
clin+geom_histogram()+facet_wrap (~sp)
str(fmd)
ggplot(fmd, aes(x=infclin_d,fill=sp))+
geom_histogram()+facet_wrap (~sp)
