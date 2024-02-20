library(PCAtools)
library(factoextra)
library(FactoMineR)
library(tidyverse)
library(readxl)
library(wesanderson)

#import PCA excel file
df_pca <- read_xlsx("./pca_table.xlsx", sheet = "norm") %>% 
  as.data.frame()

df_pca2<-df_pca %>% dplyr::select(c('beta','incubation',"infectious","fatality"))

rownames(df_pca2)<-df_pca$model

mt<-data.matrix(df_pca2)
head(mt)

gpca <- FactoMineR::PCA(mt, graph = FALSE,scale.unit = T)

#eigenvalues
get_eigenvalue(gpca)
head(gpca$var$contrib)
head(gpca$var)
cat <- length(unique(df_pca$disease))
pal <- wesanderson::wes_palette("Cavalcanti1", 
                                length(unique(df_pca$disease)), 
                                type = "continuous")[1:length(unique(df_pca$disease))]

#pal2<- (RColorBrewer::brewer.pal(cat,"Spectral"))

pal3 <- wes_palette("Zissou1", 6, type = "continuous")

RColorBrewer::display.brewer.pal(cat,'Spectral')

#correlation plot
library(corrplot)
corrplot(gpca$var$cos2, is.corr=F,col=pal) #component contribution

# Total contribution on PC1 and PC2
# fviz_contrib(gpca, choice = "ind", axes = 1:4)

fviz_pca_var(gpca, col.var = "black")
# Contributions of variables to PC1
fviz_contrib(gpca, choice = "var", axes = 1, top = 10)

# Contributions of variables to PC2
fviz_contrib(gpca, choice = "var", axes = 2, top = 10)
# Disease
fviz_pca_biplot(gpca, label ="var", col.ind="cos2",
                habillage=factor(df_pca$disease)) + 
  theme_minimal()

#PCATools
head(df_pca2)
df_t2<-data.table::transpose(df_pca2)

rownames(df_t2) <- colnames(df_pca2)
colnames(df_t2) <- rownames(df_pca2)
head(df_t2)
colnames(df_t2)

metadata<-df_pca
rownames(metadata)
rownames(metadata) <- df_pca$model
head(metadata)

pc<-pca(df_t2,metadata=metadata,
        scale=T,
        center=T)
pc

getComponents(pc)

#pairsplot(pc,
#          components = getComponents(pc, c(1:4)),
#          triangle = F, trianglelabSize = 12,
#          hline = 0, vline = 0,
#          pointSize = 0.8,
#          gridlines.major = FALSE, 
#          gridlines.minor = FALSE,
#          colby = 'Nchange_p',
#          title = 'Pairs plot', 
#          plotaxes = T,
#          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))

limits<-c(-200,-100,0,100,200)

bi_pc<-biplot(pc,
              showLoadings = TRUE,
              lab = pc$metadata$model,
              colby = 'Nchange_p',
              hline = 0, vline = 0,
              shape = 'disease',
              sizeLoadingsNames = 5,
              boxedLoadingsNames = F,
              colLoadingsArrows = 'grey30',
              legendPosition = 'right',
              legendLabSize = 11, legendIconSize = 5,
              xlim = c(-10,10),
              ylim = c(-5,5),
              title = "PCA biplot",
              subtitle="Diseases parameters contribute to the % of the population change",
              titleLabSize = 16,
              subtitleLabSize = 15,
              max.overlaps = 7 #ggrepel
)+
  scale_color_gradientn(colours = rev(pal3))+
  guides(color = guide_colorbar(limits = limits))+
  labs(color = "population(%)")
bi_pc
