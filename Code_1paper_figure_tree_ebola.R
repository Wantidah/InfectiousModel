# install.packages("phytools")
# install from Bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite("EBImage")
# biocLite("ggtree")
# https://www.molecularecologist.com/2017/02/phylogenetic-trees-in-r-using-ggtree/
# source("http://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ggtree", version = "3.8")
# install.packages("remotes")
# remotes::install_github("sckott/rphylopic")

library("ape")
# library("Biostrings")
library("ggplot2")
library("ggtree")
library('rphylopic')
library(ggimage)
library("phytools") # for sims and ASRs
library("EBImage") # for images
library("colorspace")

par(mfrow=c(1,2))


set.seed(10)
tr <- rtree(20)
ggtree(tr) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

tr <- groupClade(tr, .node=c(24,28,29))
p1 <- ggtree(tr, aes(color=group, linetype=group)) + 
  scale_color_manual(values=c("black", 'red','red','blue')) #+
# geom_hilight(node=2, fill="pink", alpha=0.5, extend = 0.4) +
# geom_hilight(node=3, fill="pink", alpha=0.5, extend = 0.4) +
# geom_hilight(node=7, fill="pink", alpha=0.5, extend = 0.4) +
# geom_hilight(node=8, fill="pink", alpha=0.5, extend = 0.4)+
# geom_hilight(node=6, fill="pink", alpha=0.5, extend = 0.4) 
d <- data.frame(node = c("1","2","3",'4','5','6',
                         "7","8","9",'10','11','12',
                         '13','14','15','16','17',
                         '18','19','20'),
                images = c("bat.png",'chimp.png',
                           'chimp.png','bat.png',
                           'bat.png','chimp.png',
                           "person.png",'person.png',
                           'bat.png','bat.png',
                           'bat.png','bat.png',
                           'bat.png','bat.png',
                           'bat.png','bat.png',
                           'bat.png','bat.png',
                           'bat.png','bat.png'))
p2<-p1 %<+% d + geom_tiplab(aes(image=images), geom="image", offset = 0.1)
p2<- p2 + xlim(0, 4)+
  geom_point2(aes(subset=(node == 18)), size=5, shape=21, 
              fill="pink")+
  geom_point2(aes(subset=(node == 17)), size=5, shape=21, 
              fill="pink")
pdf('phylo1.pdf',width = 5,height = 6)
p2
dev.off()

set.seed(10)
#tr <- rtree(12)
tr <- groupClade(tr, .node=c(24,28,29,34,38))
p3 <- ggtree(tr, aes(color=group, linetype=group)) + 
  scale_color_manual(values=c("black", 'red','red','blue','red','black')) #+
# geom_hilight(node=2, fill="pink", alpha=0.5, extend = 0.4) +
# geom_hilight(node=3, fill="pink", alpha=0.5, extend = 0.4) +
# geom_hilight(node=20, fill="pink", alpha=0.5, extend = 0.4) +
# geom_hilight(node=19, fill="pink", alpha=0.5, extend = 0.4) +
# geom_hilight(node=16, fill="pink", alpha=0.5, extend = 0.4) +
# geom_hilight(node=14, fill="pink", alpha=0.5, extend = 0.4) +
# geom_hilight(node=15, fill="pink", alpha=0.5, extend = 0.4) +
# geom_hilight(node=7, fill="pink", alpha=0.5, extend = 0.4) +
# geom_hilight(node=8, fill="pink", alpha=0.5, extend = 0.4)+
# geom_hilight(node=6, fill="pink", alpha=0.5, extend = 0.4) 
d <- data.frame(node = c("1","2","3",'4','5','6',
                         "7","8","9",'10','11','12',
                         '13','14','15','16','17',
                         '18','19','20'),
                images = c("bat.png",'chimp.png',
                           'chimp.png','bat.png',
                           'bat.png','chimp.png',
                           "person.png",'person.png',
                           'bat.png','bat.png',
                           'bat.png','bat.png',
                           'bat.png','chimp.png',
                           'chimp.png','chimp.png',
                           'bat.png','bat.png',
                           'chimp.png','chimp.png'))
p4<-p3 %<+% d + geom_tiplab(aes(image=images), geom="image", offset = 0.1)
p4<- p4 + xlim(0, 4)+
  # ggtitle('B')+
  geom_point2(aes(subset=(node == 18)), size=5, shape=21, 
              fill="pink")+
  geom_point2(aes(subset=(node == 17)), size=5, shape=21, 
              fill="pink")

pdf('phylo2.pdf',width = 5,height = 6)
p4
dev.off()
