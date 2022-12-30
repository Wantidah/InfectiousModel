#patchwork
library(patchwork)
#non-scale
png("gaur_1sim_patch.png",width = 15, height = 30, units = 'cm', res = 600)
p<-p1/p2/p3/p4/p5/p6/p7 + plot_layout(guides = 'collect') +
  plot_annotation(title = "1 simulation")
dev.off()

#scale
png("gaur_1sim_patch_scale.png",width = 15, height = 30, units = 'cm', res = 600)
ps<-p1s/p2s/p3s/p4s/p5s/p6s/p7s + plot_layout(guides = 'collect') +
  #plot_annotation(title = "100 simulation")
  labs(subtitle = "100 simulation")
dev.off()

png("test_patchwork2.png",width = 30, height = 40, units = 'cm', res = 600)
pp <- p1/p2/p3/p4/p5/p6/p7 | pl1/pl2/pl3/pl4/pl5/pl6/pl7 
pp 
dev.off()

pp1 <- p1/p2/p3/p4/p5/p6/p7 +  plot_annotation(title = '1 simulation') 
pp2<- pl1/pl2/pl3/pl4/pl5/pl6/pl7 + plot_annotation(title = "100 simulations") 
pp3<-pp1|pp2 + plot_layout(guides='collect')

png("test_patchwork3.png",width = 30, height = 40, units = 'cm', res = 600)
print(pp3)
dev.off()

