##load the table data generated from PPI network analysis 
data1<-read.csv("Growth rate of PPI network parameters.csv")

library(ggplot2)

## Reorder groups 
data1$pathology <-factor(data1$pathology,
                         levels =c("paracancer","adenoma","carcinoma") )

##Data visualization
#Growth.rate.of.avg.number.of.neighbors
ggplot(data1,aes(pathology,Growth.rate.of.avg.number.of.neighbors,group=group,color=cell.type,shape=patient))+
  geom_point(size=4)+
  geom_line(position = position_dodge(0.05),cex=1.5)+
  scale_color_manual(values = c('#00B934',"#619CFF",'#cca4e3',"#c5c56a"))+ 
  scale_shape_manual(values = c(15,16,17,18))+  #shapestyle
  labs(x='Tissue Type',y='Growth Rate of Average Number of Neighbors')+
  theme_test(base_size = 20)+  
  theme(legend.title = element_blank(),
        legend.text = element_text(family = 'serif'), 
        legend.position = "none", #c(.9,.5) 
        legend.direction = "vertical",
        axis.text = element_text(color = 'black',family = 'serif'),
        axis.title = element_text(family = 'serif',size = 14,color = 'black'))

#Growth.rate.of.avg.number.of.neighbors
ggplot(data1,aes(pathology,growth.rate.of.graphindex,group=group,color=cell.type,shape=patient))+
  geom_point(size=4)+
  geom_line(position = position_dodge(0.05),cex=1.5)+
  scale_color_manual(values = c('#00B934',"#619CFF",'#cca4e3',"#c5c56a"))+ 
  scale_shape_manual(values = c(15,16,17,18))+  #shapestyle
  labs(x='Tissue Type',y='Growth Rate of Graphindex')+
  theme_test(base_size = 20)+  
  theme(legend.title = element_blank(),
        legend.text = element_text(family = 'serif'), 
        legend.position = "none",# c(.9,.5)
        legend.direction = "vertical",
        axis.text = element_text(color = 'black',family = 'serif'),
        axis.title = element_text(family = 'serif',size = 14,color = 'black'))
