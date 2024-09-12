
library(dplyr)
library(ggplot2)
library(ggpubr)
df2 <- data.frame(condition = rep(c("con-20","con-20","20um","20um","con-50","con-50","50um","50um"),3), 
                  dose = rep(c("3G","3GL","MYC"),each=8),
                  expression = c(1,1,0.871994571,0.894672618,1,1,0.692836927,0.757292096,1,1,0.838346654,1.079488407,1,1,0.662313535,0.678724027,1,1,0.889017747,1.103610427,
                                  1,1,0.702229849,0.809838475),
                  group = c(rep(20,4),rep(50,4),rep(20,4),rep(50,4),rep(20,4),rep(50,4)))
df2$ add <- as.factor(paste(df2$condition,df2$dose, sep="_"))
df3 <- df2 %>% group_by(add) %>% summarise_at (vars(expression), list(mean= mean)) 
df4 <- df2 %>% group_by(add) %>% summarise_at (vars(expression), list(SD = sd)) 
df3 <- cbind(df3,SD=df4$SD)
df3$group <- as.factor(rep(c("3G_20","3GL_20","Myc_20","3G_50","3GL_50","Myc_50"),2))
df3$treat <- rep(c("treatment","control"),each=6)
df3$rep1 <- c(0.871994571, 0.838346654,0.889017747,0.692836927,0.662313535,0.702229849,rep(1,6))
df3$rep2 <- c(0.894672618,1.079488407,1.103610427,0.757292096,0.678724027,0.809838475,rep(1,6))

quartz(width=5.3, height = 4.5)
p<- ggplot2::ggplot(data=df3, aes(x=group,y=mean, fill=treat)) +
                geom_bar(stat="identity",color="black", position=position_dodge())+
                 geom_errorbar(aes(ymin = mean-SD,ymax = mean+SD), width= .2,
                    position=position_dodge(.9)) +
                ylim(0,1.5) +
  geom_dotplot(data = df3, aes(x=group,y=rep1,
               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=.9),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5) +
  geom_dotplot(data = df3, aes(x=group,y=rep2,
                               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=.9),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5)

p<- p+labs(x= "THP-1 + KJ-Pyr-9",y = "relative RNA level")
print(p)
ggsave("/Users/vanbortlelab/Desktop/Sihang_MacBook_Pro/kvlab/EXPT/plots/myci_barplot.png", plot = last_plot())



#WB results####
dff <- data.frame(group=c(rep("MYC",2),rep("MAX",2)), treat = rep(c("control","treatment"),2),mean=c(1,0.721022317,1,0.615511425),SD=c(0,0.038339577,0,0.046339204),rep1 = c(1,0.693912142,1,0.64827819),rep2 = c(1,0.748132492,1,0.582744659))
quartz(width=3.5, height = 4.5)
q<- ggplot2::ggplot(data=dff, aes(x=group,y=mean, fill=treat)) +
  geom_bar(stat="identity",color="black", position=position_dodge())+
  geom_errorbar(aes(ymin = mean-SD,ymax = mean+SD), width= .2,
                position=position_dodge(.9)) +
  ylim(0,1.5) +
  geom_dotplot(data = dff, aes(x=group,y=rep1,
                               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=.9),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5) +
  geom_dotplot(data = dff, aes(x=group,y=rep2,
                               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=.9),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5) 
  

q<- q+labs(x= "THP-1 + KJ-Pyr-9",y = "relative protein level")
print(q) 

figure <- ggarrange(p,q,ncol=2,nrow=1,common.legend = T, widths = c(6,3),legend = "right")
quartz(height=4.5,width= 7.3)
print(figure)
ggsave("/Users/vanbortlelab/Desktop/Sihang_MacBook_Pro/kvlab/EXPT/plots/myci_barplot.final.png", plot = last_plot())

#KJ,F4####
df5 <- data.frame(condition = rep(c("KJ-Pyr9","10058-F4","DMSO"),each=9), 
                  dose = rep(c("MYC","CDK4","PTMA"),9),
                  expression = c(1.237770378,1.519777103,1.266265698,1.324068618,2.049423015,1.309525772,1.296016927,1.899956343,1.316639387,0.522461813,1.059643669,0.542781615,0.656658971,1.316828803,
                                 0.49377767,0.500057544,1.02407597,0.769349943,rep(1,9)))
df5$ add <- as.factor(paste(df5$condition,df5$dose, sep="_"))
df6 <- df5 %>% group_by(add) %>% summarise_at (vars(expression), list(mean= mean)) 
df7 <- df5 %>% group_by(add) %>% summarise_at (vars(expression), list(SD = sd)) 
df6 <- cbind(df6,SD=df7$SD)
df6$group <- factor(rep(c("CDK4","MYC","PTMA"),3),levels = c("MYC","CDK4","PTMA"))
df6$treat <- factor(rep(c("10058-F4","DMSO","KJ-Pyr9"),each=3),levels =c("DMSO","KJ-Pyr9","10058-F4"))
df6$rep1 <- c(1.059643669, 0.522461813, 0.542781615,rep(1,3),1.519777103,1.237770378,1.266265698)
df6$rep2 <- c(1.316828803,0.656658971,0.49377767,rep(1,3),2.049423015,1.324068618,1.309525772)
df6$rep3 <- c(1.02407597,0.542781615,0.769349943,rep(1,3),1.899956343,1.296016927,1.316639387)
 quartz(width=5.3, height = 4.5)
p<- ggplot2::ggplot(data=df6, aes(x=group,y=mean, fill=treat)) +
  geom_bar(stat="identity",color="black", position=position_dodge())+
  geom_errorbar(aes(ymin = mean-SD,ymax = mean+SD), width= .2,
                position=position_dodge(width=1)) + scale_fill_brewer(palette="Set2")+
  ylim(0,2.5) +
  geom_dotplot(data = df6, aes(x=group,y=rep1,
                               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=1),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5) +
  geom_dotplot(data = df6, aes(x=group,y=rep2,
                               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=1),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5) +
  geom_dotplot(data = df6, aes(x=group,y=rep3,
                               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=1),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5) 
p<- p+labs(x= "THP1",y = "relative RNA level")
print(p)
ggsave("/Users/vanbortlelab/Desktop/Sihang_MacBook_Pro/kvlab/EXPT/plots/myci_kj_F4_barplot.png", plot = last_plot())

#ASO results####
dff <- data.frame(group=c(rep("MYC",2),rep("MAX",2)), treat = rep(c("control","treatment"),2),mean=c(1,0.721022317,1,0.615511425),SD=c(0,0.038339577,0,0.046339204),rep1 = c(1,0.693912142,1,0.64827819),rep2 = c(1,0.748132492,1,0.582744659))
quartz(width=3.5, height = 4.5)
q<- ggplot2::ggplot(data=dff, aes(x=group,y=mean, fill=treat)) +
  geom_bar(stat="identity",color="black", position=position_dodge())+
  geom_errorbar(aes(ymin = mean-SD,ymax = mean+SD), width= .2,
                position=position_dodge(.9)) +
  ylim(0,1.5) +
  geom_dotplot(data = dff, aes(x=group,y=rep1,
                               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=.9),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5) +
  geom_dotplot(data = dff, aes(x=group,y=rep2,
                               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=.9),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5) 


q<- q+labs(x= "THP-1 + KJ-Pyr-9",y = "relative protein level")
print(q) 

figure <- ggarrange(p,q,ncol=2,nrow=1,common.legend = T, widths = c(6,3),legend = "right")
quartz(height=4.5,width= 7.3)
print(figure)
ggsave("/Users/vanbortlelab/Desktop/Sihang_MacBook_Pro/kvlab/EXPT/plots/myci_barplot.final.png", plot = last_plot())

#KJ,F4####
df5 <- data.frame(condition = rep(c("KJ-Pyr9","10058-F4","DMSO"),each=9), 
                  dose = rep(c("MYC","CDK4","PTMA"),9),
                  expression = c(1.237770378,1.519777103,1.266265698,1.324068618,2.049423015,1.309525772,1.296016927,1.899956343,1.316639387,0.522461813,1.059643669,0.542781615,0.656658971,1.316828803,
                                 0.49377767,0.500057544,1.02407597,0.769349943,rep(1,9)))
df5$ add <- as.factor(paste(df5$condition,df5$dose, sep="_"))
df6 <- df5 %>% group_by(add) %>% summarise_at (vars(expression), list(mean= mean)) 
df7 <- df5 %>% group_by(add) %>% summarise_at (vars(expression), list(SD = sd)) 
df6 <- cbind(df6,SD=df7$SD)
df6$group <- factor(rep(c("CDK4","MYC","PTMA"),3),levels = c("MYC","CDK4","PTMA"))
df6$treat <- factor(rep(c("10058-F4","DMSO","KJ-Pyr9"),each=3),levels =c("DMSO","KJ-Pyr9","10058-F4"))
df6$rep1 <- c(1.059643669, 0.522461813, 0.542781615,rep(1,3),1.519777103,1.237770378,1.266265698)
df6$rep2 <- c(1.316828803,0.656658971,0.49377767,rep(1,3),2.049423015,1.324068618,1.309525772)
df6$rep3 <- c(1.02407597,0.542781615,0.769349943,rep(1,3),1.899956343,1.296016927,1.316639387)
quartz(width=5.3, height = 4.5)
p<- ggplot2::ggplot(data=df6, aes(x=group,y=mean, fill=treat)) +
  geom_bar(stat="identity",color="black", position=position_dodge())+
  geom_errorbar(aes(ymin = mean-SD,ymax = mean+SD), width= .2,
                position=position_dodge(width=1)) + scale_fill_brewer(palette="Set2")+
  ylim(0,2.5) +
  geom_dotplot(data = df6, aes(x=group,y=rep1,
                               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=1),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5) +
  geom_dotplot(data = df6, aes(x=group,y=rep2,
                               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=1),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5) +
  geom_dotplot(data = df6, aes(x=group,y=rep3,
                               group = interaction(group, treat)),
               position_jitterdodge(jitter.width=.1,dodge.width=1),binaxis = "y", alpha=.8,stackdir="center",stackratio = 1,dotsize = 0.5) 
p<- p+labs(x= "THP1",y = "relative RNA level")
print(p)
ggsave("/Users/vanbortlelab/Desktop/Sihang_MacBook_Pro/kvlab/EXPT/plots/myci_kj_F4_barplot.png", plot = last_plot())
