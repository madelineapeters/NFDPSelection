library(cowplot)
library(grid)
library(gridExtra)
library(tidyverse)

set.seed(1234)
colors = c("#FC4E07","#E7B800","#00AFBB","#A4C639")
#########################################################
#### Read in data frames for four forms of selection
#########################################################
concave.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_concave_51_sfs.joint.csv") %>% mutate(.,SIM="concave")
constant.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_constant_50_sfs.joint.csv") %>% mutate(.,SIM="constant")
linear.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_linear_52_sfs.joint.csv") %>% mutate(.,SIM="linear") 
convex.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_convex_53_sfs.joint.csv") %>% mutate(.,SIM="convex") 

joint.df = bind_rows(concave.df,constant.df,linear.df,convex.df)

joint.group = joint.df %>% group_by(copies,SIM) %>% summarize(meanFreq=mean(Freq))
joint.group$SIM = factor(joint.group$SIM,levels=c("constant","concave","linear","convex"),labels=c("Constant","Concave NFDP","Linear NFDP","Convex NFDP"))

a = ggplot()+geom_bar(data=joint.group %>% filter(.,copies<=10),aes(x=copies,y=meanFreq,fill=SIM),stat="identity",position="dodge")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10))+
  theme_classic()+labs(x="Frequency",y="Average number of sites",fill=NULL)+theme(legend.position=c(0.8,0.8),axis.title=element_text(size=14),axis.text=element_text(size=12))+scale_fill_manual(values=colors)

b = ggplot()+geom_bar(data=joint.group %>% filter(.,copies>=385),aes(x=copies,y=meanFreq,fill=SIM),stat="identity",position="dodge")+
  scale_x_continuous(breaks=seq(375,400,5))+
  theme_classic()+labs(x="Frequency",y="Average number of sites",fill=NULL)+theme(legend.position="none",axis.title=element_text(size=14),axis.text=element_text(size=12))+scale_fill_manual(values=colors)

pos.grid = ggpubr::ggarrange(a,b,ncol=1,nrow=2,labels=c("A","B"))
pos.grid
ggsave("~/Desktop/PopDataProcessing/SIM_VCFs/dummy/Figures/SFS_grid_full.png",plot=pos.grid,width=5.5,height=8,units=c("in"),dpi=600)

#### Intermediate sites histogram ####
joint.int = joint.df %>% 
  filter(.,copies>10&copies<376) %>% 
  group_by(.,SIM,ID) %>% 
  summarize(.,Count=sum(Freq))
head(joint.int)
joint.int$form = factor(joint.int$SIM,levels=c("constant","concave","linear","convex"),labels=c("Constant","Concave NFDP","Linear NFDP","Convex NFDP"))

hist.int = ggplot()+
  geom_histogram(data=joint.int,aes(x=Count,fill=form),binwidth=1,position="dodge")+
  labs(fill=NULL,y="Frequency",x="Number of intermediate sites")+
  scale_fill_manual(values=colors)+
  theme_classic()+
  theme(strip.background = element_blank(),panel.grid = element_blank(),axis.text.y=element_text(size=12),axis.title = element_text(size=14),strip.text = element_text(size=14), legend.title = element_text(size=14), legend.text = element_text(size=14),panel.spacing = unit(1, "lines"),legend.position=c(0.7,0.8),axis.text.x=element_text(size=12),axis.ticks.x=element_blank())

ggsave("~/Desktop/PopDataProcessing/SIM_VCFs/dummy/Figures/Histogram_IntSeg_full.png",plot=hist.int,width=5.5,height=5,units=c("in"),dpi=600)

#### Segregating sites boxplots ####
joint.seg.box <- joint.df |> group_by(ID,SIM) |> summarize(Count=sum(Freq))
joint.seg.box$SIM <- factor(joint.seg.box$SIM,levels=c("constant","concave","linear","convex"),labels=c("Constant","Concave NFDP","Linear NFDP","Convex NFDP"))

box.seg = ggplot()+geom_boxplot(data=joint.seg.box,aes(x=SIM,y=Count,col=SIM))+
  labs(fill=NULL,y="Segregating sites",x="")+
  scale_color_manual(values=colors)+
  #scale_y_continuous(breaks=seq(5,25,5))+
  theme_classic()+
  theme(legend.position="none",strip.background = element_blank(),panel.grid = element_blank(),axis.text.y=element_text(size=12),axis.title = element_text(size=14),strip.text = element_text(size=14), legend.title = element_text(size=14), legend.text = element_text(size=14),panel.spacing = unit(1, "lines"),axis.text.x=element_text(size=12),axis.ticks.x=element_blank())

ggsave("~/Desktop/PopDataProcessing/SIM_VCFs/dummy/Figures/BoxPlot_SegSites_full.png",plot=box.seg,width=5.5,height=5,units=c("in"),dpi=600)
