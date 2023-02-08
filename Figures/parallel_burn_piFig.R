library(tidyverse)

set.seed(1234)
colors = c("#FC4E07","#E7B800","#00AFBB","#A4C639")
#########################################################
#### Read in data frames for four forms of selection
#########################################################
concave.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_concave_51_sum.stat.joint.csv") %>% mutate(.,SIM="concave") %>% sample_n(1000)
constant.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_constant_50_sum.stat.joint.csv") %>% mutate(.,SIM="constant") %>% sample_n(1000)
linear.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_linear_52_sum.stat.joint.csv") %>% mutate(.,SIM="linear") %>% sample_n(1000)
convex.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_convex_53_sum.stat.joint.csv") %>% mutate(.,SIM="convex") %>% sample_n(1000)

joint.df <- bind_rows(concave.df,constant.df,linear.df,convex.df) |> dplyr::select(pi,SIM)
joint.df$SIM <- factor(joint.df$SIM,levels=c("constant","concave","linear","convex"),labels=c(expression(atop("Constant,",s[0]==0.025)),expression(atop("Concave NFDP,",s[0]==0.076)),expression(atop("Linear NFDP,",s[0]==0.132)),expression(atop("Convex NFDP,",s[0]==0.235))))

box.pi = ggplot()+geom_boxplot(data=joint.df,aes(x=SIM,y=pi,col=SIM))+
  labs(fill=NULL,y=expression("Average pairwise diversity"~pi),x="")+
  scale_color_manual(values=colors)+
  ggtitle(expression(s[pi]==0.025))+
  #scale_y_continuous(breaks=seq(5,25,5))+
  scale_x_discrete(labels = ggplot2:::parse_safe)+
  theme_classic()+
  geom_hline(yintercept=0.6,size=2)+
  theme(legend.position="none",strip.background = element_blank(),panel.grid = element_blank(),axis.text.y=element_text(size=12),axis.title = element_text(size=14),strip.text = element_text(size=14), legend.title = element_text(size=14), legend.text = element_text(size=14),panel.spacing = unit(1, "lines"),axis.text.x=element_text(size=12),axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5))

ggsave("~/Desktop/PopDataProcessing/SIM_VCFs/dummy/Figures/BoxPlot_pi_full.png",plot=box.pi,width=5.5,height=5,units=c("in"),dpi=600)


