library(tidyverse)
library(mixOmics)
library(cowplot)
library(grid)
library(gridExtra)

set.seed(1234)
colors = c("#FC4E07","#E7B800","#00AFBB","#A4C639")
#########################################################
#### Read in data frames for four forms of selection
#########################################################
concave.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_concave_51_sum.stat.joint.csv") %>% mutate(.,SIM="concave") %>% sample_n(1000)
constant.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_constant_50_sum.stat.joint.csv") %>% mutate(.,SIM="constant") %>% sample_n(1000)
linear.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_linear_52_sum.stat.joint.csv") %>% mutate(.,SIM="linear") %>% sample_n(1000)
convex.df = read.csv("~/Desktop/PopDataProcessing/SIM_VCFs/SIM_burn/burn_convex_53_sum.stat.joint.csv") %>% mutate(.,SIM="convex") %>% sample_n(1000)

joint.df = bind_rows(concave.df,constant.df,linear.df,convex.df)
joint.df$SIM = factor(joint.df$SIM,levels=c("constant","concave","linear","convex"),labels=c("Constant","Concave\nNFDP","Linear\nNFDP","Convex\nNFDP"))

#############################################
#### Run initial PLS-DA to consider loadings
#############################################
X = joint.df %>% dplyr::select(.,SegSites:ihh)
Y = joint.df %>% dplyr::select(.,SIM) %>% mutate_if(., is.character, as.factor) %>% c

#Choose number of components
cluster.perf = plsda(X, Y$SIM,ncomp=ncol(X))
ncomp = which(cluster.perf$prop_expl_var$X %>% cumsum > 0.95) %>% min
ncomp2 = 3 #number of classes - 1
#PLSDA with chosen number of components
cluster.res = plsda(X, Y$SIM,ncomp=ncomp2)

indivPlot = plotIndiv(cluster.res,ellipse=TRUE,legend=TRUE,title="",col.per.group=colors[c(2,1,3,4)],ind.names=FALSE,legend.title=NULL)
loadingsPlot = plotLoadings(cluster.res,comp=1,legend.color=colors[c(1,2,3,4)],legend=TRUE,contrib = 'max', method = 'mean',border=T,legend.title="Form",title="Component 1")

ggsave("~/Desktop/PopDataProcessing/SIM_VCFs/dummy/Figures/Loadings_PLSDA_full.png",plot=loadingsPlot,width=5.5,height=5,units=c("in"),dpi=600)

###############################################################
#### Classification bootstrapping for constant versus NFDP ####
###############################################################
joint.pred = data.frame()
for (i in 1:20){
  
  set.seed(i)
  
  #Set up training and testing data
  test.df = data.frame()
  train.df = data.frame()
  
  for (sim in unique(joint.df$SIM)){
    
    if (sim == "Constant"){
      train.list = sample(1:1000,900,replace=FALSE) 
      test.tmp = joint.df %>% filter(.,SIM==sim)
      train.tmp = test.tmp[train.list,]
      test.tmp = test.tmp[-train.list,] 
    } else if (sim == "Concave\nNFDP") {
      train.list = sample(1:1000,300,replace=FALSE)
      test.tmp = joint.df %>% filter(.,SIM==sim)
      train.tmp = test.tmp[train.list,]
      test.tmp = test.tmp[-train.list,] %>% sample_n(34)
      train.tmp$SIM = "NFDP"
      test.tmp$SIM = "NFDP"
    } else {
      train.list = sample(1:1000,300,replace=FALSE)
      test.tmp = joint.df %>% filter(.,SIM==sim)
      train.tmp = test.tmp[train.list,]
      test.tmp = test.tmp[-train.list,] %>% sample_n(33)
      train.tmp$SIM = "NFDP"
      test.tmp$SIM = "NFDP"
    }
    train.df = train.df %>% bind_rows(.,train.tmp)
    test.df = test.df %>% bind_rows(.,test.tmp)
  }  
  
  #Set up PLS-DA with training data
  X = train.df %>% dplyr::select(.,SegSites:ihh)
  Y = train.df %>% dplyr::select(.,SIM) %>% mutate_if(., is.character, as.factor) %>% c
  
  cluster.perf = plsda(X, Y$SIM,ncomp=ncol(X))
  ncomp = which(cluster.perf$prop_expl_var$X %>% cumsum > 0.95) %>% min
  
  cluster.res = plsda(X, Y$SIM,ncomp=ncomp)
  
  #Predict for testing data
  pred = predict(cluster.res,test.df %>% dplyr::select(.,SegSites:ihh))
  
  # store prediction for the 3rd component
  prediction = pred$class$mahalanobis.dist[,ncomp]
  cluster.pred = test.df %>% dplyr::select(.,SIM) %>% droplevels()
  cluster.pred$prediction = prediction
  
  pred.tab = table(cluster.pred$SIM,cluster.pred$prediction) %>% as.data.frame
  names(pred.tab) = c("SIM","SIM.pred","Freq")
  pred.tab$iteration = i
  
  joint.pred = bind_rows(joint.pred,pred.tab)
  
}
pred.tab$SIM.pred = factor(pred.tab$SIM.pred,levels=c("Constant","NFDP"))

a = ggplot(pred.tab, aes(fill=SIM.pred, y=Freq/100, x=SIM)) + 
  geom_bar(position="fill", stat="identity")+theme_minimal()+
  #scale_fill_manual(values=form.alpha.colors,labels=sapply(levels(pred.tab$form.alpha.pred),FUN=function(x){parse(text=x)}))+
  scale_fill_manual(values=c(colors[1],"grey"))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0),labels=c("0%","25%","50%","75%","100%"))+
  #scale_x_discrete(breaks=levels(pred.tab$form),FUN=function(x){parse(text=x)}))+
  ggtitle("")+
  labs(x=NULL,y=NULL,fill=NULL,parse=TRUE)+
  theme(strip.background = element_blank(),panel.grid = element_blank(),axis.text=element_text(size=12),axis.text.x=element_text(size=12,angle = 0),axis.title = element_text(size=12),strip.text = element_text(size=12), legend.title = element_text(size=10), legend.text = element_text(size=14),panel.spacing = unit(1, "lines"),legend.position="right")


#####################################################
#### Classification bootstrapping for all four forms
#####################################################
joint.pred = data.frame()
for (i in 1:20){
  
  set.seed(i)
  
  #Set up training and testing data
  test.df = data.frame()
  train.df = data.frame()
  
  for (sim in unique(joint.df$SIM)){
    
    train.list = sample(1:1000,900,replace=FALSE) 
    test.tmp = joint.df %>% filter(.,SIM==sim)
    train.tmp = test.tmp[train.list,]
    test.tmp = test.tmp[-train.list,]
    
    train.df = train.df %>% bind_rows(.,train.tmp)
    test.df = test.df %>% bind_rows(.,test.tmp)
  }  
  
  #Set up PLS-DA with training data
  X = train.df %>% dplyr::select(.,SegSites:ihh)
  Y = train.df %>% dplyr::select(.,SIM) %>% mutate_if(., is.character, as.factor) %>% c
  
  cluster.perf = plsda(X, Y$SIM,ncomp=ncol(X))
  ncomp = which(cluster.perf$prop_expl_var$X %>% cumsum > 0.95) %>% min
  
  cluster.res = plsda(X, Y$SIM,ncomp=ncomp)
  
  #Predict for testing data
  pred = predict(cluster.res,test.df %>% dplyr::select(.,SegSites:ihh))
  
  # store prediction for the 3rd component
  prediction = pred$class$mahalanobis.dist[,ncomp]
  cluster.pred = test.df %>% dplyr::select(.,SIM) %>% droplevels()
  cluster.pred$prediction = prediction
  
  pred.tab = table(cluster.pred$SIM,cluster.pred$prediction) %>% as.data.frame
  names(pred.tab) = c("SIM","SIM.pred","Freq")
  pred.tab$iteration = i
  
  joint.pred = bind_rows(joint.pred,pred.tab)
  
}
pred.tab$SIM.pred = factor(pred.tab$SIM.pred,levels=c("Constant","Concave\nNFDP","Linear\nNFDP","Convex\nNFDP"))

b = ggplot(pred.tab, aes(fill=SIM.pred, y=Freq/100, x=SIM)) + 
  geom_bar(position="fill", stat="identity")+theme_minimal()+
  #scale_fill_manual(values=form.alpha.colors,labels=sapply(levels(pred.tab$form.alpha.pred),FUN=function(x){parse(text=x)}))+
  scale_fill_manual(values=colors)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0),labels=c("0%","25%","50%","75%","100%"))+
  #scale_x_discrete(breaks=levels(pred.tab$form),FUN=function(x){parse(text=x)}))+
  ggtitle("")+
  labs(x=NULL,y=NULL,fill=NULL,parse=TRUE)+
  theme(strip.background = element_blank(),panel.grid = element_blank(),axis.text=element_text(size=12),axis.text.x=element_text(size=12,angle = 0),axis.title = element_text(size=12),strip.text = element_text(size=12), legend.title = element_text(size=10), legend.text = element_text(size=14),panel.spacing = unit(1, "lines"),legend.position="right")

joint.pred %>% group_by(.,SIM,SIM.pred) %>% summarize(.,Mean=mean(Freq/50),SD=sd(Freq/50)) %>% mutate(.,SE=SD/sqrt(10))

pos.grid = ggpubr::ggarrange(a,b,ncol=1,nrow=2,labels=c("A","B"))
pos.grid
ggsave("~/Desktop/PopDataProcessing/SIM_VCFs/dummy/Figures/Stacked_PLSDA_full.png",plot=pos.grid,width=5.5,height=8,units=c("in"),dpi=600)
