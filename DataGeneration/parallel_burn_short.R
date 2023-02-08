.libPaths("/project/def-mideon/peter114/Rpkgs")
#.libPaths("/home/peter114/Rpkgs")
suppressMessages(library(parallel))
suppressMessages(library(dplyr))
suppressMessages(library(pryr))

#Read in command line arguments
#SIM, dummy_linear, 5, x, i, dummy, alpha, mfChoice 
args = commandArgs(trailingOnly = TRUE)
#args=c("SIM","burn_constant",5,56,1,"burn",6000,1)
#Data type (OBS or SIM)
type = as.character(args[1])

#Specify selection form to simulate
SIM = as.character(args[2])

#Specify length of genome segment to simulate
kb = as.integer(args[3])

#Specify parameter combination index
x = as.integer(args[4])

#Specify population (dummy)
POPlist = as.character(args[6])

#Specify alpha (alpha = 2N*s)
alpha = as.numeric(args[7])

#Specify initial allele frequency
mf.choice = as.integer(args[8])
mf = c(1/60000,10/60000)[mf.choice]

#Specify ID of simulation
sim = as.integer(args[5])

FREQlist = c(98)
fGENlist = c(1)

#clusterFUN = function(i,x,mf,POPlist,minPOS,maxPOS,kb){
simFUN = function(i,x,POPlist,kb,SIM,FREQlist,fGENlist,type,alpha,mf){
  
  .libPaths("/project/def-mideon/peter114/Rpkgs")
  #.libPaths("/home/peter114/Rpkgs")
  library(vcfR)
  library(PopGenome)
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(caTools)
  library(rehh)
  library(ape) 
  #Write function to make iteration specific directory
  mkDIR = function(POP,b,i,spec,type) {
    POPpath = paste0(type,"_",POP)
    output = system2("mkdir", c(paste0("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/",POP,"_",b,"_",i,"_",spec)),stdout=T)
  }
  #Write function to sort and slice .vcf
  sliceVCF = function(POP,b,i,kb,spec,type){
    POPpath = paste0(type,"_",POP)
    output = system2("bash",c(paste0("/project/def-mideon/peter114/tabixVCF_check.sh"),POP,b,i,kb,spec,type),stdout=T)
  }
  #Write function to remove iteration specific directory
  rmDIR = function(POP,b,i,spec,type){
    POPpath = paste0(type,"_",POP)
    output = system2("rm", c("-r",paste0("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/",POP,"_",b,"_",i,"_",spec),"2>","/dev/null"),stdout=T)
  }
  
  #Write function to run SLiM from command line
  runSLiM = function(x) {
    seed = x[1]
    alpha = x[2]
    y = x[3]
    i = x[4]
    model = x[5]
    type = x[6]
    mf = x[7]
    output = system2("/home/peter114/SLiM_old/install/bin/slim", c("-d", paste0("alpha=", alpha),"-d",paste0("i=", i),"-d", paste0("x=", y),"-d", paste0("mf=", mf),"-s", seed, paste0("/project/def-mideon/peter114/",type,"_",model,"_burn.slim")), stdout=T)
  } 
  
  #Write function to list summary statistics
  sum.stat.pos = function(x,i,POP,kb,spec,type){
    
    POPpath=paste0(type,"_",POP)
    
    #Define number kb for sequence (used for IHH)
    l = kb
    
    #################################
    ##### Load (and clean) data #####
    #################################
    
    #### Data for PopGenome method ####
    #GENOME.class = readData("~/Desktop/PopDataProcessing/test/tmp",format="VCF")
    GENOME.class.sel = readData(paste("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/",POP,"_",x,"_",i,"_",spec,sep=""),format="VCF")
    
    #### Data for original method ####
    #vcf.string = "~/Desktop/PopDataProcessing/test/tmp/wAGR2_20_203875_12000_sliced.vcf"
    vcf.string = paste("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/",POP,"_",x,"_",i,"_",spec,"/",POP,"_",x,"_",i,"_",spec,"_sorted.vcf",sep="")
    
    vcf = read.vcfR(vcf.string, verbose = FALSE)
    
    #Get positions for sites variant in VCF sample
    POS.df = vcf@fix %>% as.data.frame %>% dplyr::select(.,POS)
    
    #Get genotype data frame (row = position, column = individual)
    GT.df = vcf@gt %>% as.data.frame %>% dplyr::select(.,-FORMAT)
    
    #Pivot genotype data frame longer and separate genotype info into two columns (columns = POS,ID,allele copy 1, allele copy 2)
    GT.long = bind_cols(POS.df,GT.df) %>% pivot_longer(.,cols=names(GT.df),names_to="ID",values_to="geno") %>% separate(.,col="geno",into=c("hap1","hap2"),convert=TRUE)
    
    #Pivot genotype data wider again
    GT.wide = GT.long %>% pivot_longer(.,hap1:hap2,names_to="hap",values_to="allele") %>% pivot_wider(.,names_from=POS,values_from=allele)
    
    #Calculate average pairwise differences for each site and bind data frame to position data frame
    pi.vals = sapply(3:ncol(GT.wide),FUN=function(c){
      as.numeric(dist(GT.wide[,c])) %>% mean
    }) %>% as.data.frame %>% bind_cols(POS.df,.)
    names(pi.vals) = c("POS","pi")
    
    #Calculate number of variant copies at each site and bid data frame to position data frame
    seg.vals = GT.wide[,3:ncol(GT.wide)] %>% colSums %>% as.data.frame %>% bind_cols(POS.df,.) %>% remove_rownames()
    names(seg.vals) = c("POS","seg")
    
    #Join pairwise diversity and SFS data frames
    POS.stat = full_join(pi.vals,seg.vals,by="POS")
    
    #Convert position from character to integer
    POS.stat$POS = POS.stat$POS %>% type.convert
    
    #Create list of all positions simulated
    #POS.full.list = seq(1,5000,1) %>% as.data.frame()
    #names(POS.full.list) = "POS"
    #POS.full = left_join(POS.full.list,POS.stat,by="POS")
    #POS.full[is.na(POS.full)] = 0
    POS.full = POS.stat
    ##########################################
    ##### Frequency of Duffy null allele #####
    ##########################################
    position1 = 2501
    Freq.sel = filter(POS.full,POS==position1)$seg/(2*ncol(GT.df))
    
    #################################
    ##### Neutrality statistics #####
    #################################
    GENOME.class.sel = GENOME.class.sel %>% neutrality.stats() %>% detail.stats(biallelic.structure=TRUE) 
    seqnum.sel = get.biallelic.matrix(GENOME.class.sel,1) %>% nrow
    SFS.PG.sel = table(sapply(GENOME.class.sel@region.stats@minor.allele.freqs[[1]],FUN=function(x){x*seqnum.sel}))
    
    SegSites.PG.sel = get.neutrality(GENOME.class.sel)[[1]][1,2]
    TajD.PG.sel = get.neutrality(GENOME.class.sel)[[1]][1,1]
    Singletons.PG.sel = SFS.PG.sel[1]
    Doubletons.PG.sel = SFS.PG.sel[2]
    Fixed.PG.sel = l*1000 - GENOME.class.sel@n.biallelic.sites[1]
   
    #################################
    ##### Diversity statistics #####
    #################################
    GENOME.class.sel = GENOME.class.sel %>% diversity.stats(.,pi=TRUE)
    pi.PG.sel = get.diversity(GENOME.class.sel)[[1]][1,1]

    ########################
    ##### H statistics #####
    ########################
    
    #### PopGenome method ####
    source("/project/def-mideon/peter114/get.Hstats.R")
    Hstat.PG.sel = get.Hstats(GENOME.class.sel)[1,]
    
    ################################
    #### EHH and IHH statistics ####
    ################################
    
    HapNum.PG.sel = GENOME.class.sel@region.stats@haplotype.counts[[1]] %>% length
    
    rehh.data = data2haplohh(hap_file = vcf.string,
                             polarize_vcf = FALSE, remove_multiple_markers = TRUE)
    
    mrk = which((abs(rehh.data@positions - 2501)) == min(abs(rehh.data@positions - 2501)))
    
    selIndices = which(rehh.data@positions<=5001)
    
    ehh.obj.sel = calc_ehh(rehh.data,mrk=mrk,limehh=0)
    
    ehh.D.O.sel = ehh.obj.sel$ehh$EHH_D[c(selIndices[1],selIndices[length(selIndices)])] %>% mean
    ehh.A.O.sel = ehh.obj.sel$ehh$EHH_A[c(selIndices[1],selIndices[length(selIndices)])] %>% mean
    ihh.D.O.sel = trapz(c(2501-(l*1000/2),ehh.obj.sel$ehh$POSITION[selIndices],2501+(l*1000/2)),c(first(ehh.obj.sel$ehh$EHH_D[selIndices]),ehh.obj.sel$ehh$EHH_D[selIndices],last(ehh.obj.sel$ehh$EHH_D[selIndices])))
    ihh.A.O.sel = trapz(c(2501-(l*1000/2),ehh.obj.sel$ehh$POSITION[selIndices],2501+(l*1000/2)),c(first(ehh.obj.sel$ehh$EHH_A[selIndices]),ehh.obj.sel$ehh$EHH_A[selIndices],last(ehh.obj.sel$ehh$EHH_A[selIndices])))
    
    ehh.O.sel = max(ehh.D.O.sel,ehh.A.O.sel)
    ihh.O.sel = max(ihh.D.O.sel,ihh.A.O.sel)
    
    ################
    #### Output ####
    ################
    
    out1 = c(pi.PG.sel,SegSites.PG.sel,TajD.PG.sel,HapNum.PG.sel,Fixed.PG.sel,Singletons.PG.sel,Doubletons.PG.sel,Hstat.PG.sel,Freq.sel,Singletons.PG.sel/Fixed.PG.sel,ehh.O.sel,ihh.O.sel) %>% unname
    #out1 = c(pi.PG,SegSites.PG,TajD.PG,HapNum.PG,Fixed.PG,Singletons.PG,Doubletons.PG,Hstat.PG[1],Hstat.PG[4],Singletons.PG/Fixed.PG) %>% unname
    
    out2 = POS.full
    
    return(c(out1,out2))
    
  }
  
  finalFreq = 0
  while (finalFreq < 0.95){
    seed = sample(1:10000000, 1)
  
  #Run SLiM and calculate summary statistics inside of tryCatch set-up
  tmp = try({
    
    #Calculate summary statistics
    sum.stat.joint = c()
    pos.stat.joint = data.frame()
    for (POP in POPlist){
      for (freq in FREQlist){
        sliceVCF(POP,x,i,kb,freq,type)
        sum.stat.list = sum.stat.pos(x,i,POP,kb,freq,type)
        #rmDIR(POP,x,i,gen,type)
        sum.stat.tmp = sum.stat.list[1:15] %>% unlist
        pos.stat.tmp = as.data.frame(sum.stat.list[16:18])
        pos.stat.tmp$spec = freq
        pos.stat.tmp$alpha = alpha
        pos.stat.tmp$mf = mf
        pos.stat.tmp$ID = i
        sum.stat.joint = append(sum.stat.joint,sum.stat.tmp)
        pos.stat.joint = bind_rows(pos.stat.joint,pos.stat.tmp)
        print(paste("Done",freq,sep=" "))
      }
    }
    sum.stat.joint = c(mf,alpha,sum.stat.joint,seed)
    finalFreq = sum.stat.tmp[12]
    print(paste0("Final frequency is ",finalFreq))
  })
  if(inherits(try, "try-error")) {
    recoverFromFailure()
  }  
  
  } #end of while statement
	  
  out = c(sum.stat.joint,pos.stat.joint)
  print("Exiting")
  return(out)
  
}

#Run SLiM simulation and return summary statistics
## i,x,POPlist,kb,SIM,FREQlist,fGENlist,type,alpha,mf
detest.pair = simFUN(sim,x,POPlist,kb,SIM,FREQlist,fGENlist,type,alpha,mf)
print(length(detest.pair))
detest = detest.pair[1:18]

sum.stat.df = as.data.frame(matrix(unlist(detest), nrow = 1, ncol = round(length(unlist(detest))), byrow = T))
sub.namelist = c("pi","SegSites","TajD","HapNum","Fixed","Singletons","Doubletons","H1","H2","H2.H1","H12","Freq","Singletons.Fixed","ehh","ihh")
nameslist = c("mf","s",apply(expand.grid(apply(expand.grid(sub.namelist,c("sel")), 1, paste, collapse="."),c(FREQlist)),1,paste, collapse="."),"seed")
names(sum.stat.df) = nameslist

POPpath="SIM_burn"
POP="burn"

df = sum.stat.df
df$ID = sim

if (sim > 1){
  
  print("Write as temporary .csv")
  write.csv(df,paste("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/sumstats_",x,"_",POP,"tmp.csv",sep=""),row.names=FALSE)
  
  print("Remove header")
  system2("awk",c("'FNR > 1'",paste("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/sumstats_",x,"_",POP,"tmp.csv",sep=""),">",paste("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/sumstats_",x,"_",POP,"noheading.csv",sep="")),stdout=T)
  
  print("Glue to master")
  system2("cat",c(paste("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/sumstats_",x,"_",POP,"noheading.csv",sep=""),">>",paste("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/",SIM,"_",kb,"_stats_",x,"_",POP,"Freq.csv",sep="")),stdout=T)
  
  print("Remove temporary file")
  system2("rm",c(paste("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/sumstats_",x,"_",POP,"tmp.csv",sep="")),stdout=T)
  
  print("Remove file without header")
  system2("rm",c(paste("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/sumstats_",x,"_",POP,"noheading.csv",sep="")),stdout=T)
  
} else {
  
  print("Write .csv")
  write.csv(df,paste("/scratch/peter114/PopDataProcessing/SIM_VCFs/",POPpath,"/",SIM,"_",kb,"_stats_",x,"_",POP,"Freq.csv",sep=""),row.names=FALSE)
  
}

