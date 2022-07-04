# Read Data 
library(lmerTest)
library(tidyverse)
library(ggseg3d)
rm(list=ls())
behavior_data <- read.csv("D:\\Rdata\\group_SLIM_stats\\PlasticityAnalysis\\BehaviorData.csv")
behavior_data <- behavior_data[order(behavior_data$ID),]
DataDir <- "D:\\Rdata\\group_SLIM_stats\\PlasticityAnalysis\\"
SaveDir <- "D:\\Rdata\\group_SLIM_stats\\PlasticityAnalysis\\ResultStats\\"
#######################################################
#                   Volume Analysis                   #
#######################################################
if(1){
  setwd("D:\\Rdata\\group_SLIM_stats\\vol")
  lh_vol_data <- read.csv("lh_brain_data.csv")
  rh_vol_data <- read.csv("rh_brain_data.csv")
  lh_vol_data <- lh_vol_data[order(lh_vol_data[,2]),]
  rh_vol_data <- rh_vol_data[order(rh_vol_data[,2]),]

  
  # data clean
  ModelData <- cbind(behavior_data,lh_vol_data,rh_vol_data)
  ModelData <- ModelData[,-c(1,51:54,235:240,421,422)]
  ##平局间隔
  temp <- subset(ModelData,ModelData$label=="Time3")
  
  temp <- aggregate(ModelData,by=list(ModelData$ID),max)
  temp <- subset(temp,temp$time!=0)
  mean(temp$time)
  #AgeMean(i)
  AgeMeanSub <- aggregate(ModelData$Age1,by=list(ModelData$ID),mean)
  colnames(AgeMeanSub)[1:2] <- c("ID","AgeMeanSub")
  ModelData <- inner_join(ModelData,AgeMeanSub,by="ID")
  #time(ij)=Age(ij)-AgeMean(i)
  #ModelData$timeLong <- ModelData$Age1-ModelData$AgeMeanSub
  #AgeMean
  ageMean <- mean(ModelData$Age1)
  #age_mean(i)=AgeMean(i)-AgeMean
  ModelData$age_mean <- ModelData$AgeMeanSub-ageMean
  
  ###########################################
  twoScan <- subset(ModelData,ModelData$label=="Time2")
  hist(twoScan$time)
  mean(twoScan$time)
  sd(twoScan$time)
  threeScan <- subset(ModelData,ModelData$label=="Time3")
  hist(threeScan$time)
  mean(threeScan$time)
  sd(threeScan$time)
  
  tempScan <- twoScan
  for (i in 1:nrow(threeScan)) {
    tempScan <- subset(tempScan,tempScan$ID!=threeScan[i,1])
  }
  matchTwoScan <- anti_join(twoScan,tempScan)
  
  tempScan <- threeScan
  for (i in 1:nrow(matchTwoScan)) {
    tempScan <- subset(tempScan,tempScan$ID!=matchTwoScan[i,1])
  }
  matchThreeScan <- anti_join(threeScan,tempScan)
  
  ageTwo2Three <- matchThreeScan$Age1-matchTwoScan$Age1
  hist(ageTwo2Three)
  mean(ageTwo2Three)
  sd(ageTwo2Three)
  ###########################################
  
  
  ModelData$centertime2 <- ModelData$centertime^2
  ModelData[,-c(1:3)] <- scale(ModelData[,-c(1:3)])
  # Meantemp <- aggregate(ModelData[,c(4,5,48:409)],by=list(ModelData$ID), mean)
  # colnames(Meantemp) <- c("ID",paste0("C",colnames(Meantemp)[2:365]))
  # 
  # ModelData <- inner_join(ModelData,Meantemp)
  # ModelData[,c(4,5,48:409)] <- ModelData[,c(4,5,48:409)]-ModelData[,c(410:773)]
  # write.csv(ModelData[,c(1:3,5,6,48,49)],"D:/Rdata/group_SLIM_stats/PlasticityAnalysis/Plasticity.csv",row.names = F)
  # Alpha&Beta
  # Fit Model
  Brain <- colnames(ModelData[,50:409])
  Var <- c("Alpha","E","O","Com","Agen","Ope","Int","time")
  ResultDf <- matrix(0,length(Brain),2*length(Var))
  colnames(ResultDf) <- c(paste(Var,"pvalue",sep = "-"),paste(Var,"tvalue",sep = "-"))
  rownames(ResultDf) <- colnames(ModelData[,50:409])
  
  for (vv in 1:length(Var)) {
    lmstore <- list()
    # centertime
    Inter <- paste0("centertime*",Var)
    Coef <- paste0("centertime:",Var)
    if (vv == 8) Coef[vv] <- "centertime"
    for (aa in 1:360) {
      print(paste(Var[vv],aa))
      ll <- switch(vv,
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Beta*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+E*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+O*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Agen*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Com*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Int*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Ope*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean*centertime+(1|ID)"))
      )
      lmstore[[aa]] <- lmer(ll,data = ModelData)
    }
    if (vv!=2&&vv!=3) {
      PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients[Coef[vv],]))
    }
    if (vv==2) PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients["E:centertime",]))
    if (vv==3) PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients["O:centertime",]))
    # PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients["AgeMean",]))
    # PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients["AgeMean:centertime",]))
    ResultDf[,vv] <- PVector[5,]
    ResultDf[,vv+8] <- PVector[4,]
  }
  ResultDf <- data.frame(ResultDf)
  PlotLabel <- str_replace(rownames(ResultDf),"_ROI_volume","")
  ResultDf$label <- PlotLabel
  ResultDf$label <- str_replace(ResultDf$label,"\\.","-")
  # write.csv(ResultDf,paste(SaveDir,"VolumeResults.csv",sep = ""),row.names = F)
  # FDR-correct
  StatsName <- paste0("Vol",Var,"Stats")
  StatsNameF <- paste0("Vol",Var,"StatsFDR")
  for (vv in 1:length(Var)) {
    cmd <- paste0(StatsName[vv]," <- data.frame(pvalue=p.adjust(ResultDf[,",vv,"],method = 'BH'),tvalue=ResultDf[,",vv+8,"],label=ResultDf$label)")
    eval(parse(text = cmd))
  }
  # Extract sig regions
  for (ss in 1:length(Var)) {
    if (nrow(get(StatsName[ss]))>=1) {
      cmd <- paste0(StatsNameF[ss]," <- subset(",StatsName[ss],",",StatsName[ss],"[,1]<0.05)")
      eval(parse(text = cmd))
    }
  }
  # Can not survival after FDR-correct 
  temp <- data.frame(BehaviorDim=StatsNameF,FDR=1)
  for (ss in 1:length(Var)) {
    if(nrow(get(StatsNameF[ss]))==0) {
      # Print regions can not pass FDR-correct
      temp[ss,2] <- 0
      # Threshold the original p-value at 0.001
      cmd1 <- paste0(StatsNameF[ss]," <- subset(ResultDf",",ResultDf[",ss,"]<0.001)")
      # Select p&t values and labels
      cmd2 <- paste0(StatsNameF[ss]," <- ",StatsNameF[ss],"[,c(",ss,",",ss,"+8,","17)]")
      # Rename colnames
      cmd3 <- paste0("colnames(",StatsNameF[ss],") <- ","colnames(",StatsName[ss],")")
      eval(parse(text = cmd1))
      eval(parse(text = cmd2))
      eval(parse(text = cmd3))
    }
  }
  write.table(temp,paste0(SaveDir,"VolFDR.txt"),quote = F,row.names = F)
}
##########################################################
#                   Thickness Analysis                   #
##########################################################
if (1) {
  setwd("D:\\Rdata\\group_SLIM_stats\\thickness")
  lh_thick_data <- read.csv("lh_brain_data.csv")
  rh_thick_data <- read.csv("rh_brain_data.csv")
  lh_thick_data <- lh_thick_data[order(lh_thick_data[,2]),]
  rh_thick_data <- rh_thick_data[order(rh_thick_data[,1]),]
  
  # data clean
  ModelData <- cbind(behavior_data,lh_thick_data,rh_thick_data)
  ModelData <- ModelData[,-c(1,51:53,234)]
  
  AgeMeanSub <- aggregate(ModelData$Age1,by=list(ModelData$ID),mean)
  colnames(AgeMeanSub)[1:2] <- c("ID","AgeMeanSub")
  ModelData <- inner_join(ModelData,AgeMeanSub,by="ID")
  #time(ij)=Age(ij)-AgeMean(i)
  #ModelData$timeLong <- ModelData$Age1-ModelData$AgeMeanSub
  #AgeMean
  ageMean <- mean(ModelData$Age1)
  #age_mean(i)=AgeMean(i)-AgeMean
  ModelData$age_mean <- ModelData$AgeMeanSub-ageMean
  
  
  ModelData[,-c(1:3)] <- scale(ModelData[,-c(1:3)])
  # Alpha&Beta
  # Fit Model
  Brain <- colnames(ModelData[,50:409])
  Var <- c("Alpha","E","O","Com","Agen","Ope","Int","centertime")
  ResultDf <- matrix(0,length(Brain),2*length(Var))
  colnames(ResultDf) <- c(paste(Var,"pvalue",sep = "-"),paste(Var,"tvalue",sep = "-"))
  rownames(ResultDf) <- colnames(ModelData[,50:409])
  
  for (vv in 1:length(Var)) {
    lmstore <- list()
    Inter <- paste0("centertime*",Var)
    Coef <- paste0("centertime:",Var)
    if (vv == 8) Coef[vv] <- "centertime"
    for (aa in 1:360) {
      print(paste(Var[vv],aa))
      ll <- switch(vv,
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Beta*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+E*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+O*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Agen*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Com*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Int*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Ope*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean*centertime+(1|ID)"))
      )
      lmstore[[aa]] <- lmer(ll,data = ModelData)
    }
    if (vv!=2&&vv!=3) {
      PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients[Coef[vv],]))
    }
    if (vv==2) PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients["E:centertime",]))
    if (vv==3) PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients["O:centertime",]))
    ResultDf[,vv] <- PVector[5,]
    ResultDf[,vv+8] <- PVector[4,]
  }
  ResultDf <- data.frame(ResultDf)
  ResultDf$label <- PlotLabel
  ResultDf$label <- str_replace(ResultDf$label,"\\.","-")
  # write.csv(ResultDf,paste(SaveDir,"ThicknessResults.csv",sep = ""),row.names = F)
  # FDR-correct
  StatsName <- paste0("Thick",Var,"Stats")
  StatsNameF <- paste0("Thick",Var,"StatsFDR")
  for (vv in 1:length(Var)) {
    cmd <- paste0(StatsName[vv]," <- data.frame(pvalue=p.adjust(ResultDf[,",vv,"],method = 'BH'),tvalue=ResultDf[,",vv+8,"],label=ResultDf$label)")
    eval(parse(text = cmd))
  }
  # Extract sig regions
  for (ss in 1:length(Var)) {
    if (nrow(get(StatsName[ss]))>=1) {
      cmd <- paste0(StatsNameF[ss]," <- subset(",StatsName[ss],",",StatsName[ss],"[,1]<0.05)")
      eval(parse(text = cmd))
    }
  }
  # Can not survival after FDR-correct 
  temp <- data.frame(BehaviorDim=StatsNameF,FDR=1)
  for (ss in 1:length(Var)) {
    if(nrow(get(StatsNameF[ss]))==0) {
      # Print regions can not pass FDR-correct
      temp[ss,2] <- 0
      # Threshold the original p-value at 0.001
      cmd1 <- paste0(StatsNameF[ss]," <- subset(ResultDf",",ResultDf[",ss,"]<0.001)")
      # Select p&t values and labels
      cmd2 <- paste0(StatsNameF[ss]," <- ",StatsNameF[ss],"[,c(",ss,",",ss,"+8,","17)]")
      # Rename colnames
      cmd3 <- paste0("colnames(",StatsNameF[ss],") <- ","colnames(",StatsName[ss],")")
      eval(parse(text = cmd1))
      eval(parse(text = cmd2))
      eval(parse(text = cmd3))
    }
  }
  write.table(temp,paste0(SaveDir,"ThickFDR.txt"),quote = F,row.names = F)
}
#####################################################
#                   Area Analysis                   #
#####################################################
if (1) {
  star_time <- Sys.time()
  setwd("D:\\Rdata\\group_SLIM_stats\\area")
  lh_area_data <- read.csv("lh_brain_data.csv")
  rh_area_data <- read.csv("rh_brain_data.csv")
  lh_area_data <- lh_area_data[order(lh_area_data[,2]),]
  rh_area_data <- rh_area_data[order(rh_area_data[,2]),]
  
  # data clean
  ModelData <- cbind(behavior_data,lh_area_data,rh_area_data)
  ModelData <- ModelData[,-c(1,51:54,235:242,423:426)]
  
  AgeMeanSub <- aggregate(ModelData$Age1,by=list(ModelData$ID),mean)
  colnames(AgeMeanSub)[1:2] <- c("ID","AgeMeanSub")
  ModelData <- inner_join(ModelData,AgeMeanSub,by="ID")
  #time(ij)=Age(ij)-AgeMean(i)
  #ModelData$timeLong <- ModelData$Age1-ModelData$AgeMeanSub
  #AgeMean
  ageMean <- mean(ModelData$Age1)
  #age_mean(i)=AgeMean(i)-AgeMean
  ModelData$age_mean <- ModelData$AgeMeanSub-ageMean

  ModelData[,-c(1:3)] <- scale(ModelData[,-c(1:3)],scale = F)
  # Alpha&Beta
  # Fit Model
  Brain <- colnames(ModelData[,50:409])
  Var <- c("Alpha","E","O","Com","Agen","Ope","Int","centertime")
  ResultDf <- matrix(0,length(Brain),2*length(Var))
  colnames(ResultDf) <- c(paste(Var,"pvalue",sep = "-"),paste(Var,"tvalue",sep = "-"))
  rownames(ResultDf) <- colnames(ModelData[,50:409])
  
  for (vv in 1:length(Var)) {
    lmstore <- list()
    Inter <- paste0("centertime*",Var)
    Coef <- paste0("centertime:",Var)
    if (vv == 8) Coef[vv] <- "centertime"
    for (aa in 1:360) {
      print(paste(Var[vv],aa))
      ll <- switch(vv,
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Beta*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+E*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+O*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Agen*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Com*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Int*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean+Ope*centertime+",Inter[vv],"+(1|ID)")),
                   as.formula(paste0(Brain[aa], "~ Sex+age_mean*centertime+(1|ID)"))
      )
      lmstore[[aa]] <- lmer(ll,data = ModelData)
    }
    if (vv!=2&&vv!=3) {
      PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients[Coef[vv],]))
    }
    if (vv==2) PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients["E:centertime",]))
    if (vv==3) PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients["O:centertime",]))
    ResultDf[,vv] <- PVector[5,]
    ResultDf[,vv+8] <- PVector[4,]
  }
  ResultDf <- data.frame(ResultDf)
  ResultDf$label <- PlotLabel
  ResultDf$label <- str_replace(ResultDf$label,"\\.","-")
  # write.csv(ResultDf,paste(SaveDir,"AreaResults.csv",sep = ""))
  # FDR-correct
  StatsName <- paste0("Area",Var,"Stats")
  StatsNameF <- paste0("Area",Var,"StatsFDR")
  for (vv in 1:length(Var)) {
    cmd <- paste0(StatsName[vv]," <- data.frame(pvalue=p.adjust(ResultDf[,",vv,"],method = 'BH'),tvalue=ResultDf[,",vv+8,"],label=ResultDf$label)")
    eval(parse(text = cmd))
  }
  # Extract sig regions
  for (ss in 1:length(Var)) {
    if (nrow(get(StatsName[ss]))>=1) {
      cmd <- paste0(StatsNameF[ss]," <- subset(",StatsName[ss],",",StatsName[ss],"[,1]<0.05)")
      eval(parse(text = cmd))
    }
  }
  # Can not survival after FDR-correct 
  temp <- data.frame(BehaviorDim=StatsNameF,FDR=1)
  for (ss in 1:length(Var)) {
    if(nrow(get(StatsNameF[ss]))==0) {
      # Print regions can not pass FDR-correct
      temp[ss,2] <- 0
      # Threshold the original p-value at 0.001
      cmd1 <- paste0(StatsNameF[ss]," <- subset(ResultDf",",ResultDf[",ss,"]<0.001)")
      # Select p&t values and labels
      cmd2 <- paste0(StatsNameF[ss]," <- ",StatsNameF[ss],"[,c(",ss,",",ss,"+8,","17)]")
      # Rename colnames
      cmd3 <- paste0("colnames(",StatsNameF[ss],") <- ","colnames(",StatsName[ss],")")
      eval(parse(text = cmd1))
      eval(parse(text = cmd2))
      eval(parse(text = cmd3))
    }
  }
  write.table(temp,paste0(SaveDir,"AreaFDR.txt"),quote = F,row.names = F)
  end_time <- Sys.time()
  run_time <- end_time - star_time
  print(run_time)
}

##########################################################
#                   Mediation Analysis                   #
##########################################################
library(mediation)
library(lavaan)
MediationData <- read.csv("D:\\Rdata\\group_SLIM_stats\\PlasticityAnalysis\\MediationTest\\MediationData.csv")
M <- colnames(MediationData)[9:18]
X <- colnames(MediationData)[19]
Y <- colnames(MediationData)[6:8]

ResultDf1 <- matrix(0,length(M),length(Y))
colnames(ResultDf1) <- Y
rownames(ResultDf1) <- M

ResultDf2 <- matrix(0,length(M),length(Y))
colnames(ResultDf2) <- Y
rownames(ResultDf2) <- M
# Fit model
set.seed(773)
MediationData[,-c(1:3)] <- scale(MediationData[,-c(1:3)])
Med <- list()
for (yy in 1:length(Y)) {
  for (mm in 1:length(M)) {
    print(paste(Y[yy],mm))
    l1 <- as.formula(paste0(M[mm],"~AlphaChange+age+Sex"))
    l2 <- as.formula(paste0(Y[yy],"~AlphaChange+",M[mm],"+age+Sex"))
    lm1 <- lm(l1,MediationData)
    lm2 <- lm(l2,MediationData)
    Med[[mm]] <- mediate(lm1,lm2,treat = "AlphaChange",mediator = M[mm],boot = T)
    ResultDf1[mm,yy] <- Med[[mm]][["d0.p"]]
  }
}
ResultDf1 <- data.frame(ResultDf1)
ResultDf1$label <- rownames(ResultDf1)
write.csv(ResultDf1,"D:\\Rdata\\group_SLIM_stats\\PlasticityAnalysis\\MediationTest\\MediationRes.csv",row.names = F)
# fit <- list()
# for (yy in 1:length(Y)) {
#   for (mm in 1:length(M)) {
#     print(paste(Y[yy],mm))
#     myModel <- paste0(paste0(M[mm],"~a*Alpha+age+Sex\n"),
#                       paste0(Y[yy],"~c*Alpha+b*",M[mm],"+age+Sex\n"),
#                       "indirect:=a*b\ndirect:=c\ntotal:=c+(a*b)")
#     fit[[mm]] <- sem(model=myModel,data=MediationData,se="bootstrap")
#     temp <- data.frame(summary(fit[[mm]]))
#     ResultDf2[mm,yy] <- temp[16,9]
#   }
# }

fit1 <- list()
for (mm in 1:length(M)) {
  lm <- as.formula(paste0(M[mm],"~age+Sex+AlphaChange"))
  fit1[[mm]] <- lm(lm,MediationData)
}

AlphaCdf <-  data.frame(t(do.call(cbind,lapply(fit1, function(x) summary(x)$coefficients["AlphaChange",]))))
AlphaCdf$label <- M
PlotSig <- subset(AlphaCdf[AlphaCdf$Pr...t..<0.05,])

SigRegionPvaluePlot <- ggseg3d (ResultDf1,atlas = glasser_3d, 
                              colour = "总分",text = "总分",
                              show.legend = T, 
                              position = "stacked")

############################################################
#                   Subcortical Analysis                   #
############################################################
SubcorticalData <- read.csv(paste0(DataDir,"SubcorticalData.csv"))
Behaviors <- colnames(SubcorticalData)[c(7,39,40,43,44,46:49)]
Regions <- colnames(SubcorticalData)[50:63]
ResultDfSub <- data.frame(matrix(0,length(Regions),2*length(Behaviors)))
rownames(ResultDfSub) <- Regions
colnames(ResultDfSub) <- c(paste0(Behaviors,"-pvalue"),paste0(Behaviors,"-tvalue"))
ResultDfSub$label <- str_replace(rownames(ResultDfSub),"\\.","-")%>%
  str_replace("\\.","-") #why should transform twice?

for (vv in 1:length(Behaviors)) {
  if(vv!=1){
  Inter <- paste0(Behaviors[[vv]],"*centertime")
  Coef <- paste0(Behaviors[[vv]],":centertime")
  fit <- list()
  }
  for (rr in 1:length(Regions)) {
    ll <- switch (vv,
                  as.formula(paste0(Regions[rr], "~ Sex+age+centertime+CSF+ICV+(1|ID)")),
                  as.formula(paste0(Regions[rr], "~ Sex+age+O+CSF+ICV+",Inter,"+(1|ID)")),
                  as.formula(paste0(Regions[rr], "~ Sex+age+E+CSF+ICV+",Inter,"+(1|ID)")),
                  as.formula(paste0(Regions[rr], "~ Sex+age+Int+CSF+ICV+",Inter,"+(1|ID)")),
                  as.formula(paste0(Regions[rr], "~ Sex+age+Ope+CSF+ICV+",Inter,"+(1|ID)")),
                  as.formula(paste0(Regions[rr], "~ Sex+age+Beta+CSF+ICV+",Inter,"+(1|ID)")),
                  as.formula(paste0(Regions[rr], "~ Sex+age+Alpha+CSF+ICV+",Inter,"+(1|ID)")),
                  as.formula(paste0(Regions[rr], "~ Sex+age+Agen+CSF+ICV+",Inter,"+(1|ID)")),
                  as.formula(paste0(Regions[rr], "~ Sex+age+Com+CSF+ICV+",Inter,"+(1|ID)"))
    )
    fit[[rr]] <- lmer(ll,SubcorticalData)
  }
  PVector <- do.call(cbind,lapply(fit, function(x) summary(x)$coefficients["centertime",]))
  if(vv!=1) PVector <- do.call(cbind,lapply(fit, function(x) summary(x)$coefficients[Coef,]))
  ResultDfSub[,vv] <- PVector[5,]
  ResultDfSub[,vv+9] <- PVector[4,]
}
# FDR correct
Stats <- paste0("Stats",Behaviors)
for (vv in 1:length(Behaviors)) {
  cmd <- paste0(Stats[vv]," <- data.frame(pfdr=p.adjust(ResultDfSub[,",vv,"]),tvalue=ResultDfSub[,",vv+9,"],label=ResultDfSub$label)")
  eval(parse(text = cmd))
}
# Subset pfdr < 0.05
for (vv in 1:length(Behaviors)) {
  StatsFdr <- paste0(Stats,"Fdr")
  cmd <- paste0(StatsFdr[vv]," <- subset(",Stats[vv],",",Stats[vv],"$pfdr<0.05)")
  eval(parse(text=cmd))
}

###################################################
#                   FA Analysis                   #
###################################################
setwd("D:\\Rdata\\group_SLIM_stats\\PlasticityAnalysis\\FA分析")
FAdata <- read.csv("FAdata.csv")
ModelData <- FAdata
AgeMean <- aggregate(ModelData$Age1,by=list(ModelData$ID),mean)
colnames(AgeMean)[1:2] <- c("ID","AgeMean")
ModelData <- inner_join(ModelData,AgeMean,by="ID")
ModelData[,-c(1:3)] <- scale(ModelData[,-c(1:3)],scale = F)
Brain <- colnames(ModelData[,50:409])
Var <- c("Alpha","E","O","Com","Agen","Ope","Int","centertime")
ResultDf <- matrix(0,length(Brain),2*length(Var))
colnames(ResultDf) <- c(paste(Var,"pvalue",sep = "-"),paste(Var,"tvalue",sep = "-"))
rownames(ResultDf) <- colnames(ModelData[,50:409])

for (vv in 1:length(Var)) {
  lmstore <- list()
  Inter <- paste0("centertime*",Var)
  Coef <- paste0("centertime:",Var)
  if (vv == 8) Coef[vv] <- "centertime"
  for (aa in 1:360) {
    print(paste(Var[vv],aa))
    ll <- switch(vv,
                 as.formula(paste0(Brain[aa], "~ Sex+AgeMean+Beta:centertime+Beta+",Inter[vv],"+(1|ID)")),
                 as.formula(paste0(Brain[aa], "~ Sex+AgeMean+",Inter[vv],"+(1|ID)")),
                 as.formula(paste0(Brain[aa], "~ Sex+AgeMean+",Inter[vv],"+(1|ID)")),
                 as.formula(paste0(Brain[aa], "~ Sex+AgeMean+Agen:centertime+Agen+",Inter[vv],"+(1|ID)")),
                 as.formula(paste0(Brain[aa], "~ Sex+AgeMean+Com:centertime+Com+",Inter[vv],"+(1|ID)")),
                 as.formula(paste0(Brain[aa], "~ Sex+AgeMean+Int:centertime+Int+",Inter[vv],"+(1|ID)")),
                 as.formula(paste0(Brain[aa], "~ Sex+AgeMean+Ope:centertime+Ope+",Inter[vv],"+(1|ID)")),
                 as.formula(paste0(Brain[aa], "~ Sex+AgeMean*centertime+(1|ID)"))
    )
    lmstore[[aa]] <- lmer(ll,data = ModelData)
  }
  PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients[Coef[vv],]))
  ResultDf[,vv] <- PVector[5,]
  ResultDf[,vv+8] <- PVector[4,]
}
ResultDf <- data.frame(ResultDf)
ResultDf$label <- PlotLabel
ResultDf$label <- str_replace(ResultDf$label,"\\.","-")
write.csv(ResultDf,paste(SaveDir,"FAResults.csv",sep = ""),row.names = F)
# FDR-correct
StatsName <- paste0("FA",Var,"Stats")
StatsNameF <- paste0("FA",Var,"StatsFDR")
for (vv in 1:length(Var)) {
  cmd <- paste0(StatsName[vv]," <- data.frame(pvalue=p.adjust(ResultDf[,",vv,"],method = 'BH'),tvalue=ResultDf[,",vv+8,"],label=ResultDf$label)")
  eval(parse(text = cmd))
}
# Extract sig regions
for (ss in 1:length(Var)) {
  if (nrow(get(StatsName[ss]))>=1) {
    cmd <- paste0(StatsNameF[ss]," <- subset(",StatsName[ss],",",StatsName[ss],"[,1]<0.05)")
    eval(parse(text = cmd))
  }
}

##########################################################
#                   Cell-type Analysis                   #
##########################################################
setwd("D:\\Rdata\\group_SLIM_stats\\PlasticityAnalysis\\EnrichmentAnalysis\\")
VolComFDR <- read.csv("VolComFDR.csv")
ThickComFDR <- read.csv("ThickComFDR.csv")
Celltype <- read.xlsx("celltype.xlsx")
Celltype <- Celltype[,-(1:3)]
Celltype$Class <- str_remove(Celltype$Class,"-")
Cell <- unique(Celltype$Class)[-c(8,9)]

CellAll <- data.frame()
# 
for (cc in 1:length(Cell)) {
  if(length(grep("TRUE",duplicated(CellAll)))!=0) {
    print("ERROR:Duplicate in CellAll")
    break
  }
  print(Cell[cc])
  Temp <- subset(Celltype,Celltype$Class==Cell[cc])
  TempNaRm <- Temp[apply(Temp,c(1,2),function(x) !any(is.na(x)))]
  cmd <- paste0(Cell[cc]," <- data.frame(Celltype=Cell[cc],Gene=unique(TempNaRm)[-1])")
  eval(parse(text=cmd))
  if (ncol(get(Cell[cc]))>2) rm(list = Cell)
  CellAll <- rbind(CellAll,get(Cell[cc]))
  if(cc==length(Cell)){ 
    VolComFDR <- inner_join(VolComFDR,CellAll)
    ThickComFDR <- inner_join(ThickComFDR,CellAll)
  }
}

