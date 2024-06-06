

wd = "D:/Biodome Publication/MarchReviews/Peer Review/Reviews/Replication/AlignedData"
setwd(wd)
SKOV = read.csv("Figure3_Clean_400sFiltered_AbundanceOnly.csv")

#Confirmed by mean first- and second-dimension retention time
#Unknown 11 = n-Heptadecanol-1
#Unknown 8 = Analyte 1216 (Normal SKOV) = Formamide, N,N-dibutyl- (RNAi)
#Succinic acid diisopropyl ester = Succinic acid diisopropyl ester
#Naphthalene = Naphthalene
#Acetophenone =Acetophenone
#Unknown 7 = Analyte 985 (normal SKOV) = Analyte 799 (RNAi)
#Unknown 5 = Undecane, 5,7-dimethyl-
#Phenol = Phenol
#1-Hexanol, 2-ethyl- = 1-Hexanol, 2-ethyl-
#Aniline = Aniline
#Dimethyl trisulfide = Dimethyl trisulfide
#Unknown 4 = Vanillin, TBDMS derivative
#Furfural = Furfural
#Acetic acid, butyl ester (n-Butyl acetate) = Butanal, 3-methyl-
#Pyrrole = Pyrrole


#Significant VOCs (unique to SKOV3 relative to instrument and media controls)
SKOVvocs = c("Pyrrole","Butanal, 3-methyl-","Furfural","Aniline","Acetophenone",
             "Naphthalene","Succinic acid diisopropyl ester","Dimethyl trisulfide","Phenol","1-Hexanol, 2-ethyl-",
             "Vanillin, TBDMS derivative","Undecane, 5,7-dimethyl-","n-Heptadecanol-1","Formamide, N,N-dibutyl-","Analyte 799") #Unknowns - not to be used as peak name 


PlotSKOV = SKOV[SKOV$Peak %in% SKOVvocs,]

#First Butanal, 3-methyl- is not the same peak as Table 1 n-Butyl acetate
PlotSKOV = PlotSKOV[-which(rownames(PlotSKOV) == 134),]

#First Undecane, 5,7-dimethyl- is not the same peak as Table 1 Unknown 5
PlotSKOV = PlotSKOV[-which(rownames(PlotSKOV) == 298),]

#Second 1-Hexanol, 2-ethyl- is not the same peak as Table 1 1-Hexanol, 2-ethyl-
PlotSKOV = PlotSKOV[-which(rownames(PlotSKOV) == 312),]


#Split Furfural peak
CombinedFurfural = PlotSKOV[which(PlotSKOV$Peak == "Furfural"),]
CombinedFurfural = CombinedFurfural[,-1:-2]

CombinedFurfural2 = as.numeric(colSums(CombinedFurfural,na.rm = T))
CombinedFurfural2[length(CombinedFurfural2)-1] = mean(CombinedFurfural$AvgDim1,na.rm = T)
CombinedFurfural2[length(CombinedFurfural2)] = mean(CombinedFurfural$AvgDim2,na.rm = T)

PlotSKOV[which(rownames(PlotSKOV) == 197),3:ncol(PlotSKOV)] = CombinedFurfural2
PlotSKOV = PlotSKOV[-which(rownames(PlotSKOV) == 204),]

#Split Phenol peak
CombinedPhenol = PlotSKOV[which(PlotSKOV$Peak == "Phenol"),]
CombinedPhenol = CombinedPhenol[,-1:-2]

CombinedPhenol2 = as.numeric(colSums(CombinedPhenol,na.rm = T))
CombinedPhenol2[length(CombinedPhenol2)-1] = mean(CombinedPhenol$AvgDim1,na.rm = T)
CombinedPhenol2[length(CombinedPhenol2)] = mean(CombinedPhenol$AvgDim2,na.rm = T)

PlotSKOV[which(rownames(PlotSKOV) == 295),3:ncol(PlotSKOV)] = CombinedPhenol2
PlotSKOV = PlotSKOV[-which(rownames(PlotSKOV) == 300),]

#Convert 0 from split peaks to NA
for(i in 1:nrow(PlotSKOV)){
  for(j in 1:ncol(PlotSKOV)){
    if(!is.na(PlotSKOV[i,j]) && PlotSKOV[i,j] == 0){
      PlotSKOV[i,j] = NA
    }
  }
}

rownames(PlotSKOV) = PlotSKOV$Peak
PlotSKOV = PlotSKOV[,-1:-7] #Remove blanks too

PlotSKOV2 = PlotSKOV

#log10 transform

tmpSKOV = log10(PlotSKOV[,1:24])

PlotSKOV[,1:24] = tmpSKOV


#Convert missing values to 0
for(i in 1:nrow(PlotSKOV)){
  for(j in 1:ncol(PlotSKOV)){
    if(is.na(PlotSKOV[i,j])){
      PlotSKOV[i,j] = 0
    }
  }
}


#Write supplemental dataset 5
wd = "D:/Biodome Publication/MarchReviews/Peer Review/Reviews/Replication/RNAi"
setwd(wd)
#write.csv(PlotSKOV,"Supplemental_Dataset_5_PeerReview_RNAi_5-27-24.csv")

#############Stats and Plots

Plot_SKOV = PlotSKOV[,1:12]
Plot_RNAi = PlotSKOV[,13:24]


#Get adjusted p-values 
FinalResults = data.frame(matrix(NA,nrow(PlotSKOV),ncol=3))
colnames(FinalResults) = c("VOC","pval","adjp")
FinalResults$VOC = rownames(PlotSKOV)

for(i in 1:nrow(PlotSKOV)){
  
  a = as.numeric(Plot_SKOV[i,])
  b = as.numeric(Plot_RNAi[i,])
  
  res = t.test(a,b,paired = F,alternative = "two.sided")
  
  FinalResults$pval[i] = res$p.value
  
}

FinalResults$adjp = p.adjust(FinalResults$pval,method = "BH")


##########################################  Publication plots  #######################################################################

library(beeswarm)
for(i in 1:nrow(PlotSKOV)){
  
  a = as.numeric(Plot_SKOV[i,])
  b = as.numeric(Plot_RNAi[i,])
  
  boxplot(a,b,main=paste(rownames(PlotSKOV)[i]),col=c("#BFBFBF","#C00000"),xaxt = "n",ylab="Integrated Unique Mass, Log10 Normalized",cex.main=2,cex.lab=1.5,cex.axis=1.5, ylim = c(0,8),outline = F)
  #Variable y
  #boxplot(a,b,c,main=paste(rownames(PlotSKOV)[i]),col=c("#d1641b","#8ebade"),xaxt = "n",ylab="Integrated Unique Mass, Log10 Normalized",cex.main=2,cex.lab=1.5,cex.axis=1.5, ylim = c(min(c(a,b,c))-0.1*min(c(a,b,c)),max(c(a,b,c)+0.2*max(c(a,b,c)))))
  
  if(max(a) > max(b)){
    legend(2.7,8,legend = c("SK-OV-3","RPMI1640 Media"),col=c("#BFBFBF","#C00000"),pch=15,cex = 1.25)
  }else{
    legend(0.5,8,legend = c("SK-OV-3","RPMI1640 Media"),col=c("#BFBFBF","#C00000"),pch=15,cex=1.25)
  }
  
  #segments(1,max(c(a,b))+0.2*max(c(a,b)),2,max(c(a,b))+0.2*max(c(a,b)),lwd=3)
  #text(1.5,max(c(a,b))+0.25*max(c(a,b)),labels = paste("Adjusted p-value: ",round(adjpval[i],3),sep=""))
  
  beeswarm(a,method = "center",col="black",pch=19,cex=1.2,add=T,at = 1)
  beeswarm(b,method = "center",col="black",pch=19,cex=1.2,add=T,at = 2)
}


days = seq(1,4)
nrep=3
daysig = data.frame(matrix(NA,nrow(PlotSKOV),ncol=5))
colnames(daysig) = c("Peak","Day1p","Day2p","Day3p","Day4p")
daysig$Peak = rownames(PlotSKOV)

for(i in 1:nrow(PlotSKOV)){
  
  Day1_SKOV_mean = mean(c(Plot_SKOV$Light_Day1_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day1_R2.1[i],Plot_SKOV$Light_SKOV3_Day1_R3.1[i]))
  Day2_SKOV_mean = mean(c(Plot_SKOV$Light_Day2_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day2_R2.1[i],Plot_SKOV$Light_SKOV3_Day2_R3.1[i]))
  Day3_SKOV_mean = mean(c(Plot_SKOV$Light_Day3_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day3_R2.1[i],Plot_SKOV$Light_SKOV3_Day3_R3.1[i]))
  Day4_SKOV_mean = mean(c(Plot_SKOV$Light_Day4_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day4_R2.1[i],Plot_SKOV$Light_SKOV3_Day4_R3.1[i]))
  
  Day1_SKOV_se = sd(c(Plot_SKOV$Light_Day1_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day1_R2.1[i],Plot_SKOV$Light_SKOV3_Day1_R3.1[i]))/sqrt(nrep)
  Day2_SKOV_se = sd(c(Plot_SKOV$Light_Day2_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day2_R2.1[i],Plot_SKOV$Light_SKOV3_Day2_R3.1[i]))/sqrt(nrep)
  Day3_SKOV_se = sd(c(Plot_SKOV$Light_Day3_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day3_R2.1[i],Plot_SKOV$Light_SKOV3_Day3_R3.1[i]))/sqrt(nrep)
  Day4_SKOV_se = sd(c(Plot_SKOV$Light_Day4_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day4_R2.1[i],Plot_SKOV$Light_SKOV3_Day4_R3.1[i]))/sqrt(nrep)
  
  SKOV_Avg = c(Day1_SKOV_mean,Day2_SKOV_mean,Day3_SKOV_mean,Day4_SKOV_mean)
  SKOV_SE = c(Day1_SKOV_se,Day2_SKOV_se,Day3_SKOV_se,Day4_SKOV_se)
  
  
  Day1_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day1_3.31.23.1[i],Plot_RNAi$RNAi_Day1_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day1_2.9.24.1[i]))
  Day2_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day2_3.31.23.1[i],Plot_RNAi$RNAi_Day2_R2_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day2_2.9.24.1[i]))
  Day3_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day3_3.31.23.1[i],Plot_RNAi$RNAi_Day3_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day3_2.9.24.1[i]))
  Day4_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day4_3.31.23.1[i],Plot_RNAi$RNAi_Day4_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day4_2.9.24.1[i]))
  
  Day1_RNAi_se = sd(c(Plot_RNAi$CPT2_Rep2_Day1_3.31.23.1[i],Plot_RNAi$RNAi_Day1_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day1_2.9.24.1[i]))/sqrt(nrep)
  Day2_RNAi_se = sd(c(Plot_RNAi$CPT2_Rep2_Day2_3.31.23.1[i],Plot_RNAi$RNAi_Day2_R2_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day2_2.9.24.1[i]))/sqrt(nrep)
  Day3_RNAi_se = sd(c(Plot_RNAi$CPT2_Rep2_Day3_3.31.23.1[i],Plot_RNAi$RNAi_Day3_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day3_2.9.24.1[i]))/sqrt(nrep)
  Day4_RNAi_se = sd(c(Plot_RNAi$CPT2_Rep2_Day4_3.31.23.1[i],Plot_RNAi$RNAi_Day4_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day4_2.9.24.1[i]))/sqrt(nrep)
  
  RNAi_Avg = c(Day1_RNAi_mean,Day2_RNAi_mean,Day3_RNAi_mean,Day4_RNAi_mean)
  RNAi_SE = c(Day1_RNAi_se,Day2_RNAi_se,Day3_RNAi_se,Day4_RNAi_se)
  
  daysig$Day1p[i] = t.test(c(Plot_SKOV$Light_Day1_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day1_R2.1[i],Plot_SKOV$Light_SKOV3_Day1_R3.1[i]),c(Plot_RNAi$CPT2_Rep2_Day1_3.31.23.1[i],Plot_RNAi$RNAi_Day1_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day1_2.9.24.1[i]),alternative = "greater")$p.value
  daysig$Day2p[i] = t.test(c(Plot_SKOV$Light_Day2_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day2_R2.1[i],Plot_SKOV$Light_SKOV3_Day2_R3.1[i]),c(Plot_RNAi$CPT2_Rep2_Day2_3.31.23.1[i],Plot_RNAi$RNAi_Day2_R2_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day2_2.9.24.1[i]),alternative = "greater")$p.value
  daysig$Day3p[i] = t.test(c(Plot_SKOV$Light_Day3_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day3_R2.1[i],Plot_SKOV$Light_SKOV3_Day3_R3.1[i]),c(Plot_RNAi$CPT2_Rep2_Day3_3.31.23.1[i],Plot_RNAi$RNAi_Day3_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day3_2.9.24.1[i]),alternative = "greater")$p.value
  daysig$Day4p[i] = t.test(c(Plot_SKOV$Light_Day4_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day4_R2.1[i],Plot_SKOV$Light_SKOV3_Day4_R3.1[i]),c(Plot_RNAi$CPT2_Rep2_Day4_3.31.23.1[i],Plot_RNAi$RNAi_Day4_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day4_2.9.24.1[i]),alternative = "greater")$p.value
  
  
  plot(days,SKOV_Avg,col="black",pch=19,ylim = c(0,8),xlab = "Time (days)",ylab="Log10 Normalized Abundance (Unique Mass)",main=paste(rownames(PlotSKOV)[i]),cex.main=2,cex.lab=1.2,cex.axis=1.3,cex=1.2)
  #Variable y
  #plot(days,a,col="#961403",pch=19,ylim = c(min(c(a,b))-0.1*min(c(a,b)),max(c(a,b)+0.1*max(c(a,b)))),xlab = "Time (days)",ylab="Log10 Normalized Abundance (Unique Mass)",main=paste(rownames(PlotSKOV)[i]),cex.main=2,cex.lab=1.2,cex.axis=1.3,cex=1.2)
  lines(days,SKOV_Avg,col="#BFBFBF",lwd=2)
  lines(days,RNAi_Avg,col="#C00000",lwd=2)
  points(days,RNAi_Avg,col="black",pch=19)
  
  legend(3,8,legend = c("Control SK-OV-3","RNAi SK-OV-3"),col = c("#BFBFBF","#C00000"),lty=1,lwd=4,cex=1.2)
  
  
  arrows(days, SKOV_Avg+SKOV_SE, days, SKOV_Avg-SKOV_SE, angle=90, code=3, length=0.08, col="#BFBFBF")
  arrows(days, RNAi_Avg+RNAi_SE, days, RNAi_Avg-RNAi_SE, angle=90, code=3, length=0.08, col="#C00000")
  
}

rownames(daysig) = daysig$Peak
daysig = daysig[,-1]

adjdaysig = data.frame(matrix(p.adjust(as.numeric(unlist(daysig)),method = "BH"),nrow(daysig),ncol(daysig)))
colnames(adjdaysig) = colnames(daysig)
rownames(adjdaysig) = rownames(daysig)



################################## Average % decreases

#Convert missing values to 0
for(i in 1:nrow(PlotSKOV2)){
  for(j in 1:ncol(PlotSKOV2)){
    if(is.na(PlotSKOV2[i,j])){
      PlotSKOV2[i,j] = 0
    }
  }
}
Plot_SKOV2 = PlotSKOV2[,1:12]
Plot_RNAi2 = PlotSKOV2[,13:24]
ctr = mean(as.numeric(Plot_SKOV2[which(rownames(Plot_SKOV2) == "1-Hexanol, 2-ethyl-"),]))
rnai = mean(as.numeric(Plot_RNAi2[which(rownames(Plot_RNAi2) == "1-Hexanol, 2-ethyl-"),]))
#Overall Average
1-rnai/ctr #% decrease


i = which(rownames(Plot_SKOV2) == "1-Hexanol, 2-ethyl-")

Plot_SKOV = Plot_SKOV2
Day1_SKOV_mean = mean(c(Plot_SKOV$Light_Day1_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day1_R2.1[i],Plot_SKOV$Light_SKOV3_Day1_R3.1[i]))
Day2_SKOV_mean = mean(c(Plot_SKOV$Light_Day2_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day2_R2.1[i],Plot_SKOV$Light_SKOV3_Day2_R3.1[i]))
Day3_SKOV_mean = mean(c(Plot_SKOV$Light_Day3_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day3_R2.1[i],Plot_SKOV$Light_SKOV3_Day3_R3.1[i]))
Day4_SKOV_mean = mean(c(Plot_SKOV$Light_Day4_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day4_R2.1[i],Plot_SKOV$Light_SKOV3_Day4_R3.1[i]))

Plot_RNAi = Plot_RNAi2
Day1_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day1_3.31.23.1[i],Plot_RNAi$RNAi_Day1_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day1_2.9.24.1[i]))
Day2_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day2_3.31.23.1[i],Plot_RNAi$RNAi_Day2_R2_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day2_2.9.24.1[i]))
Day3_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day3_3.31.23.1[i],Plot_RNAi$RNAi_Day3_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day3_2.9.24.1[i]))
Day4_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day4_3.31.23.1[i],Plot_RNAi$RNAi_Day4_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day4_2.9.24.1[i]))

#Day1 Average
Day1_RNAi_mean/Day1_SKOV_mean 
1-Day1_RNAi_mean/Day1_SKOV_mean #% decrease
#Day2 Average
Day2_RNAi_mean/Day2_SKOV_mean
1-Day2_SKOV_mean/Day2_RNAi_mean #% decrease
#Day2 Average
1-Day3_RNAi_mean/Day3_SKOV_mean #% decrease
#Day2 Average
Day4_SKOV_mean/Day4_RNAi_mean
1-Day4_SKOV_mean/Day4_RNAi_mean #% decrease






ctr = mean(as.numeric(Plot_SKOV2[which(rownames(Plot_SKOV2) == "Furfural"),]))
rnai = mean(as.numeric(Plot_RNAi2[which(rownames(Plot_RNAi2) == "Furfural"),]))
#Overall Average
1-rnai/ctr #% decrease


i = which(rownames(Plot_SKOV2) == "Furfural")

Plot_SKOV = Plot_SKOV2
Day1_SKOV_mean = mean(c(Plot_SKOV$Light_Day1_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day1_R2.1[i],Plot_SKOV$Light_SKOV3_Day1_R3.1[i]))
Day2_SKOV_mean = mean(c(Plot_SKOV$Light_Day2_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day2_R2.1[i],Plot_SKOV$Light_SKOV3_Day2_R3.1[i]))
Day3_SKOV_mean = mean(c(Plot_SKOV$Light_Day3_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day3_R2.1[i],Plot_SKOV$Light_SKOV3_Day3_R3.1[i]))
Day4_SKOV_mean = mean(c(Plot_SKOV$Light_Day4_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day4_R2.1[i],Plot_SKOV$Light_SKOV3_Day4_R3.1[i]))

Plot_RNAi = Plot_RNAi2
Day1_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day1_3.31.23.1[i],Plot_RNAi$RNAi_Day1_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day1_2.9.24.1[i]))
Day2_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day2_3.31.23.1[i],Plot_RNAi$RNAi_Day2_R2_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day2_2.9.24.1[i]))
Day3_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day3_3.31.23.1[i],Plot_RNAi$RNAi_Day3_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day3_2.9.24.1[i]))
Day4_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day4_3.31.23.1[i],Plot_RNAi$RNAi_Day4_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day4_2.9.24.1[i]))

#Day1 Average
Day1_RNAi_mean/Day1_SKOV_mean 
1-Day1_RNAi_mean/Day1_SKOV_mean #% decrease
#Day2 Average
Day2_RNAi_mean/Day2_SKOV_mean
#1-Day2_SKOV_mean/Day2_RNAi_mean #% decrease
#Day2 Average
1-Day3_RNAi_mean/Day3_SKOV_mean #% decrease
#Day2 Average
1-Day4_RNAi_mean/Day4_SKOV_mean #% decrease






ctr = mean(as.numeric(Plot_SKOV2[which(rownames(Plot_SKOV2) == "Pyrrole"),]))
rnai = mean(as.numeric(Plot_RNAi2[which(rownames(Plot_RNAi2) == "Pyrrole"),]))
#Overall Average
rnai / ctr
#1-rnai/ctr #% decrease


i = which(rownames(Plot_SKOV2) == "Pyrrole")

Plot_SKOV = Plot_SKOV2
Day1_SKOV_mean = mean(c(Plot_SKOV$Light_Day1_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day1_R2.1[i],Plot_SKOV$Light_SKOV3_Day1_R3.1[i]))
Day2_SKOV_mean = mean(c(Plot_SKOV$Light_Day2_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day2_R2.1[i],Plot_SKOV$Light_SKOV3_Day2_R3.1[i]))
Day3_SKOV_mean = mean(c(Plot_SKOV$Light_Day3_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day3_R2.1[i],Plot_SKOV$Light_SKOV3_Day3_R3.1[i]))
Day4_SKOV_mean = mean(c(Plot_SKOV$Light_Day4_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day4_R2.1[i],Plot_SKOV$Light_SKOV3_Day4_R3.1[i]))

Plot_RNAi = Plot_RNAi2
Day1_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day1_3.31.23.1[i],Plot_RNAi$RNAi_Day1_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day1_2.9.24.1[i]))
Day2_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day2_3.31.23.1[i],Plot_RNAi$RNAi_Day2_R2_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day2_2.9.24.1[i]))
Day3_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day3_3.31.23.1[i],Plot_RNAi$RNAi_Day3_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day3_2.9.24.1[i]))
Day4_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day4_3.31.23.1[i],Plot_RNAi$RNAi_Day4_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day4_2.9.24.1[i]))

#Day1 Average
Day1_RNAi_mean/Day1_SKOV_mean 
1-Day1_RNAi_mean/Day1_SKOV_mean #% decrease
#Day2 Average
Day2_RNAi_mean/Day2_SKOV_mean
#1-Day2_SKOV_mean/Day2_RNAi_mean #% decrease
#Day2 Average
1-Day3_RNAi_mean/Day3_SKOV_mean #% decrease
#Day2 Average
1-Day4_RNAi_mean/Day4_SKOV_mean #% decrease








ctr = mean(as.numeric(Plot_SKOV2[which(rownames(Plot_SKOV2) == "Formamide, N,N-dibutyl-"),]))
rnai = mean(as.numeric(Plot_RNAi2[which(rownames(Plot_RNAi2) == "Formamide, N,N-dibutyl-"),]))
#Overall Average
rnai / ctr
1-rnai/ctr #% decrease


i = which(rownames(Plot_SKOV2) == "Formamide, N,N-dibutyl-")

Plot_SKOV = Plot_SKOV2
Day1_SKOV_mean = mean(c(Plot_SKOV$Light_Day1_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day1_R2.1[i],Plot_SKOV$Light_SKOV3_Day1_R3.1[i]))
Day2_SKOV_mean = mean(c(Plot_SKOV$Light_Day2_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day2_R2.1[i],Plot_SKOV$Light_SKOV3_Day2_R3.1[i]))
Day3_SKOV_mean = mean(c(Plot_SKOV$Light_Day3_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day3_R2.1[i],Plot_SKOV$Light_SKOV3_Day3_R3.1[i]))
Day4_SKOV_mean = mean(c(Plot_SKOV$Light_Day4_3.3.23.1[i],Plot_SKOV$Light_SKOV3_Day4_R2.1[i],Plot_SKOV$Light_SKOV3_Day4_R3.1[i]))

Plot_RNAi = Plot_RNAi2
Day1_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day1_3.31.23.1[i],Plot_RNAi$RNAi_Day1_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day1_2.9.24.1[i]))
Day2_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day2_3.31.23.1[i],Plot_RNAi$RNAi_Day2_R2_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day2_2.9.24.1[i]))
Day3_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day3_3.31.23.1[i],Plot_RNAi$RNAi_Day3_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day3_2.9.24.1[i]))
Day4_RNAi_mean = mean(c(Plot_RNAi$CPT2_Rep2_Day4_3.31.23.1[i],Plot_RNAi$RNAi_Day4_R3_1.14.24.1[i],Plot_RNAi$RNAi_R3_Day4_2.9.24.1[i]))

#Day1 Average
Day1_RNAi_mean/Day1_SKOV_mean 
1-Day1_RNAi_mean/Day1_SKOV_mean #% decrease
#Day2 Average
Day2_RNAi_mean/Day2_SKOV_mean
#1-Day2_SKOV_mean/Day2_RNAi_mean #% decrease
#Day2 Average
1-Day3_RNAi_mean/Day3_SKOV_mean #% decrease
#Day2 Average
1-Day4_RNAi_mean/Day4_SKOV_mean #% decrease




  