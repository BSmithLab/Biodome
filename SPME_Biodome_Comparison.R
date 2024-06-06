#SPME in 20mL vs Biodome w TDU

#Manually filtered peaks < 400 RT dim1
#summed peaks between 2094-2110 for indole

# Read the data in
wd = "D:/Biodome Publication/MarchReviews/Peer Review/Reviews/Replication/AlignedData"

setwd(wd)
SPME = read.csv("SPMEComparison_v4.csv")

SPME2 = SPME[-which(!is.na(SPME$Blank)),]

rownames(SPME2) = make.unique(SPME2$Peak)
SPME3 = SPME2[,-1:-3]

########################################## Pre-Processing and Filtering  #######################################################################

reprobiodome = reprovial = rep(NA,nrow(SPME3))
for(i in 1:nrow(SPME3)){
  
  Biodome = SPME3[i,1:9]
  Vial = SPME3[i,10:18]
  
  count1 = sum(is.na(Biodome))
  count2 = sum(is.na(Vial))
  
  if(count1 <= 3){
    reprobiodome[i] = i
  }
  
  if(count2 <= 3){
    reprovial[i] = i
  }
  
  
}

reprobiodome = reprobiodome[!is.na(reprobiodome)]
reprovial = reprovial[!is.na(reprovial)]

reproduciblevocs = as.numeric(names(table(c(reprobiodome,reprovial))))

ComparisonVOCs = SPME3[reproduciblevocs,]
ComparisonVOCs$Mass = SPME2$Mass[reproduciblevocs]

######################################################

#Remove contaminants

#Cells
CellContam1 = which(ComparisonVOCs$Mass == 44)
CellContam2 = which(ComparisonVOCs$Mass == 73)
CellContam3 = which(ComparisonVOCs$Mass == 207)
CellContam4 = which(ComparisonVOCs$Mass == 281)
CellContam = as.numeric(names(table(c(CellContam1,CellContam2,CellContam3,CellContam4))))
Filt_ComparisonVOCs = ComparisonVOCs[-CellContam,]

#Total Signal Increase
BDsum = sum(Filt_ComparisonVOCs$EColi_Day1_R1,Filt_ComparisonVOCs$EColi_Day2_R1,Filt_ComparisonVOCs$EColi_Day3_R1,Filt_ComparisonVOCs$EColi_Day1_R2,Filt_ComparisonVOCs$EColi_Day2_R2,Filt_ComparisonVOCs$EColi_Day3_R2,Filt_ComparisonVOCs$EColi_Day1_R3,Filt_ComparisonVOCs$EColi_Day2_R3,Filt_ComparisonVOCs$EColi_Day3_R3,na.rm = T)
SPMEsum = sum(Filt_ComparisonVOCs$SPME_Day1_R3,Filt_ComparisonVOCs$SPME_Day1_R2,Filt_ComparisonVOCs$SPME_Day1_R3.1,Filt_ComparisonVOCs$SPME_Day2_R1,Filt_ComparisonVOCs$SPME_Day2_R2,Filt_ComparisonVOCs$SPME_Day2_R3,Filt_ComparisonVOCs$SPME_Day3_R1,Filt_ComparisonVOCs$SPME_Day3_R2,Filt_ComparisonVOCs$SPME_Day3_R3,na.rm=T)
BDsum/SPMEsum


FinalComparison = Filt_ComparisonVOCs[,-ncol(Filt_ComparisonVOCs)]

FinalComparison = log10(FinalComparison[,1:18])
FinalComparison$AvgDim1 = Filt_ComparisonVOCs$AvgDim1
FinalComparison$AvgDim2 = Filt_ComparisonVOCs$AvgDim2

#Impute missing values to 0
for(i in 1:nrow(FinalComparison)){
  for(j in 1:ncol(FinalComparison)){
    if(is.na(FinalComparison[i,j])){
      FinalComparison[i,j] = 0
    }
  }
}

par(mar=c(5.1, 5.5, 4.1, 2.1))
nrep=3
time = c(24,48,72)


index_ecoli = c(which(rownames(FinalComparison) == "2-Nonanone.1"),which(rownames(FinalComparison) == "2-Decanone"),
                which(rownames(FinalComparison) == "2-Dodecanone"),which(rownames(FinalComparison) == "3-Tridecanone.1"),
                which(rownames(FinalComparison) == "3-Octanone"),which(rownames(FinalComparison) == "Heptane, 2,2,4,6,6-pentamethyl-"),
                which(rownames(FinalComparison) == "Benzene, 1-ethyl-3-methyl-"),which(rownames(FinalComparison) == "2-Pentadecanol"),
                which(rownames(FinalComparison) == "Pentane, 2,2,4-trimethyl-"),which(rownames(FinalComparison) == "Octane, 4-methyl-"),
                which(rownames(FinalComparison) == "Octadecane"),which(rownames(FinalComparison) == "2-Dodecanone.3"),
                which(rownames(FinalComparison) == "Indole"),which(rownames(FinalComparison) == "Acetic acid, undec-2-enyl ester"),
                which(rownames(FinalComparison) == "Ethanol.11"),which(rownames(FinalComparison) == "Disulfide, dimethyl.1"),
                which(rownames(FinalComparison) == "Hexane, 2,4-dimethyl-"),which(rownames(FinalComparison) == "Dimethyl trisulfide"),
                which(rownames(FinalComparison) == "Ethylbenzene"),which(rownames(FinalComparison) == "Benzene, 1-ethyl-4-methyl-"),
                which(rownames(FinalComparison) == "Analyte 352"),which(rownames(FinalComparison) == "o-Xylene"),
                which(rownames(FinalComparison) == "Pyrazine, trimethyl-"),which(rownames(FinalComparison) == "Benzene, 1,2,3-trimethyl-"),
                which(rownames(FinalComparison) == "Undecane, 4-methyl-"),which(rownames(FinalComparison) == "Nonane, 2,2,4,4,6,8,8-heptamethyl-"),
                which(rownames(FinalComparison) == "2-Dodecanone.2"))



#FinalComparison[which(rownames(FinalComparison) == "2-Nonanone.1"),]
#which(rownames(FinalComparison) == "Styrene")
PlasticContams = FinalComparison[64:65,]
filtnames = c("2-Nonanone","2-Decanone","2-Undecanone","3-Tridecanone","3-Octanone","2,2,4,6,6-Pentamethylheptane","Aromatic Hydrocarbon 1","Alcohol 1","Hydrocarbon 2","2,3-Dimethylheptane","Unknown 3","Ketone 2","Indole",
              "Unknown 2","Unknown 1","Dimethyl disulfide","Hydrocarbon 1","Dimethyl trisulfide","Ethylbenzene","1-Ethyl-4-methyl-benzene","Ketone 1","o-Xylene","Heteroaromatic 1","1,2,3-Trimethylbenzene","4-Methylundecane","4,6-Dimethyldodecane","2-Tridecanone")

tmp = FinalComparison[index_ecoli,]
rownames(tmp) = filtnames

pval1 = pval2 = pval3 = rep(NA,nrow(tmp))

#Stats
for(i in 1:nrow(tmp)){
  
  Biodome = tmp[i,1:9]
  Vial = tmp[i,10:18]
  
  a = c(Biodome$EColi_Day1_R1,Biodome$EColi_Day1_R2,Biodome$EColi_Day1_R3)
  b = c(Vial$SPME_Day1_R2,Vial$SPME_Day1_R3,Vial$SPME_Day1_R3.1)
  
  resd1 = t.test(a,b,paired = F,alternative = "two.sided")
  
  c = c(Biodome$EColi_Day2_R1,Biodome$EColi_Day2_R2,Biodome$EColi_Day2_R3)
  d = c(Vial$SPME_Day2_R1,Vial$SPME_Day2_R2,Vial$SPME_Day2_R3)
  
  resd2 = t.test(c,d,paired = F,alternative = "two.sided")
  
  e = c(Biodome$EColi_Day3_R1,Biodome$EColi_Day3_R2,Biodome$EColi_Day3_R3)
  f = c(Vial$SPME_Day3_R1,Vial$SPME_Day3_R2,Vial$SPME_Day3_R3)
  
  resd3 = t.test(e,f,paired = F,alternative = "two.sided")
  
  pval1[i] = resd1$p.value
  pval2[i] = resd2$p.value
  pval3[i] = resd3$p.value
  
}

adjpval = p.adjust(c(pval1,pval2,pval3),method = "BH")

adjpvalmat = data.frame(matrix(adjpval,nrow(tmp),ncol=3))
colnames(adjpvalmat) = c("Day1","Day2","Day3")
rownames(adjpvalmat) = rownames(tmp)


for(i in 1:length(index_ecoli)){
  
  Biodome = tmp[i,1:9]
  Vial = tmp[i,10:18]
  
  Biodome_avg1 = mean(c(Biodome$EColi_Day1_R1,Biodome$EColi_Day1_R2,Biodome$EColi_Day1_R3))
  Biodome_se1 = sd(c(Biodome$EColi_Day1_R1,Biodome$EColi_Day1_R2,Biodome$EColi_Day1_R3))/sqrt(nrep)
  
  Biodome_avg2 = mean(c(Biodome$EColi_Day2_R1,Biodome$EColi_Day2_R2,Biodome$EColi_Day2_R3))
  Biodome_se2 = sd(c(Biodome$EColi_Day2_R1,Biodome$EColi_Day2_R2,Biodome$EColi_Day2_R3))/sqrt(nrep)
  
  Biodome_avg3 = mean(c(Biodome$EColi_Day3_R1,Biodome$EColi_Day3_R2,Biodome$EColi_Day3_R3))
  Biodome_se3 = sd(c(Biodome$EColi_Day3_R1,Biodome$EColi_Day3_R2,Biodome$EColi_Day3_R3))/sqrt(nrep)
  
  Bavgs = c(Biodome_avg1,Biodome_avg2,Biodome_avg3)
  Bses = c(Biodome_se1,Biodome_se2,Biodome_se3)
  
  
  Vial_avg1 = mean(c(Vial$SPME_Day1_R1,Vial$SPME_Day1_R2,Vial$SPME_Day1_R3))
  Vial_se1 = sd(c(Vial$SPME_Day1_R1,Vial$SPME_Day1_R2,Vial$SPME_Day1_R3))/sqrt(nrep)
  
  Vial_avg2 = mean(c(Vial$SPME_Day2_R1,Vial$SPME_Day2_R2,Vial$SPME_Day2_R3))
  Vial_se2 = sd(c(Vial$SPME_Day2_R1,Vial$SPME_Day2_R2,Vial$SPME_Day2_R3))/sqrt(nrep)
  
  Vial_avg3 = mean(c(Vial$SPME_Day3_R1,Vial$SPME_Day3_R2,Vial$SPME_Day3_R3))
  Vial_se3 = sd(c(Vial$SPME_Day3_R1,Vial$SPME_Day3_R2,Vial$SPME_Day3_R3))/sqrt(nrep)
  
  Vavgs = c(Vial_avg1,Vial_avg2,Vial_avg3)
  Vses = c(Vial_se1,Vial_se2,Vial_se3)
  
  par(mar=c(5.1, 5.5, 4.1, 2.1))
  plot(time,Bavgs,pch=20,col="darkred",ylim=c(0,8),cex=2,main=paste(rownames(Biodome)),xlab="Time (hours)",xaxt='n',ylab="Log10(Abundance)",cex.main=2,cex.axis=1.5,cex.lab=1.5)
  axis(1,at=c(24,48,72),labels = c("24","48","72"),cex.axis = 1.5)
  points(time,Vavgs,pch=20,col="#05557a",cex=2)
  lines(time,Bavgs,col="darkred")
  lines(time,Vavgs,col="#05557a")
  
  if(rownames(Biodome) == "Indole"){
    legend(24,4,legend = c("Biodome","Vial"),col = c("darkred","#05557a"),lty=1,lwd=4,cex=2)
  }else{
    legend(24,8,legend = c("Biodome","Vial"),col = c("darkred","#05557a"),lty=1,lwd=4,cex=2)
  }
  
  arrows(time, Bavgs+Bses, time, Bavgs-Bses, angle=90, code=3, length=0.08, col="#ba687e")
  arrows(time, Vavgs+Vses, time, Vavgs-Vses, angle=90, code=3, length=0.08, col="#88b7bd")
}
