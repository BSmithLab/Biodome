
#E Coli Data Analysis
#Written by: Jarrett Eshima

#Before you begin, combine split indole peaks

wd = "D:/Biodome Publication/MarchReviews/Peer Review/Reviews/Replication/AlignedData"
setwd(wd)

tmp = read.csv("EColi_Filtered_DA_Reviews_v3.csv")

EColiContam1 = which(tmp$Mass == 44)
EColiContam2 = which(tmp$Mass == 73)
EColiContam3 = which(tmp$Mass == 207)
EColiContam4 = which(tmp$Mass == 281)
EColiContam = as.numeric(names(table(c(EColiContam1,EColiContam2,EColiContam3,EColiContam4))))
Filt_EColiVOCs = tmp[-EColiContam,]

EColi = Filt_EColiVOCs

rownames(EColi) = make.unique(EColi$Peak)
EColi = EColi[,-1]

#Remove VOCs observed in instrument blanks

inblanks = rep(NA,nrow(EColi))
for(i in 1:nrow(EColi)){
  
  if(sum(is.na(EColi[i,2:7])) < 4){ #analysis timeline left carryover in 2 blanks during replicate 3. If both blank observations occurred during Replicate 3 then the VOC is kept 
    inblanks[i] = i
  }
  
}
inblanks = inblanks[!is.na(inblanks)]

#indole was detected in one of the six blanks due to carryover (< 2*10^3 abundance). We elected to retain it the downstream analysis.

indole = which(rownames(EColi) == "Indole")

inblanks = inblanks[!inblanks %in% indole]

EColi = EColi[-inblanks,]


#write.csv(EColi,"EColi_Reviews_Filtered.csv")

######## MANUALLY REMOVE KNOWN CONTAMINANTS (Siloxanes, Ethanol, single mass peaks) #########
Clean_Coli = read.csv("EColi_Reviews_Filtered_Post.csv") #Supplemental Dataset 4
rownames(Clean_Coli) = Clean_Coli$X
Clean_Coli = Clean_Coli[,-1]

#Remove blank columns
Clean_Coli = Clean_Coli[,-1:-7]

RTref = Clean_Coli

Clean_Coli = Clean_Coli[,-19:-20] #Remove average RTs for stats

################################# Shared volatilome

remind = rep(NA,nrow(Clean_Coli))
for(i in 1:nrow(Clean_Coli)){
  
  tmp = table(is.na(Clean_Coli[i,]))
  
  if(length(tmp)>1){
    if(tmp[[2]]>4){ #at least 70% of all chromatograms
      remind[i] = i
    }
  }
  
  
}

remind = remind[!is.na(remind)]
Final_Coli = Clean_Coli[-remind,]


Final_Coli = log10(Final_Coli)

for(i in 1:nrow(Final_Coli)){
  for(j in 1:ncol(Final_Coli)){
    if(is.na(Final_Coli[i,j])){
      Final_Coli[i,j] = 0
    }
  }
}

SharedVOCs = Final_Coli
SharedVOCs$AvgDim1 = RTref$AvgDim1[-remind]
SharedVOCs$AvgDim2 = RTref$AvgDim2[-remind]

Bacteria = Final_Coli[,1:9]
Broth = Final_Coli[,10:18]

Final_Coli$pval = rep(NA,nrow(Final_Coli))

pval = rep(NA,nrow(Final_Coli))
for(i in 1:nrow(Final_Coli)){
  
  a = as.numeric(Bacteria[i,])
  b = as.numeric(Broth[i,])
  
  res = t.test(a,b,paired = F,alternative = "two.sided")
  
  pval[i] = res$p.value
  
  print(rownames(Bacteria)[i])
  print(pval[i])
  
  Final_Coli$pval[i] = pval[i]
}

adjpval =p.adjust(pval,method = "BH")


##############################################LB VOCs only

Bacteria = Clean_Coli[,1:9]
Broth = Clean_Coli[,10:18]

ECind = rep(NA,nrow(Clean_Coli))
for(i in 1:nrow(Clean_Coli)){
  
  if(length(which(is.na(Broth[i,]))) <= 3 && length(which(is.na(Bacteria[i,]))) >= 6){ #Allow for one spurious and missing observation
    ECind[i] = i
  }
  
}

ECind = ECind[!is.na(ECind)]

LB_VOCs = Clean_Coli[ECind,]
rownames(LB_VOCs) = rownames(Clean_Coli)[ECind]
colnames(LB_VOCs) = colnames(Clean_Coli)


#MScontam = c(3,11,13,16,24)
#LB_VOCs = LB_VOCs[-MScontam,]

LB_VOCs = log10(LB_VOCs)

for(i in 1:nrow(LB_VOCs)){
  for(j in 1:ncol(LB_VOCs)){
    if(is.na(LB_VOCs[i,j])){
      LB_VOCs[i,j] = 0
    }
  }
}

#Convert to Numeric
LB_VOCs2 = matrix(as.numeric(unlist(LB_VOCs)),nrow(LB_VOCs))
colnames(LB_VOCs2) = colnames(LB_VOCs)
rownames(LB_VOCs2) = rownames(LB_VOCs)

Bacteria = LB_VOCs2[,1:9]
Broth = LB_VOCs2[,10:18]

tmp_res_LB = data.frame(matrix(NA,nrow(LB_VOCs2),3))
colnames(tmp_res_LB) = c("VOC","pval","adjp")
tmp_res_LB$VOC = rownames(LB_VOCs2)

pval = rep(NA,nrow(LB_VOCs2))
for(i in 1:nrow(LB_VOCs2)){
  
  a = as.numeric(Bacteria)
  b = as.numeric(Broth)
  
  res = t.test(a,b,paired = F,alternative = "two.sided")
  
  pval[i] = res$p.value
  
  tmp_res_LB$pval[i] = pval[i]
  
}

LB_VOCs3 = data.frame(LB_VOCs2)
LB_VOCs3$AvgDim1 = RTref$AvgDim1[ECind]
LB_VOCs3$AvgDim2 = RTref$AvgDim2[ECind]

##############################################EColi VOCs only

Bacteria = Clean_Coli[,1:9]
Broth = Clean_Coli[,10:18]

ECind = rep(NA,nrow(Clean_Coli))
for(i in 1:nrow(Clean_Coli)){
  
  if(length(which(is.na(Bacteria[i,]))) <= 3 && length(which(is.na(Broth[i,]))) >= 6){ #Allow for one spurious and missing observation
    ECind[i] = i
  }
  
}

ECind = ECind[!is.na(ECind)]

EColi_VOCs = Clean_Coli[ECind,]
rownames(EColi_VOCs) = rownames(Clean_Coli)[ECind]
colnames(EColi_VOCs) = colnames(Clean_Coli)


EColi_VOCs = log10(EColi_VOCs)

for(i in 1:nrow(EColi_VOCs)){
  for(j in 1:ncol(EColi_VOCs)){
    if(is.na(EColi_VOCs[i,j])){
      EColi_VOCs[i,j] = 0
    }
  }
}

#Convert to Numeric
EColi_VOCs2 = matrix(as.numeric(unlist(EColi_VOCs)),nrow(EColi_VOCs))
colnames(EColi_VOCs2) = colnames(EColi_VOCs)
rownames(EColi_VOCs2) = rownames(EColi_VOCs)

Bacteria = EColi_VOCs2[,1:9]
Broth = EColi_VOCs2[,10:18]

tmp_res = data.frame(matrix(NA,nrow(EColi_VOCs2),3))
colnames(tmp_res) = c("VOC","pval","adjp")
tmp_res$VOC = rownames(EColi_VOCs2)

pval = rep(NA,nrow(EColi_VOCs2))
for(i in 1:nrow(EColi_VOCs2)){
  
  a = as.numeric(Bacteria[i,])
  b = as.numeric(Broth[i,])
  
  res = t.test(a,b,paired = F,alternative = "two.sided")
  
  pval[i] = res$p.value
  
  tmp_res$pval[i] = pval[i]
  
}

EColi_VOCs3 = data.frame(EColi_VOCs2)
EColi_VOCs3$AvgDim1 = RTref$AvgDim1[ECind]
EColi_VOCs3$AvgDim2 = RTref$AvgDim2[ECind]

FinalVOCs = rbind(SharedVOCs,LB_VOCs3,EColi_VOCs3)

#Clean Expression Matrix

write.csv(FinalVOCs,"Supplemental_Table1_PeerReview_LogTransformedAbundance_4-29-24.csv")

################## Combine results

AllRes  = Final_Coli
AllRes2 = data.frame(matrix(NA,nrow(AllRes),ncol=3))
colnames(AllRes2) = colnames(tmp_res)
AllRes2$VOC = rownames(AllRes)
AllRes2$pval = AllRes$pval

AllResfinal = rbind(AllRes2,tmp_res_LB,tmp_res)

AllResfinal$adjp = p.adjust(AllResfinal$pval,method = "BH") #must adjust all p-values at the same time

Dim1 = rep(NA,nrow(AllResfinal))
Dim2 = rep(NA,nrow(AllResfinal))
for(i in 1:nrow(EColi)){
  for(j in 1:nrow(AllResfinal)){
    if(rownames(EColi)[i] == AllResfinal$VOC [j]){
      Dim1[j] = EColi$AvgDim1[i]
      Dim2[j] = EColi$AvgDim2[i]
    }
  }
}

FinalResults = cbind(AllResfinal,Dim1,Dim2)

write.csv(AllResfinal,"EColi_StatResults_PeerReview_4-29-24.csv")


##################### Plot results
Plot_Coli = FinalVOCs
Plot_Bacteria = Plot_Coli[,1:9]
Plot_Broth = Plot_Coli[,10:18]

adjpval = FinalResults$adjp

#Publication plots
for(i in 1:nrow(Plot_Coli)){
  
  a = as.numeric(Plot_Bacteria[i,])
  b = as.numeric(Plot_Broth[i,])
  
  boxplot(a,b,main=paste(rownames(Plot_Bacteria)[i]),col=c("tomato1","skyblue1"),xaxt = "n",ylab="Integrated Unique Mass, PQN Normalized",cex.main=2,cex.lab=1.2,cex.axis=1.3, ylim = c(min(c(a,b))-0.1*min(c(a,b)),max(c(a,b)+0.3*max(c(a,b)))))
  #axis(1,at=c(1,2),labels = c("E. Coli","LB Broth"))
  
  if(max(a) > max(b)){
    legend(2.05,max(a)+0.3*max(a),legend = c("E.Coli","LB Broth"),col=c("tomato1","skyblue1"),pch=15,cex = 1.25)
  }else{
    legend(0.45,max(b)+0.3*max(b),legend = c("E.Coli","LB Broth"),col=c("tomato1","skyblue1"),pch=15,cex=1.25)
  }
  
  segments(1,max(c(a,b))+0.2*max(c(a,b)),2,max(c(a,b))+0.2*max(c(a,b)),lwd=3)
  text(1.5,max(c(a,b))+0.25*max(c(a,b)),labels = paste("Adjusted p-value: ",formatC(adjpval[i],format = "e",digits=2),sep=""))
  
  points(rep(1,length(a)),a,col="darkred",pch=19)
  points(rep(2,length(b)),b,col="darkblue",pch=19)
}

days = seq(1,3)
nrep=3
for(i in 1:nrow(Plot_Coli)){
  a = as.numeric(Plot_Bacteria[i,])
  b = as.numeric(Plot_Broth[i,])
  
  day1 = c(1,4,7)
  day2 = c(2,5,8)
  day3 = c(3,6,9)
  
  Day1_Coli_mean = mean(a[day1])
  Day2_Coli_mean = mean(a[day2])
  Day3_Coli_mean = mean(a[day3])
  
  Day1_Coli_se = sd(a[day1])/sqrt(nrep)
  Day2_Coli_se = sd(a[day2])/sqrt(nrep)
  Day3_Coli_se = sd(a[day3])/sqrt(nrep)
  
  Coli_Avg = c(Day1_Coli_mean,Day2_Coli_mean,Day3_Coli_mean)
  Coli_SE = c(Day1_Coli_se,Day2_Coli_se,Day3_Coli_se)
  
  Day1_Broth_mean = mean(b[day1])
  Day2_Broth_mean = mean(b[day2])
  Day3_Broth_mean = mean(b[day3])
  
  Day1_Broth_se = sd(b[day1])/sqrt(nrep)
  Day2_Broth_se = sd(b[day2])/sqrt(nrep)
  Day3_Broth_se = sd(b[day3])/sqrt(nrep)
  
  Broth_Avg = c(Day1_Broth_mean,Day2_Broth_mean,Day3_Broth_mean)
  Broth_SE = c(Day1_Broth_se,Day2_Broth_se,Day3_Broth_se)
  
  # plot(days,a,col="darkred",pch=19,ylim = c(min(c(a,b))-0.1*min(c(a,b)),max(c(a,b)+0.1*max(c(a,b)))),xlab = "Time (days)",ylab="Normalized Abundance (Unique Mass)",main=paste(rownames(Plot_Bacteria)[i]),cex.main=2,cex.lab=1.2,cex.axis=1.3,cex=1.2)
  # lines(days,a,col="tomato1",lwd=2)
  # lines(days,b,col="skyblue1",lwd=2)
  # points(days,b,col="darkblue",pch=19)
  
  par(mar=c(5.1, 5.5, 4.1, 2.1))
  plot(days,Coli_Avg,pch=20,col="darkred",ylim=c(0,8),cex=2,main=paste(rownames(Plot_Bacteria)[i]),xlab="Time (hours)",xaxt='n',ylab="Log10(Abundance)",cex.main=2,cex.axis=1.5,cex.lab=1.5)
  axis(1,at=c(24,48,72),labels = c("24","48","72"),cex.axis = 1.5)
  points(days,Broth_Avg,pch=20,col="darkblue",cex=2)
  lines(days,Coli_Avg,col="tomato1")
  lines(days,Broth_Avg,col="skyblue1")
  legend(24,8,legend = c("E. Coli","Broth"),col = c("darkred","darkblue"),lty=1,lwd=4,cex=2)
  
  arrows(days, Coli_Avg+Coli_SE, days, Coli_Avg-Coli_SE, angle=90, code=3, length=0.08, col="#ba687e")
  arrows(days, Broth_Avg+Broth_SE, days, Broth_Avg-Broth_SE, angle=90, code=3, length=0.08, col="#88b7bd")
  
}

####### PCA

FullVolatilome = Clean_Coli
FullVolatilome[23,12] = NA #Set 0 abundance to NA to avoid errors below
FullVolatilome = log10(FullVolatilome)
for(i in 1:nrow(FullVolatilome)){
  for(j in 1:ncol(FullVolatilome)){
    if(is.na(FullVolatilome[i,j])){
      FullVolatilome[i,j] = 0
    }
  }
}

FinalVOCs2 = FinalVOCs[,-19:-20]

pca = prcomp(t(FinalVOCs2))
PC1 = c(pca$x[,1])
PC2 = c(pca$x[,2])

#Get axis bounds
min1 = min(PC1)
max1 = max(PC1)
min2 = min(PC2)
max2 = max(PC2)

mycolors = c(rep("tomato1",9),rep("skyblue1",9))

plot(PC1,PC2,xlim=c((min1 + (0.2*min1)),(max1 + (0.6*max1))),ylim=c((min2 + (0.2*min2)),(max2 + (0.6*max2))),main="E. Coli Broth PCA",col=mycolors,pch=20,cex=2,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5)
legend(max1 + (0.1*max1),max2 + (0.55*max2),legend = c("E Coli","LB Broth"),col = c("tomato1","skyblue1"),pch=20,pt.cex = 1.5)


pca = prcomp(t(FullVolatilome))
PC1 = c(pca$x[,1])
PC2 = c(pca$x[,2])

#Get axis bounds
min1 = min(PC1)
max1 = max(PC1)
min2 = min(PC2)
max2 = max(PC2)

mycolors = c(rep("tomato1",9),rep("skyblue1",9))

plot(PC1,PC2,xlim=c((min1 + (0.2*min1)),(max1 + (0.6*max1))),ylim=c((min2 + (0.2*min2)),(max2 + (0.6*max2))),main="E. Coli Broth PCA - Full Volatilome",col=mycolors,pch=20,cex=2,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5)
legend(max1 + (0.2*max1),max2 + (0.55*max2),legend = c("E Coli","LB Broth"),col = c("tomato1","skyblue1"),pch=20,pt.cex = 1.5)


####### E. Coli Functional Groups
EColi_FG = c("Unknown","Alcohol","Hydrocarbon","Ketone","Hydrocarbon","Hydrocarbon","Hydrocarbon","Ketone","Aromatic Hydrocarbon","Ketone","Ketone","Ketone","Unknown","Ketone","Alcohol","Ketone","Unknown","Ketone")

#Colors for the piechart
color1 = "darksalmon"
color2 = "darkred"
piecolorfnc = colorRampPalette(c(color1,color2))

pie(table(EColi_FG)/length(EColi_FG),col=piecolorfnc(length(table(EColi_FG))),main="E. Coli specific functional groups",init.angle = -40)

####### LB Functional Groups
LB_FG = c("Other")

#Colors for the piechart
color1 = "aquamarine2"
color2 = "cadetblue4"
piecolorfnc = colorRampPalette(c(color1,color2))

pie(table(LB_FG)/length(LB_FG),col=piecolorfnc(length(table(LB_FG))),main="LB specific functional groups")

####### Shared Functional Groups
Shared_FG = c("Other","Aromatic Hydrocarbon","Hydrocarbon","Aromatic Hydrocarbon","Aromatic Hydrocarbon","Other","Heteroaromatic","Aromatic Hydrocarbon","Hydrocarbon","Hydrocarbon","Heteroaromatic")
length(Shared_FG)

color1 = "chartreuse"
color2 = "darkgreen"
piecolorfnc = colorRampPalette(c(color1,color2))

pie(table(Shared_FG)/length(Shared_FG),col=piecolorfnc(length(table(Shared_FG))),main="Shared VOCs - functional groups",init.angle = -67)



##Check

count1 = count2 = rep(NA,nrow(Clean_Coli))
for(i in 1:nrow(Clean_Coli)){
  if(sum(is.na(Clean_Coli[i,1:9]))>=3 && sum(is.na(Clean_Coli[i,1:9]))<=6){
    count1[i] = i
  }
  if(sum(is.na(Clean_Coli[i,10:18]))>=3 && sum(is.na(Clean_Coli[i,10:18]))<=6){
    count2[i] = i
  }
}
count1 = count1[!is.na(count1)]
count2 = count2[!is.na(count2)]
length(table(c(count1,count2)))
