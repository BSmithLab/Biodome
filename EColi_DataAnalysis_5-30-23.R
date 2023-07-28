
#E Coli Data Analysis
#Written by: Jarrett Eshima

wd = "C:/Users/jeshima/Documents/Biodome Manuscript/EColiFigures"
setwd(wd)

EColi = read.csv("EColi_v3_5-30-23_Filtered_DA.csv")


rownames(EColi) = make.unique(EColi$Peak)
EColi = EColi[,-1]

#Remove blank columns (performed manually)
EColi = EColi[,-1:-2]

Clean_Coli = EColi[,-7:-8] #Remove average RTs

################################# Shared volatilome

remind = rep(NA,nrow(Clean_Coli))
for(i in 1:nrow(Clean_Coli)){
  
  tmp = table(is.na(Clean_Coli[i,]))
  
  if(length(tmp)>1){
    if(tmp[[2]]>1){
      remind[i] = i
    }
  }
  
}



remind = remind[!is.na(remind)]
Final_Coli = Clean_Coli[-remind,]

MScontam = c(2,14,54,55)
Final_Coli = Final_Coli[-MScontam,]

Final_Coli = log10(Final_Coli)

for(i in 1:nrow(Final_Coli)){
  for(j in 1:ncol(Final_Coli)){
    if(is.na(Final_Coli[i,j])){
      Final_Coli[i,j] = 0
    }
  }
}

SharedVOCs = Final_Coli

Bacteria = Final_Coli[,1:3]
Broth = Final_Coli[,4:6]

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

Bacteria = Clean_Coli[,1:3]
Broth = Clean_Coli[,4:6]

ECind = rep(NA,nrow(Clean_Coli))
for(i in 1:nrow(Clean_Coli)){
  
  if(length(which(is.na(Broth[i,]))) == 0 && length(which(is.na(Bacteria[i,]))) == 3){
    ECind[i] = i
  }
  
}

ECind = ECind[!is.na(ECind)]

LB_VOCs = Clean_Coli[ECind,]
rownames(LB_VOCs) = rownames(Clean_Coli)[ECind]
colnames(LB_VOCs) = colnames(Clean_Coli)

MScontam = c(3,11,13,16,24)
LB_VOCs = LB_VOCs[-MScontam,]

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

Bacteria = LB_VOCs2[,1:3]
Broth = LB_VOCs2[,4:6]

tmp_res_LB = data.frame(matrix(NA,nrow(LB_VOCs2),3))
colnames(tmp_res_LB) = c("VOC","pval","adjp")
tmp_res_LB$VOC = rownames(LB_VOCs2)

pval = rep(NA,nrow(LB_VOCs2))
for(i in 1:nrow(LB_VOCs2)){
  
  a = as.numeric(Bacteria[i,])
  b = as.numeric(Broth[i,])
  
  res = t.test(a,b,paired = F,alternative = "two.sided")
  
  pval[i] = res$p.value
  
  tmp_res_LB$pval[i] = pval[i]
  
}


##############################################EColi VOCs only

Bacteria = Clean_Coli[,1:3]
Broth = Clean_Coli[,4:6]

ECind = rep(NA,nrow(Clean_Coli))
for(i in 1:nrow(Clean_Coli)){
  
  if(length(which(is.na(Bacteria[i,]))) == 0 && length(which(is.na(Broth[i,]))) == 3){
    ECind[i] = i
  }
  
}

ECind = ECind[!is.na(ECind)]

EColi_VOCs = Clean_Coli[ECind,]
rownames(EColi_VOCs) = rownames(Clean_Coli)[ECind]
colnames(EColi_VOCs) = colnames(Clean_Coli)

MScontam = 1
EColi_VOCs = EColi_VOCs[-MScontam,]

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

Bacteria = EColi_VOCs2[,1:3]
Broth = EColi_VOCs2[,4:6]

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



FinalVOCs = rbind(SharedVOCs,LB_VOCs2,EColi_VOCs2)

#Clean Expression Matrix

write.csv(FinalVOCs,"SupplementalData_Biodome_EColi_LB_LogTransformedAbundance_6-2-23.csv")

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

write.csv(FinalResults,"EColi_Table2Results_v2_6-2-23.csv")


##################### Plot results
PubVOCs = c(23,33,36,66,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,91,92,93,95,96,99,103,105,106) #After manually filtering by ID confidence level
Plot_Coli = FinalVOCs[PubVOCs,]
Plot_Bacteria = Plot_Coli[,1:3]
Plot_Broth = Plot_Coli[,4:6]


#Publication plots


for(i in 1:nrow(Plot_Coli)){
  
  a = as.numeric(Plot_Bacteria[i,])
  b = as.numeric(Plot_Broth[i,])
  
  boxplot(a,b,main=paste(rownames(Plot_Bacteria)[i]),col=c("tomato1","skyblue1"),xaxt = "n",ylab="Integrated Unique Mass, PQN Normalized",cex.main=2,cex.lab=1.2,cex.axis=1.3, ylim = c(min(c(a,b))-0.1*min(c(a,b)),max(c(a,b)+0.3*max(c(a,b)))))
  #axis(1,at=c(1,2),labels = c("E. Coli","LB Broth"))
  
  if(max(a) > max(b)){
    legend(2.15,max(a)+0.3*max(a),legend = c("E.Coli","LB Broth"),col=c("tomato1","skyblue1"),pch=15,cex = 1.25)
  }else{
    legend(0.5,max(b)+0.3*max(b),legend = c("E.Coli","LB Broth"),col=c("tomato1","skyblue1"),pch=15,cex=1.25)
  }
  
  segments(1,max(c(a,b))+0.2*max(c(a,b)),2,max(c(a,b))+0.2*max(c(a,b)),lwd=3)
  text(1.5,max(c(a,b))+0.25*max(c(a,b)),labels = paste("Adjusted p-value: ",round(adjpval[i],3),sep=""))
  
  points(rep(1,length(a)),a,col="darkred",pch=19)
  points(rep(2,length(b)),b,col="darkblue",pch=19)
}


for(i in 1:nrow(Plot_Coli)){
  a = as.numeric(Plot_Bacteria[i,])
  b = as.numeric(Plot_Broth[i,])
  days = seq(1,3)
  plot(days,a,col="darkred",pch=19,ylim = c(min(c(a,b))-0.1*min(c(a,b)),max(c(a,b)+0.1*max(c(a,b)))),xlab = "Time (days)",ylab="Normalized Abundance (Unique Mass)",main=paste(rownames(Plot_Bacteria)[i]),cex.main=2,cex.lab=1.2,cex.axis=1.3,cex=1.2)
  lines(days,a,col="tomato1",lwd=2)
  lines(days,b,col="skyblue1",lwd=2)
  points(days,b,col="darkblue",pch=19)
}

####### PCA

pca = prcomp(t(FinalVOCs))
PC1 = c(pca$x[,1])
PC2 = c(pca$x[,2])

#Get axis bounds
min1 = min(PC1)
max1 = max(PC1)
min2 = min(PC2)
max2 = max(PC2)

mycolors = c(rep("tomato1",3),rep("skyblue1",3))

plot(PC1,PC2,xlim=c((min1 + (0.2*min1)),(max1 + (0.6*max1))),ylim=c((min2 + (0.2*min2)),(max2 + (0.6*max2))),main="E. Coli Broth PCA",col=mycolors,pch=20,cex=2,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5)
legend(max1 + (0.2*max1),max2 + (0.55*max2),legend = c("E Coli","LB Broth"),col = c("tomato1","skyblue1"),pch=20,pt.cex = 1.5)


####### E. Coli Functional Groups
EColi_FG = c("Aldehyde","Alcohol","Hydrocarbon","Unknown","Unknown","Ketone","Unknown","Unknown","Unknown","Functionalized Benzene","Unknown","Ketone","Functionalized Benzene","Unknown","Unknown")

#Colors for the piechart
color1 = "darksalmon"
color2 = "darkred"
piecolorfnc = colorRampPalette(c(color1,color2))

pie(table(EColi_FG)/length(EColi_FG),col=piecolorfnc(length(table(EColi_FG))),main="E. Coli specific functional groups")

####### LB Functional Groups
LB_FG = c("Unknown","Heteroaromatic","Unknown","Unknown","Unknown","Other","Unknown","Unknown","Unknown","Unknown","Ketone","Unknown","Unknown","Functionalized Benzene","Unknown","Unknown","Unknown","Ester","Functionalized Benzene","Unknown","Ketone","Unknown","Unknown","Unknown")

#Colors for the piechart
color1 = "aquamarine2"
color2 = "cadetblue4"
piecolorfnc = colorRampPalette(c(color1,color2))

pie(table(LB_FG)/length(LB_FG),col=piecolorfnc(length(table(LB_FG))),main="LB specific functional groups")

####### Shared Functional Groups
Shared_FG = c(
  "Ketone","Unknown","Unknown","Hydrocarbon","Heteroaromatic","Ketone","Hydrocarbon","Hydrocarbon","Hydrocarbon","Ketone","Hydrocarbon","Unknown","Heteroaromatic","Hydrocarbon","Ketone","Hydrocarbon","Unknown","Unknown","Hydrocarbon","Unknown","Aromatic Hydrocarbon","Unknown","Hydrocarbon",
  "Unknown","Other","Hydrocarbon","Heteroaromatic","Aromatic Hydrocarbon","Heteroaromatic","Hydrocarbon","Hydrocarbon","Hydrocarbon","Unknown","Hydrocarbon","Unknown","Heteroaromatic","Hydrocarbon",
  "Aromatic Hydrocarbon","Hydrocarbon","Hydrocarbon","Azole","Unknown","Heteroaromatic","Hydrocarbon","Heteroaromatic","Hydrocarbon","Other","Unknown","Alcohol","Heteroaromatic","Unknown","Functionalized Benzene",
  "Hydrocarbon","Hydrocarbon","Hydrocarbon","Unknown","Heteroaromatic","Aldehyde","Heteroaromatic","Functionalized Benzene","Functionalized Benzene","Unknown","Heteroaromatic","Heteroaromatic","Unknown","Heteroaromatic","Ketone","Other","Unknown"
)
length(Shared_FG)

color1 = "chartreuse"
color2 = "darkgreen"
piecolorfnc = colorRampPalette(c(color1,color2))

pie(table(Shared_FG)/length(Shared_FG),col=piecolorfnc(length(table(Shared_FG))),main="Shared VOCs - functional groups")

