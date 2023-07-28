
## Figure 2 in Biodome Publication

#SK-OV-3 Data Analysis
#Written by: Jarrett Eshima
#Summer 2023

# Read the data in
wd = "D:/Biodome Publication/Biodome Manuscript/Data/GCxGC/Endogenous"
setwd(wd)

SKOV = read.csv("Figure2_Main_ready.csv")

########################################## Pre-Processing and Filtering  #######################################################################

repromedia = reproempty = reprocells = rep(NA,nrow(SKOV))
for(i in 1:nrow(SKOV)){
  
  media = c(SKOV$Biodome_MediaBlank_Day1_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day2_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day3_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day4_7.12.22.1[i])
  empty = c(SKOV$EmptyBiodome_Day1_7.18.22.1[i],SKOV$EmptyBiodome_Day2_7.18.22.1[i],SKOV$EmptyBiodome_Day3_7.18.22.1[i],SKOV$EmptyBiodome_Day4_7.18.22.1[i])
  cells = c(SKOV$Light_Day1_3.3.23.1[i],SKOV$Light_Day2_3.3.23.1[i],SKOV$Light_Day3_3.3.23.1[i],SKOV$Light_Day4_3.3.23.1[i])
  
  mediacount = sum(is.na(media))
  emptycount = sum(is.na(empty))
  cellcount = sum(is.na(cells))
  
  if(mediacount == 0 || mediacount == 1){
    repromedia[i] = i
  }
  
  if(emptycount == 0 || emptycount == 1){
    reproempty[i] = i
  }
  
  if(cellcount == 0 || cellcount == 1){
    reprocells[i] = i
  }
  
}

repromedia = repromedia[!is.na(repromedia)]
reproempty = reproempty[!is.na(reproempty)]
reprocells = reprocells[!is.na(reprocells)]

MediaVOCs = SKOV[repromedia,]
EmptyVOCs = SKOV[reproempty,]
CellsVOCs = SKOV[reprocells,]

######################################################

#Remove contaminants

#Cells
CellContam1 = which(CellsVOCs$Mass == 44)
CellContam2 = which(CellsVOCs$Mass == 73)
CellContam3 = which(CellsVOCs$Mass == 207)
CellContam4 = which(CellsVOCs$Mass == 281)
CellContam = as.numeric(names(table(c(CellContam1,CellContam2,CellContam3,CellContam4))))
Filt_CellsVOCs = CellsVOCs[-CellContam,]

#Media
MediaContam1 = which(MediaVOCs$Mass == 44)
MediaContam2 = which(MediaVOCs$Mass == 73)
MediaContam3 = which(MediaVOCs$Mass == 207)
MediaContam4 = which(MediaVOCs$Mass == 281)
MediaContam = as.numeric(names(table(c(MediaContam1,MediaContam2,MediaContam3,MediaContam4))))
Filt_MediaVOCs = MediaVOCs[-MediaContam,]

#Empty
EmptyContam1 = which(EmptyVOCs$Mass == 44)
EmptyContam2 = which(EmptyVOCs$Mass == 73)
EmptyContam3 = which(EmptyVOCs$Mass == 207)
EmptyContam4 = which(EmptyVOCs$Mass == 281)
EmptyContam = as.numeric(names(table(c(EmptyContam1,EmptyContam2,EmptyContam3,EmptyContam4))))
Filt_EmptyVOCs = EmptyVOCs[-EmptyContam,]

ReproIndex = as.numeric(names(table(c(rownames(Filt_CellsVOCs),rownames(Filt_MediaVOCs),rownames(Filt_EmptyVOCs)))))
ReproIndex = ReproIndex[order(ReproIndex,decreasing = F)]
ReproducibleVOCs = SKOV[ReproIndex,]
NonRepro = SKOV[-ReproIndex,]

#write.csv(ReproducibleVOCs,"BiodomeVOCs_Figure2_75percentthresh.csv")


#######################################

#Global VOCs

tmp = rownames(Filt_CellsVOCs)[rownames(Filt_CellsVOCs) %in% rownames(Filt_EmptyVOCs)]
tmp = as.numeric(tmp)
tmp2 = tmp[tmp %in% rownames(Filt_MediaVOCs)]

GlobalVOCs = Filt_MediaVOCs[rownames(Filt_MediaVOCs) %in% tmp2,]

#write.csv(GlobalVOCs,"Biodome_GlobalVOCs_Figure2.csv")


SD1 = GlobalVOCs

SD1$Biodome_MediaBlank_Day1_7.12.22.1 = log10(SD1$Biodome_MediaBlank_Day1_7.12.22.1)
SD1$Biodome_MediaBlank_Day2_7.12.22.1 = log10(SD1$Biodome_MediaBlank_Day2_7.12.22.1)
SD1$Biodome_MediaBlank_Day3_7.12.22.1 = log10(SD1$Biodome_MediaBlank_Day3_7.12.22.1)
SD1$Biodome_MediaBlank_Day4_7.12.22.1 = log10(SD1$Biodome_MediaBlank_Day4_7.12.22.1)

SD1$EmptyBiodome_Day1_7.18.22.1 = log10(SD1$EmptyBiodome_Day1_7.18.22.1)
SD1$EmptyBiodome_Day2_7.18.22.1 = log10(SD1$EmptyBiodome_Day2_7.18.22.1)
SD1$EmptyBiodome_Day3_7.18.22.1 = log10(SD1$EmptyBiodome_Day3_7.18.22.1)
SD1$EmptyBiodome_Day4_7.18.22.1 = log10(SD1$EmptyBiodome_Day4_7.18.22.1)

SD1$Light_Day1_3.3.23.1 = log10(SD1$Light_Day1_3.3.23.1)
SD1$Light_Day2_3.3.23.1 = log10(SD1$Light_Day2_3.3.23.1)
SD1$Light_Day3_3.3.23.1 = log10(SD1$Light_Day3_3.3.23.1)
SD1$Light_Day4_3.3.23.1 = log10(SD1$Light_Day4_3.3.23.1)

#write.csv(SD1,"Supplemental_Dataset_1.csv")

#######################################

#Media and cell culture

keep = rep(NA,nrow(SKOV))
for(i in 1:nrow(SKOV)){
  
  media = c(SKOV$Biodome_MediaBlank_Day1_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day2_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day3_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day4_7.12.22.1[i])
  empty = c(SKOV$EmptyBiodome_Day1_7.18.22.1[i],SKOV$EmptyBiodome_Day2_7.18.22.1[i],SKOV$EmptyBiodome_Day3_7.18.22.1[i],SKOV$EmptyBiodome_Day4_7.18.22.1[i])
  cells = c(SKOV$Light_Day1_3.3.23.1[i],SKOV$Light_Day2_3.3.23.1[i],SKOV$Light_Day3_3.3.23.1[i],SKOV$Light_Day4_3.3.23.1[i])
  
  mediacount = sum(is.na(media))
  emptycount = sum(is.na(empty))
  cellcount = sum(is.na(cells))
  
  if(mediacount == 0 || mediacount == 1){
    if(cellcount == 0 || cellcount == 1){
        keep[i] = i
    }
  }
  
}
keep = keep[!is.na(keep)]

MediaOrigin = SKOV[keep,]

Clean_MediaOrigin = MediaOrigin[! rownames(MediaOrigin) %in% rownames(GlobalVOCs),]

MContam1 = which(Clean_MediaOrigin$Mass == 44)
MContam2 = which(Clean_MediaOrigin$Mass == 73)
MContam3 = which(Clean_MediaOrigin$Mass == 207)
MContam4 = which(Clean_MediaOrigin$Mass == 281)
MContam = as.numeric(names(table(c(MContam1,MContam2,MContam3,MContam4))))
Clean_MediaOrigin2 = Clean_MediaOrigin[-MContam,]

#Media only (possibly indicates consumption)
keep2 = rep(NA,nrow(SKOV))
for(i in 1:nrow(SKOV)){
  
  media = c(SKOV$Biodome_MediaBlank_Day1_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day2_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day3_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day4_7.12.22.1[i])
  empty = c(SKOV$EmptyBiodome_Day1_7.18.22.1[i],SKOV$EmptyBiodome_Day2_7.18.22.1[i],SKOV$EmptyBiodome_Day3_7.18.22.1[i],SKOV$EmptyBiodome_Day4_7.18.22.1[i])
  cells = c(SKOV$Light_Day1_3.3.23.1[i],SKOV$Light_Day2_3.3.23.1[i],SKOV$Light_Day3_3.3.23.1[i],SKOV$Light_Day4_3.3.23.1[i])
  
  mediacount = sum(is.na(media))
  emptycount = sum(is.na(empty))
  cellcount = sum(is.na(cells))
  
  if(mediacount == 0 || mediacount == 1){
    if(cellcount == 3 || cellcount == 4){
      if(emptycount == 3 || emptycount == 4){
        keep2[i] = i
      }
    }
  }
  
}
keep2 = keep2[!is.na(keep2)]

MediaOnly = SKOV[keep2,]

Clean_MediaOnly = MediaOnly[! rownames(MediaOnly) %in% rownames(GlobalVOCs),]
Clean_MediaOnly2 = Clean_MediaOnly[! rownames(Clean_MediaOnly) %in% rownames(Clean_MediaOrigin2),]

MOContam1 = which(Clean_MediaOnly2$Mass == 44)
MOContam2 = which(Clean_MediaOnly2$Mass == 73)
MOContam3 = which(Clean_MediaOnly2$Mass == 207)
MOContam4 = which(Clean_MediaOnly2$Mass == 281)
MOContam = as.numeric(names(table(c(MOContam1,MOContam2,MOContam3,MOContam4))))
Clean_MediaOnly3 = Clean_MediaOnly2[-MOContam,]

table(colnames(Clean_MediaOnly3) == colnames(Clean_MediaOrigin2))
SD2 = rbind(Clean_MediaOrigin2,Clean_MediaOnly3) 

#write.csv(SD2,"Supplemental_Dataset_2.csv")

######################################################################

#SK-OV-3 Only VOCs
keep3 = rep(NA,nrow(SKOV))
for(i in 1:nrow(SKOV)){
  
  media = c(SKOV$Biodome_MediaBlank_Day1_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day2_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day3_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day4_7.12.22.1[i])
  empty = c(SKOV$EmptyBiodome_Day1_7.18.22.1[i],SKOV$EmptyBiodome_Day2_7.18.22.1[i],SKOV$EmptyBiodome_Day3_7.18.22.1[i],SKOV$EmptyBiodome_Day4_7.18.22.1[i])
  cells = c(SKOV$Light_Day1_3.3.23.1[i],SKOV$Light_Day2_3.3.23.1[i],SKOV$Light_Day3_3.3.23.1[i],SKOV$Light_Day4_3.3.23.1[i])
  
  mediacount = sum(is.na(media))
  emptycount = sum(is.na(empty))
  cellcount = sum(is.na(cells))
  
  if(mediacount == 4){
    if(cellcount == 0 || cellcount == 1){
      if(emptycount == 4){
        keep3[i] = i
      }
    }
  }
  
}
keep3 = keep3[!is.na(keep3)]

CellsOnly = SKOV[keep3,]

keep4 = rep(NA,nrow(SKOV))
for(i in 1:nrow(SKOV)){
  
  media = c(SKOV$Biodome_MediaBlank_Day1_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day2_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day3_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day4_7.12.22.1[i])
  empty = c(SKOV$EmptyBiodome_Day1_7.18.22.1[i],SKOV$EmptyBiodome_Day2_7.18.22.1[i],SKOV$EmptyBiodome_Day3_7.18.22.1[i],SKOV$EmptyBiodome_Day4_7.18.22.1[i])
  cells = c(SKOV$Light_Day1_3.3.23.1[i],SKOV$Light_Day2_3.3.23.1[i],SKOV$Light_Day3_3.3.23.1[i],SKOV$Light_Day4_3.3.23.1[i])
  
  mediacount = sum(is.na(media))
  emptycount = sum(is.na(empty))
  cellcount = sum(is.na(cells))
  
  if(mediacount == 3 || mediacount == 4){
    if(cellcount == 0 || cellcount == 1){
      if(emptycount == 3 || emptycount == 4){
        keep4[i] = i
      }
    }
  }
  
}
keep4 = keep4[!is.na(keep4)]

CellsOnly2 = SKOV[keep4,]
CellsOnly3 = CellsOnly2[!rownames(CellsOnly2) %in% rownames(CellsOnly),]

CellsOnlyFinal = rbind(CellsOnly,CellsOnly3)

CContam1 = which(CellsOnlyFinal$Mass == 44)
CContam2 = which(CellsOnlyFinal$Mass == 73)
CContam3 = which(CellsOnlyFinal$Mass == 207)
CContam4 = which(CellsOnlyFinal$Mass == 281)
CContam = as.numeric(names(table(c(CContam1,CContam2,CContam3,CContam4))))
ST1 = CellsOnlyFinal[-CContam,]

ST1$Biodome_MediaBlank_Day1_7.12.22.1 = log10(ST1$Biodome_MediaBlank_Day1_7.12.22.1)
ST1$Biodome_MediaBlank_Day2_7.12.22.1 = log10(ST1$Biodome_MediaBlank_Day2_7.12.22.1)
ST1$Biodome_MediaBlank_Day3_7.12.22.1 = log10(ST1$Biodome_MediaBlank_Day3_7.12.22.1)
ST1$Biodome_MediaBlank_Day4_7.12.22.1 = log10(ST1$Biodome_MediaBlank_Day4_7.12.22.1)

ST1$EmptyBiodome_Day1_7.18.22.1 = log10(ST1$EmptyBiodome_Day1_7.18.22.1)
ST1$EmptyBiodome_Day2_7.18.22.1 = log10(ST1$EmptyBiodome_Day2_7.18.22.1)
ST1$EmptyBiodome_Day3_7.18.22.1 = log10(ST1$EmptyBiodome_Day3_7.18.22.1)
ST1$EmptyBiodome_Day4_7.18.22.1 = log10(ST1$EmptyBiodome_Day4_7.18.22.1)

ST1$Light_Day1_3.3.23.1 = log10(ST1$Light_Day1_3.3.23.1)
ST1$Light_Day2_3.3.23.1 = log10(ST1$Light_Day2_3.3.23.1)
ST1$Light_Day3_3.3.23.1 = log10(ST1$Light_Day3_3.3.23.1)
ST1$Light_Day4_3.3.23.1 = log10(ST1$Light_Day4_3.3.23.1)

ST1$BLANK_3.3.23.1 = log10(ST1$BLANK_3.3.23.1)
ST1$Blank_7.18.22.2 = log10(ST1$Blank_7.18.22.2)


#write.csv(ST1,"Table_S1.csv")

##########################################  STATS  #######################################################################

PlotSKOV = read.csv("SKOV_VOCs_Abundance.csv")
rownames(PlotSKOV) = PlotSKOV$Peak
PlotSKOV = PlotSKOV[,-1]

FinalResults = data.frame(matrix(NA,nrow(PlotSKOV),ncol=5))
colnames(FinalResults) = c("VOC","pval_SvM","pval_SvE","adjp_SvM","adjp_SvE")
FinalResults$VOC = rownames(PlotSKOV)

Plot_Media = PlotSKOV[,1:4]
Plot_Empty = PlotSKOV[,5:8]
Plot_SKOV = PlotSKOV[,9:12]

for(i in 1:nrow(PlotSKOV)){
  
  a = as.numeric(Plot_SKOV[i,])
  b = as.numeric(Plot_Media[i,])
  
  res = t.test(a,b,paired = F,alternative = "two.sided")
  
  FinalResults$pval_SvM[i] = res$p.value
  
  c = as.numeric(Plot_SKOV[i,])
  d = as.numeric(Plot_Empty[i,])
  
  res2 = t.test(c,d,paired = F,alternative = "two.sided")
  
  FinalResults$pval_SvE[i] = res2$p.value
  
}

pvals = c(FinalResults$pval_SvM,FinalResults$pval_SvE)
adjpval = p.adjust(pvals,method = "BH")
FinalResults$adjp_SvM = adjpval[1:13]
FinalResults$adjp_SvE = adjpval[14:26]

#write.csv(FinalResults,"Table1_pvalues_7-2-23.csv")

##########################################  Publication plots  #######################################################################

library(beeswarm)
for(i in 1:nrow(PlotSKOV)){
  
  a = as.numeric(Plot_SKOV[i,])
  b = as.numeric(Plot_Media[i,])
  c = as.numeric(Plot_Empty[i,])
  
  boxplot(a,b,c,main=paste(rownames(PlotSKOV)[i]),col=c("#d1641b","#8ebade","#845d9c"),xaxt = "n",ylab="Integrated Unique Mass, Log10 Normalized",cex.main=2,cex.lab=1.5,cex.axis=1.5, ylim = c(0,8))
  axis(1,at=c(1,2,3),labels=c("SK-OV-3","RPMI1640 Media","System Control"))
  #Variable y
  #boxplot(a,b,c,main=paste(rownames(PlotSKOV)[i]),col=c("#d1641b","#8ebade","#845d9c"),xaxt = "n",ylab="Integrated Unique Mass, Log10 Normalized",cex.main=2,cex.lab=1.5,cex.axis=1.5, ylim = c(min(c(a,b,c))-0.1*min(c(a,b,c)),max(c(a,b,c)+0.2*max(c(a,b,c)))))
  
  if(max(a) > max(c)){
    legend(2.7,8,legend = c("SK-OV-3","RPMI1640 Media","Flow System"),col=c("#d1641b","#8ebade","#845d9c"),pch=15,cex = 1.25)
  }else{
    legend(0.5,8,legend = c("SK-OV-3","RPMI1640 Media","Flow System"),col=c("#d1641b","#8ebade","#845d9c"),pch=15,cex=1.25)
  }
  
  #segments(1,max(c(a,b))+0.2*max(c(a,b)),2,max(c(a,b))+0.2*max(c(a,b)),lwd=3)
  #text(1.5,max(c(a,b))+0.25*max(c(a,b)),labels = paste("Adjusted p-value: ",round(adjpval[i],3),sep=""))
  
  beeswarm(a,method = "center",col="#961403",pch=19,cex=1.2,add=T,at = 1)
  beeswarm(b,method = "center",col="#053166",pch=19,cex=1.2,add=T,at = 2)
  beeswarm(c,method = "center",col="#280247",pch=19,cex=1.2,add=T,at = 3)
}


for(i in 1:nrow(PlotSKOV)){
  a = as.numeric(Plot_SKOV[i,])
  b = as.numeric(Plot_Media[i,])
  c = as.numeric(Plot_Empty[i,])
  days = seq(1,4)
  plot(days,a,col="#961403",pch=19,ylim = c(0,8),xlab = "Time (days)",ylab="Log10 Normalized Abundance (Unique Mass)",main=paste(rownames(PlotSKOV)[i]),cex.main=2,cex.lab=1.2,cex.axis=1.3,cex=1.2)
  #Variable y
  #plot(days,a,col="#961403",pch=19,ylim = c(min(c(a,b))-0.1*min(c(a,b)),max(c(a,b)+0.1*max(c(a,b)))),xlab = "Time (days)",ylab="Log10 Normalized Abundance (Unique Mass)",main=paste(rownames(PlotSKOV)[i]),cex.main=2,cex.lab=1.2,cex.axis=1.3,cex=1.2)
  lines(days,a,col="#d1641b",lwd=2)
  lines(days,b,col="#8ebade",lwd=2)
  points(days,b,col="#053166",pch=19)
  lines(days,c,col="#845d9c",lwd=2)
  points(days,c,col="#280247",pch=19)
}

####### PCA

QualityVolatilome = read.csv("Figure2A_PCAdata.csv")
rownames(QualityVolatilome) = make.unique(QualityVolatilome$Peak)
QualityVolatilome = QualityVolatilome[,-1]
QualityVolatilome = log10(QualityVolatilome)

for(i in 1:nrow(QualityVolatilome)){
  for(j in 1:ncol(QualityVolatilome)){
    if(is.na(QualityVolatilome[i,j])){
      QualityVolatilome[i,j] = 0
    }
  }
}

pca = prcomp(t(QualityVolatilome))
PC1 = c(pca$x[,1])
PC2 = c(pca$x[,2])

#Get axis bounds
min1 = min(PC1)
max1 = max(PC1)
min2 = min(PC2)
max2 = max(PC2)

mycolors = c(rep("#8ebade",4),rep("#845d9c",4),rep("#d1641b",4))

plot(PC1,PC2,xlim=c(-40,40),ylim=c(-40,40),main="SK-OV-3 and Control Volatilomes",col=mycolors,pch=20,cex=3.5,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5)
legend(25,40,legend = c("SK-OV-3","RPMI1640 Media","Flow System"),col = c("#d1641b","#8ebade","#845d9c"),pch=20,pt.cex = 1.5)

###### Barplot

VOCnumber = c(13,86,88,197)
bpcolors = c("#d1641b","#8ebade","#845d9c","#bcbbbf")
barplot(VOCnumber,col = bpcolors,ylim=c(0,250),main="Endogenous and Exogenous VOCs in the Biodome In Vitro System",ylab = "Number of VOCs",cex.axis = 1.5,cex.main=1.5,cex.lab=1.5)
abline(h=0,col="black")
axis(1,at=c(0.68,1.9,3.1,4.35),labels = c("Cellular Origin (Putative)","Media Origin","System Origin","Variably Present"),cex=1.5,cex.axis=1.25)

##########################################  Less Consistent VOCs  #######################################################################

Pass = read.csv("BiodomeVOCs_Figure2_75percentthresh_FilteredList_FINAL.csv")
Pass$Biodome_MediaBlank_Day1_7.12.22.1 = log10(Pass$Biodome_MediaBlank_Day1_7.12.22.1)

#Filter Cell Only VOCs
CellFilter = read.csv("Table_S1.csv")

removeindex = rep(NA,nrow(CellFilter))
for(i in 1:nrow(CellFilter)){
  tmpind = which(CellFilter$Peak[i] == Pass$Peak)
  
  if(length(tmpind)>1){
    ind = which(CellFilter$X1st.Dimension.Time..s.[i] == Pass$X1st.Dimension.Time..s.[tmpind])
    
    if(length(ind)>0){
      tmpind = tmpind[ind]
    }else{
      if(is.na(CellFilter$X1st.Dimension.Time..s.[i])){
        ind2 = which(is.na(Pass$X1st.Dimension.Time..s.[tmpind]))
        tmpind = tmpind[ind2]
      }else{
        tmpind = NA
      }
    }
    
    
  }
  
  removeindex[i] = tmpind
}

Pass$Peak[removeindex]

Pass = Pass[-removeindex,]

#Filter Cell Only VOCs
MediaFilter = read.csv("Supplemental_Dataset_2.csv")
colnames(MediaFilter)[1] = "Peak"
MediaFilter = MediaFilter[-1,]

removeindex2 = rep(NA,nrow(MediaFilter))
for(i in 1:nrow(MediaFilter)){
  tmpind = which(MediaFilter$Peak[i] == Pass$Peak)
  
  if(length(tmpind)>1){
    if(MediaFilter$Media.Day.1[i] == ""){
      ind = which(MediaFilter$Media.Day.2[i] == Pass$X1st.Dimension.Time..s..1[tmpind])
    }else{
      ind = which(MediaFilter$Media.Day.1[i] == Pass$X1st.Dimension.Time..s.[tmpind])
    }
    
    if(length(ind)>1){
      indception = which(MediaFilter$X.3[i] == Pass$Biodome_MediaBlank_Day1_7.12.22.1[tmpind])
      ind = ind[indception]
    }
    
    if(length(ind)==1){
      tmpind = tmpind[ind]
    }else{
      if(is.na(MediaFilter$Media.Day.1[i])){
        ind2 = which(is.na(Pass$X1st.Dimension.Time..s.[tmpind]))
        tmpind = tmpind[ind2]
      }else{
        tmpind = NA
      }
    }
    
    
  }
  
  removeindex2[i] = tmpind
}

which(is.na(removeindex2))
removeindex2[50] = 375
removeindex2[88] = 378

Pass$Peak[removeindex2]

Pass = Pass[-removeindex2,]

#Filter Empty Only VOCs
EmptyFilter = read.csv("Supplemental_Dataset_1.csv")
colnames(EmptyFilter)[1] = "Peak"
EmptyFilter = EmptyFilter[-1,]

removeindex3 = rep(NA,nrow(EmptyFilter))
for(i in 1:nrow(EmptyFilter)){
  tmpind = which(EmptyFilter$Peak[i] == Pass$Peak)
  
  if(length(tmpind)>1){
    if(EmptyFilter$Media.Day.1[i] == ""){
      ind = which(EmptyFilter$Media.Day.2[i] == Pass$X1st.Dimension.Time..s..1[tmpind])
    }else{
      ind = which(EmptyFilter$Media.Day.1[i] == Pass$X1st.Dimension.Time..s.[tmpind])
    }
    
    if(length(ind)>1){
      indception = which(EmptyFilter$X.3[i] == Pass$Biodome_MediaBlank_Day1_7.12.22.1[tmpind])
      ind = ind[indception]
    }
    
    if(length(ind)==1){
      tmpind = tmpind[ind]
    }else{
      if(is.na(EmptyFilter$Media.Day.1[i])){
        ind2 = which(is.na(Pass$X1st.Dimension.Time..s.[tmpind]))
        tmpind = tmpind[ind2]
      }else{
        tmpind = NA
      }
    }
    
    
  }
  
  removeindex3[i] = tmpind
}

Pass$Peak[removeindex3]

Pass = Pass[-removeindex3,]

InterestingIndex = c(removeindex,removeindex2,removeindex3)

#write.csv(Pass,"InvestigateVOCs.csv")
weakempty = weakmedia = weakcell = CellEmpty = MediaEmpty = rep(NA,nrow(Pass))
for(i in 1:nrow(Pass)){
  
  media = c(Pass$Biodome_MediaBlank_Day1_7.12.22.1[i],Pass$Biodome_MediaBlank_Day2_7.12.22.1[i],Pass$Biodome_MediaBlank_Day3_7.12.22.1[i],Pass$Biodome_MediaBlank_Day4_7.12.22.1[i])
  empty = c(Pass$EmptyBiodome_Day1_7.18.22.1[i],Pass$EmptyBiodome_Day2_7.18.22.1[i],Pass$EmptyBiodome_Day3_7.18.22.1[i],Pass$EmptyBiodome_Day4_7.18.22.1[i])
  cells = c(Pass$Light_Day1_3.3.23.1[i],Pass$Light_Day2_3.3.23.1[i],Pass$Light_Day3_3.3.23.1[i],Pass$Light_Day4_3.3.23.1[i])
  
  mediacount = sum(is.na(media))
  emptycount = sum(is.na(empty))
  cellcount = sum(is.na(cells))
  
  if(emptycount == 0 || emptycount == 1){
    if(cellcount > 1 && mediacount > 1){
      weakempty[i] = i
    }
  }
  
  if(mediacount == 0 || mediacount == 1){
    if(emptycount > 1 && cellcount > 1)
    weakmedia[i] = i
  }
  
  if(cellcount == 0 || cellcount == 1){
    if(emptycount > 1 && mediacount > 1)
    weakcell[i] = i
  }
  
  if(emptycount == 0 || emptycount == 1){
    if(cellcount == 0 || cellcount == 1){
      CellEmpty[i] = i
    }
  }
  
  if(emptycount == 0 || emptycount == 1){
    if(mediacount == 0 || mediacount == 1){
      MediaEmpty[i] = i
    }
  }
  
}

weakempty = weakempty[!is.na(weakempty)]
weakmedia = weakmedia[!is.na(weakmedia)]
weakcell = weakcell[!is.na(weakcell)]
CellEmpty = CellEmpty[!is.na(CellEmpty)]
MediaEmpty = MediaEmpty[!is.na(MediaEmpty)]

View(Pass[weakcell,])
View(Pass[weakmedia,])
View(Pass[weakempty,])
View(Pass[CellEmpty,])
View(Pass[MediaEmpty,])

#tmpp = c(weakcell,weakmedia,weakempty,CellEmpty,MediaEmpty)
#FinalCheck = Pass[-tmpp,]

#write.csv(Pass,"VariableVOCs_OtherHalf_6-17-23.csv")
