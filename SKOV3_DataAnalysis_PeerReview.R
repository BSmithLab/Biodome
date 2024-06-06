
## Figure 2 in Biodome Publication

#SK-OV-3 Data Analysis
#Written by: Jarrett Eshima
#Summer 2023

# Read the data in
wd = "D:/Biodome Publication/MarchReviews/Peer Review/Reviews/Replication/AlignedData"
setwd(wd)
SKOV = read.csv("Figure2_Clean.csv")

########################################## Pre-Processing and Filtering  #######################################################################

repromedia = reproempty = reprocells = rep(NA,nrow(SKOV))
for(i in 1:nrow(SKOV)){
  
  media = c(SKOV$Biodome_MediaBlank_Day1_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day2_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day3_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day4_7.12.22.1[i],SKOV$Media_Day1_R2_1.25.24.1[i],SKOV$Media_Day4_R2_1.25.24.1[i],SKOV$Media_R3_Day1_2.4.24.1[i],SKOV$Media_R3_Day2_2.4.24.1[i],SKOV$Media_R3_Day3_2.4.24.1[i],SKOV$Media_R3_Day4_2.4.24.1[i]) #Instrument failed to run two media controls
  empty = c(SKOV$EmptyBiodome_Day1_7.18.22.1[i],SKOV$EmptyBiodome_Day2_7.18.22.1[i],SKOV$EmptyBiodome_Day3_7.18.22.1[i],SKOV$EmptyBiodome_Day4_7.18.22.1[i],SKOV$Empty_Day1_R2_1.30.24.1[i],SKOV$Empty_Day2_R2_1.30.24.1[i],SKOV$Empty_Day3_R2_1.30.24.1[i],SKOV$Empty_Day4_R2_1.30.24.1[i],SKOV$Empty_Day1_R3_1.30.24.1[i],SKOV$Empty_Day2_R3_1.30.24.1[i],SKOV$Empty_Day3_R3_1.30.24.1[i],SKOV$Empty_Day4_R3_1.30.24.1[i])
  cells = c(SKOV$Light_Day1_3.3.23.1[i],SKOV$Light_Day2_3.3.23.1[i],SKOV$Light_Day3_3.3.23.1[i],SKOV$Light_Day4_3.3.23.1[i],SKOV$Light_SKOV3_Day1_R2.1[i],SKOV$Light_SKOV3_Day2_R2.1[i],SKOV$Light_SKOV3_Day3_R2.1[i],SKOV$Light_SKOV3_Day4_R2.1[i],SKOV$Light_SKOV3_Day1_R3.1[i],SKOV$Light_SKOV3_Day2_R3.1[i],SKOV$Light_SKOV3_Day3_R3.1[i],SKOV$Light_SKOV3_Day4_R3.1[i])
  
  mediacount = sum(is.na(media))
  emptycount = sum(is.na(empty))
  cellcount = sum(is.na(cells))
  
  if(mediacount <= 5){
    repromedia[i] = i
  }
  
  if(emptycount <= 6){
    reproempty[i] = i
  }
  
  if(cellcount <= 6){
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

#Remove features found in all blanks
inblank = rep(NA,nrow(ReproducibleVOCs))
for(i in 1:nrow(ReproducibleVOCs)){
  if(sum(is.na(c(ReproducibleVOCs$Blank2_1.30.24.1[i],ReproducibleVOCs$Blank_1.8.24.1[i],ReproducibleVOCs$Blank_7.18.22.2[i],ReproducibleVOCs$Blank_2.4.24.1[i],ReproducibleVOCs$BLANK_3.3.23.1[i]))) == 0){
    inblank[i] = i
  }
}

inblank = inblank[!is.na(inblank)]

ReproducibleVOCs2 = ReproducibleVOCs[-inblank,] 

#Remove less than 400s
ReproducibleVOCs2$AvgDim1 = NA
ReproducibleVOCs2$AvgDim2 = NA

ReproducibleVOCs2$AvgDim1 = rowMeans(ReproducibleVOCs2[seq(4,ncol(ReproducibleVOCs2),4)],na.rm = T)
ReproducibleVOCs2$AvgDim2 = rowMeans(ReproducibleVOCs2[seq(5,ncol(ReproducibleVOCs2),4)],na.rm = T)

ReproducibleVOCs3 = ReproducibleVOCs2[order(ReproducibleVOCs2$AvgDim1),]

poorres = rep(NA,nrow(ReproducibleVOCs3))
for(i in 1:nrow(ReproducibleVOCs3)){
  if(ReproducibleVOCs3$AvgDim1[i] < 400){
    poorres[i] = i
  }
}
poorres = poorres[!is.na(poorres)]

ReproducibleVOCs3 = ReproducibleVOCs3[-poorres,]

setwd("D:/Biodome Publication/MarchReviews/Peer Review/Reviews/Replication/Light SKOV")
#write.csv(ReproducibleVOCs2,"BiodomeVOCs_Figure2_50percentthresh_PeerReview_5-21-24.csv")


#######################################

SKOV = ReproducibleVOCs3

#Empty only
keep3 = rep(NA,nrow(SKOV))
for(i in 1:nrow(SKOV)){
  
  media = c(SKOV$Biodome_MediaBlank_Day1_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day2_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day3_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day4_7.12.22.1[i],SKOV$Media_Day1_R2_1.25.24.1[i],SKOV$Media_Day4_R2_1.25.24.1[i],SKOV$Media_R3_Day1_2.4.24.1[i],SKOV$Media_R3_Day2_2.4.24.1[i],SKOV$Media_R3_Day3_2.4.24.1[i],SKOV$Media_R3_Day4_2.4.24.1[i])
  empty = c(SKOV$EmptyBiodome_Day1_7.18.22.1[i],SKOV$EmptyBiodome_Day2_7.18.22.1[i],SKOV$EmptyBiodome_Day3_7.18.22.1[i],SKOV$EmptyBiodome_Day4_7.18.22.1[i],SKOV$Empty_Day1_R2_1.30.24.1[i],SKOV$Empty_Day2_R2_1.30.24.1[i],SKOV$Empty_Day3_R2_1.30.24.1[i],SKOV$Empty_Day4_R2_1.30.24.1[i],SKOV$Empty_Day1_R3_1.30.24.1[i],SKOV$Empty_Day2_R3_1.30.24.1[i],SKOV$Empty_Day3_R3_1.30.24.1[i],SKOV$Empty_Day4_R3_1.30.24.1[i])
  cells = c(SKOV$Light_Day1_3.3.23.1[i],SKOV$Light_Day2_3.3.23.1[i],SKOV$Light_Day3_3.3.23.1[i],SKOV$Light_Day4_3.3.23.1[i],SKOV$Light_SKOV3_Day1_R2.1[i],SKOV$Light_SKOV3_Day2_R2.1[i],SKOV$Light_SKOV3_Day3_R2.1[i],SKOV$Light_SKOV3_Day4_R2.1[i],SKOV$Light_SKOV3_Day1_R3.1[i],SKOV$Light_SKOV3_Day2_R3.1[i],SKOV$Light_SKOV3_Day3_R3.1[i],SKOV$Light_SKOV3_Day4_R3.1[i])
  
  mediacount = sum(is.na(media))
  emptycount = sum(is.na(empty))
  cellcount = sum(is.na(cells))
  
  if(emptycount <= 6){
    if(cellcount >= 7){
        keep3[i] = i
    }
  }
}
keep3 = keep3[!is.na(keep3)]

EmptyOnly = SKOV[keep3,]


SD1 = EmptyOnly
SD1[,seq(6,ncol(GlobalVOCs),4)] = log10(SD1[,seq(6,ncol(GlobalVOCs),4)])

#write.csv(SD1,"Supplemental_Dataset_1_PeerReview_5-21-24.csv")

#######################################

#Media Origin

keep = rep(NA,nrow(SKOV))
for(i in 1:nrow(SKOV)){
  
  media = c(SKOV$Biodome_MediaBlank_Day1_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day2_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day3_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day4_7.12.22.1[i],SKOV$Media_Day1_R2_1.25.24.1[i],SKOV$Media_Day4_R2_1.25.24.1[i],SKOV$Media_R3_Day1_2.4.24.1[i],SKOV$Media_R3_Day2_2.4.24.1[i],SKOV$Media_R3_Day3_2.4.24.1[i],SKOV$Media_R3_Day4_2.4.24.1[i])
  empty = c(SKOV$EmptyBiodome_Day1_7.18.22.1[i],SKOV$EmptyBiodome_Day2_7.18.22.1[i],SKOV$EmptyBiodome_Day3_7.18.22.1[i],SKOV$EmptyBiodome_Day4_7.18.22.1[i],SKOV$Empty_Day1_R2_1.30.24.1[i],SKOV$Empty_Day2_R2_1.30.24.1[i],SKOV$Empty_Day3_R2_1.30.24.1[i],SKOV$Empty_Day4_R2_1.30.24.1[i],SKOV$Empty_Day1_R3_1.30.24.1[i],SKOV$Empty_Day2_R3_1.30.24.1[i],SKOV$Empty_Day3_R3_1.30.24.1[i],SKOV$Empty_Day4_R3_1.30.24.1[i])
  cells = c(SKOV$Light_Day1_3.3.23.1[i],SKOV$Light_Day2_3.3.23.1[i],SKOV$Light_Day3_3.3.23.1[i],SKOV$Light_Day4_3.3.23.1[i],SKOV$Light_SKOV3_Day1_R2.1[i],SKOV$Light_SKOV3_Day2_R2.1[i],SKOV$Light_SKOV3_Day3_R2.1[i],SKOV$Light_SKOV3_Day4_R2.1[i],SKOV$Light_SKOV3_Day1_R3.1[i],SKOV$Light_SKOV3_Day2_R3.1[i],SKOV$Light_SKOV3_Day3_R3.1[i],SKOV$Light_SKOV3_Day4_R3.1[i])
  
  mediacount = sum(is.na(media))
  emptycount = sum(is.na(empty))
  cellcount = sum(is.na(cells))
  
  if(mediacount <= 5){
        keep[i] = i
  }
  
}
keep = keep[!is.na(keep)]

MediaOrigin = SKOV[keep,]

MediaOrigin = MediaOrigin[! rownames(MediaOrigin) %in% rownames(EmptyOnly),]


#Media only (possibly indicates consumption)
keep2 = rep(NA,nrow(SKOV))
for(i in 1:nrow(SKOV)){
  
  media = c(SKOV$Biodome_MediaBlank_Day1_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day2_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day3_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day4_7.12.22.1[i],SKOV$Media_Day1_R2_1.25.24.1[i],SKOV$Media_Day4_R2_1.25.24.1[i],SKOV$Media_R3_Day1_2.4.24.1[i],SKOV$Media_R3_Day2_2.4.24.1[i],SKOV$Media_R3_Day3_2.4.24.1[i],SKOV$Media_R3_Day4_2.4.24.1[i])
  empty = c(SKOV$EmptyBiodome_Day1_7.18.22.1[i],SKOV$EmptyBiodome_Day2_7.18.22.1[i],SKOV$EmptyBiodome_Day3_7.18.22.1[i],SKOV$EmptyBiodome_Day4_7.18.22.1[i],SKOV$Empty_Day1_R2_1.30.24.1[i],SKOV$Empty_Day2_R2_1.30.24.1[i],SKOV$Empty_Day3_R2_1.30.24.1[i],SKOV$Empty_Day4_R2_1.30.24.1[i],SKOV$Empty_Day1_R3_1.30.24.1[i],SKOV$Empty_Day2_R3_1.30.24.1[i],SKOV$Empty_Day3_R3_1.30.24.1[i],SKOV$Empty_Day4_R3_1.30.24.1[i])
  cells = c(SKOV$Light_Day1_3.3.23.1[i],SKOV$Light_Day2_3.3.23.1[i],SKOV$Light_Day3_3.3.23.1[i],SKOV$Light_Day4_3.3.23.1[i],SKOV$Light_SKOV3_Day1_R2.1[i],SKOV$Light_SKOV3_Day2_R2.1[i],SKOV$Light_SKOV3_Day3_R2.1[i],SKOV$Light_SKOV3_Day4_R2.1[i],SKOV$Light_SKOV3_Day1_R3.1[i],SKOV$Light_SKOV3_Day2_R3.1[i],SKOV$Light_SKOV3_Day3_R3.1[i],SKOV$Light_SKOV3_Day4_R3.1[i])
  
  mediacount = sum(is.na(media))
  emptycount = sum(is.na(empty))
  cellcount = sum(is.na(cells))
  
  if(mediacount <= 5){
    if(cellcount >= 11){
      if(emptycount >= 11){
        keep2[i] = i
      }
    }
  }
  
}
keep2 = keep2[!is.na(keep2)]

MediaOnly = SKOV[keep2,]


#table(colnames(Clean_MediaOnly3) == colnames(Clean_MediaOrigin2))
MediaOrigin = MediaOrigin[!rownames(MediaOrigin) %in% rownames(MediaOnly),]
MediaOrigin = MediaOrigin[!rownames(MediaOrigin) %in% rownames(EmptyOnly),]
SD2 = MediaOrigin

SD2[,seq(6,ncol(SD2),4)] = log10(SD2[,seq(6,ncol(SD2),4)])

#write.csv(SD2,"Supplemental_Dataset_2_PeerReview_5-21-24.csv")

SD3 = MediaOnly

SD3[,seq(6,ncol(SD3),4)] = log10(SD3[,seq(6,ncol(SD3),4)])

#write.csv(SD3,"Supplemental_Dataset_3_PeerReview_5-21-24.csv")

######################################################################

#SK-OV-3 Only VOCs

#keep4 = FC_Cells_Media = FC_Cells_Empty = rep(NA,nrow(SKOV))
keep4 = rep(NA,nrow(SKOV))
for(i in 1:nrow(SKOV)){
  
  media = c(SKOV$Biodome_MediaBlank_Day1_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day2_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day3_7.12.22.1[i],SKOV$Biodome_MediaBlank_Day4_7.12.22.1[i],SKOV$Media_Day1_R2_1.25.24.1[i],SKOV$Media_Day4_R2_1.25.24.1[i],SKOV$Media_R3_Day1_2.4.24.1[i],SKOV$Media_R3_Day2_2.4.24.1[i],SKOV$Media_R3_Day3_2.4.24.1[i],SKOV$Media_R3_Day4_2.4.24.1[i])
  empty = c(SKOV$EmptyBiodome_Day1_7.18.22.1[i],SKOV$EmptyBiodome_Day2_7.18.22.1[i],SKOV$EmptyBiodome_Day3_7.18.22.1[i],SKOV$EmptyBiodome_Day4_7.18.22.1[i],SKOV$Empty_Day1_R2_1.30.24.1[i],SKOV$Empty_Day2_R2_1.30.24.1[i],SKOV$Empty_Day3_R2_1.30.24.1[i],SKOV$Empty_Day4_R2_1.30.24.1[i],SKOV$Empty_Day1_R3_1.30.24.1[i],SKOV$Empty_Day2_R3_1.30.24.1[i],SKOV$Empty_Day3_R3_1.30.24.1[i],SKOV$Empty_Day4_R3_1.30.24.1[i])
  cells = c(SKOV$Light_Day1_3.3.23.1[i],SKOV$Light_Day2_3.3.23.1[i],SKOV$Light_Day3_3.3.23.1[i],SKOV$Light_Day4_3.3.23.1[i],SKOV$Light_SKOV3_Day1_R2.1[i],SKOV$Light_SKOV3_Day2_R2.1[i],SKOV$Light_SKOV3_Day3_R2.1[i],SKOV$Light_SKOV3_Day4_R2.1[i],SKOV$Light_SKOV3_Day1_R3.1[i],SKOV$Light_SKOV3_Day2_R3.1[i],SKOV$Light_SKOV3_Day3_R3.1[i],SKOV$Light_SKOV3_Day4_R3.1[i])
  
  mediacount = sum(is.na(media))
  emptycount = sum(is.na(empty))
  cellcount = sum(is.na(cells))
  
  if(cellcount <= 6){
        keep4[i] = i
    }
  
  
}
keep4 = keep4[!is.na(keep4)]

CellsOnly = SKOV[keep4,]


CellsOnly = CellsOnly[!rownames(CellsOnly) %in% rownames(MediaOrigin),]
CellsOnly = CellsOnly[!rownames(CellsOnly) %in% rownames(MediaOnly),]
CellsOnly = CellsOnly[!rownames(CellsOnly) %in% rownames(EmptyOnly),]

ST1 = CellsOnly

#ST1 = ST1[! rownames(ST1) %in% rownames(GlobalVOCs),]
#ST1 = ST1[! rownames(ST1) %in% rownames(Clean_MediaOrigin2),]

ST1[,seq(6,ncol(ST1),4)] = log10(ST1[,seq(6,ncol(ST1),4)])

#write.csv(ST1,"Table1_PeerReview_5-21-24_v2.csv")

##########################################  STATS  #######################################################################

PlotSKOV = read.csv("SKOV_VOCs_Abundances_5-21-24.csv")
rownames(PlotSKOV) = PlotSKOV$Peak
PlotSKOV = PlotSKOV[,-1:-2]

for(i in 1:nrow(PlotSKOV)){
  for(j in 1:ncol(PlotSKOV)){
    if(is.na(PlotSKOV[i,j])){
      PlotSKOV[i,j] = 0
    }
  }
}

FinalResults = data.frame(matrix(NA,nrow(PlotSKOV),ncol=5))
colnames(FinalResults) = c("VOC","pval_SvM","pval_SvE","adjp_SvM","adjp_SvE")
FinalResults$VOC = rownames(PlotSKOV)

Plot_Media = PlotSKOV[,18:27]
Plot_Empty = PlotSKOV[,28:39]
Plot_SKOV = PlotSKOV[,6:17]

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
FinalResults$adjp_SvM = adjpval[1:nrow(FinalResults)]
FinalResults$adjp_SvE = adjpval[(1+nrow(FinalResults)):length(adjpval)]

#write.csv(FinalResults,"Table1_pvalues_PeerReview_5-23-24.csv")

##########################################  Publication plots  #######################################################################

library(beeswarm)
for(i in 1:nrow(PlotSKOV)){
  
  a = as.numeric(Plot_SKOV[i,])
  b = as.numeric(Plot_Media[i,])
  c = as.numeric(Plot_Empty[i,])
  
  boxplot(a,b,c,main=paste(rownames(PlotSKOV)[i]),col=c("#d1641b","#8ebade","#845d9c"),xaxt = "n",ylab="Integrated Unique Mass, Log10 Normalized",cex.main=2,cex.lab=1.5,cex.axis=1.5, ylim = c(0,8))
  #Variable y
  #boxplot(a,b,c,main=paste(rownames(PlotSKOV)[i]),col=c("#d1641b","#8ebade","#845d9c"),xaxt = "n",ylab="Integrated Unique Mass, Log10 Normalized",cex.main=2,cex.lab=1.5,cex.axis=1.5, ylim = c(min(c(a,b,c))-0.1*min(c(a,b,c)),max(c(a,b,c)+0.2*max(c(a,b,c)))))
  
  if(max(a) > max(c)){
    legend(2.7,8,legend = c("SK-OV-3","RPMI1640 Media","System"),col=c("#d1641b","#8ebade","#845d9c"),pch=15,cex = 1.25)
  }else{
    legend(0.5,8,legend = c("SK-OV-3","RPMI1640 Media","System"),col=c("#d1641b","#8ebade","#845d9c"),pch=15,cex=1.25)
  }
  
  #segments(1,max(c(a,b))+0.2*max(c(a,b)),2,max(c(a,b))+0.2*max(c(a,b)),lwd=3)
  #text(1.5,max(c(a,b))+0.25*max(c(a,b)),labels = paste("Adjusted p-value: ",round(adjpval[i],3),sep=""))
  
  beeswarm(a,method = "center",col="#961403",pch=19,cex=1.2,add=T,at = 1)
  beeswarm(b,method = "center",col="#053166",pch=19,cex=1.2,add=T,at = 2)
  beeswarm(c,method = "center",col="#280247",pch=19,cex=1.2,add=T,at = 3)
}


days = seq(1,4)
nrep=3
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
  
  
  Day1_media_mean = mean(c(Plot_Media$Biodome_MediaBlank_Day1_7.12.22.1[i],Plot_Media$Media_Day1_R2_1.25.24.1[i],Plot_Media$Media_R3_Day1_2.4.24.1[i]))
  Day2_media_mean = mean(c(Plot_Media$Biodome_MediaBlank_Day2_7.12.22.1[i],Plot_Media$Media_R3_Day2_2.4.24.1[i]))
  Day3_media_mean = mean(c(Plot_Media$Biodome_MediaBlank_Day3_7.12.22.1[i],Plot_Media$Media_R3_Day3_2.4.24.1[i]))
  Day4_media_mean = mean(c(Plot_Media$Biodome_MediaBlank_Day4_7.12.22.1[i],Plot_Media$Media_Day4_R2_1.25.24.1[i],Plot_Media$Media_R3_Day4_2.4.24.1[i]))
  
  Day1_media_se = sd(c(Plot_Media$Biodome_MediaBlank_Day1_7.12.22.1[i],Plot_Media$Media_Day1_R2_1.25.24.1[i],Plot_Media$Media_R3_Day1_2.4.24.1[i]))/sqrt(nrep)
  Day2_media_se = sd(c(Plot_Media$Biodome_MediaBlank_Day2_7.12.22.1[i],Plot_Media$Media_R3_Day2_2.4.24.1[i]))/sqrt(2)
  Day3_media_se = sd(c(Plot_Media$Biodome_MediaBlank_Day3_7.12.22.1[i],Plot_Media$Media_R3_Day3_2.4.24.1[i]))/sqrt(2)
  Day4_media_se = sd(c(Plot_Media$Biodome_MediaBlank_Day4_7.12.22.1[i],Plot_Media$Media_Day4_R2_1.25.24.1[i],Plot_Media$Media_R3_Day4_2.4.24.1[i]))/sqrt(nrep)
  
  Media_Avg = c(Day1_media_mean,Day2_media_mean,Day3_media_mean,Day4_media_mean)
  Media_SE = c(Day1_media_se,Day2_media_se,Day3_media_se,Day4_media_se)
  
  
  Day1_empty_mean = mean(c(Plot_Empty$EmptyBiodome_Day1_7.18.22.1[i],Plot_Empty$Empty_Day1_R2_1.30.24.1[i],Plot_Empty$Empty_Day1_R3_1.30.24.1[i]))
  Day2_empty_mean = mean(c(Plot_Empty$EmptyBiodome_Day2_7.18.22.1[i],Plot_Empty$Empty_Day2_R2_1.30.24.1[i],Plot_Empty$Empty_Day2_R3_1.30.24.1[i]))
  Day3_empty_mean = mean(c(Plot_Empty$EmptyBiodome_Day3_7.18.22.1[i],Plot_Empty$Empty_Day3_R2_1.30.24.1[i],Plot_Empty$Empty_Day3_R3_1.30.24.1[i]))
  Day4_empty_mean = mean(c(Plot_Empty$EmptyBiodome_Day4_7.18.22.1[i],Plot_Empty$Empty_Day4_R2_1.30.24.1[i],Plot_Empty$Empty_Day4_R3_1.30.24.1[i]))
  
  Day1_empty_se = sd(c(Plot_Empty$EmptyBiodome_Day1_7.18.22.1[i],Plot_Empty$Empty_Day1_R2_1.30.24.1[i],Plot_Empty$Empty_Day1_R3_1.30.24.1[i]))/sqrt(nrep)
  Day2_empty_se = sd(c(Plot_Empty$EmptyBiodome_Day2_7.18.22.1[i],Plot_Empty$Empty_Day2_R2_1.30.24.1[i],Plot_Empty$Empty_Day2_R3_1.30.24.1[i]))/sqrt(nrep)
  Day3_empty_se = sd(c(Plot_Empty$EmptyBiodome_Day3_7.18.22.1[i],Plot_Empty$Empty_Day3_R2_1.30.24.1[i],Plot_Empty$Empty_Day3_R3_1.30.24.1[i]))/sqrt(nrep)
  Day4_empty_se = sd(c(Plot_Empty$EmptyBiodome_Day4_7.18.22.1[i],Plot_Empty$Empty_Day4_R2_1.30.24.1[i],Plot_Empty$Empty_Day4_R3_1.30.24.1[i]))/sqrt(nrep)
  
  Empty_Avg = c(Day1_empty_mean,Day2_empty_mean,Day3_empty_mean,Day4_empty_mean)
  Empty_SE = c(Day1_empty_se,Day2_empty_se,Day3_empty_se,Day4_empty_se)
  
  
  
  plot(days,SKOV_Avg,col="#961403",pch=19,ylim = c(0,8),xlab = "Time (days)",ylab="Log10 Normalized Abundance (Unique Mass)",main=paste(rownames(PlotSKOV)[i]),cex.main=2,cex.lab=1.2,cex.axis=1.3,cex=1.2)
  #Variable y
  #plot(days,a,col="#961403",pch=19,ylim = c(min(c(a,b))-0.1*min(c(a,b)),max(c(a,b)+0.1*max(c(a,b)))),xlab = "Time (days)",ylab="Log10 Normalized Abundance (Unique Mass)",main=paste(rownames(PlotSKOV)[i]),cex.main=2,cex.lab=1.2,cex.axis=1.3,cex=1.2)
  lines(days,SKOV_Avg,col="#d1641b",lwd=2)
  lines(days,Media_Avg,col="#8ebade",lwd=2)
  points(days,Media_Avg,col="#053166",pch=19)
  lines(days,Empty_Avg,col="#845d9c",lwd=2)
  points(days,Empty_Avg,col="#280247",pch=19)
  
  legend(3.3,8,legend = c("SK-OV-3","Media","Flow System"),col = c("#d1641b","#8ebade","#845d9c"),lty=1,lwd=4,cex=1.2)
  
  
  arrows(days, SKOV_Avg+SKOV_SE, days, SKOV_Avg-SKOV_SE, angle=90, code=3, length=0.08, col="#d1641b")
  arrows(days, Media_Avg+Media_SE, days, Media_Avg-Media_SE, angle=90, code=3, length=0.08, col="#8ebade")
  arrows(days, Empty_Avg+Empty_SE, days, Empty_Avg-Empty_SE, angle=90, code=3, length=0.08, col="#845d9c")
  
}
  

####### PCA - All VOCs

QualityVolatilome = SKOV[,seq(6,ncol(SKOV),4)] = log10(SKOV[,seq(6,ncol(SKOV),4)])
rownames(QualityVolatilome) = make.unique(SKOV$Peak)

blanks = c(5,10,15,26,35)
QualityVolatilome = QualityVolatilome[,-blanks]

for(i in 1:nrow(QualityVolatilome)){
  for(j in 1:ncol(QualityVolatilome)){
    if(is.na(QualityVolatilome[i,j])){
      QualityVolatilome[i,j] = 0
    }
  }
}

for(i in 1:nrow(QualityVolatilome)){
  for(j in 1:ncol(QualityVolatilome)){
    if(is.infinite(QualityVolatilome[i,j])){
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

mycolors = c("#8ebade","#8ebade","#8ebade","#8ebade","#845d9c","#845d9c","#845d9c","#845d9c",rep("#d1641b",12),"#8ebade","#8ebade",rep("#845d9c",8),"#8ebade","#8ebade","#8ebade","#8ebade")

plot(PC1,PC2,xlim=c(-40,40),ylim=c(-40,40),main="SK-OV-3 and Control Volatilomes",col=mycolors,pch=20,cex=3,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5)
legend(20,40,legend = c("SK-OV-3","RPMI1640 Media","System"),col = c("#d1641b","#8ebade","#845d9c"),pch=20,pt.cex = 1.5)


####### PCA - SKOV3 VOCs

QualityVolatilome2 = PlotSKOV
blanks = c(1,2,3,4,5,40,41)
QualityVolatilome2 = QualityVolatilome2[,-blanks]

pca = prcomp(t(QualityVolatilome2))
PC1 = c(pca$x[,1])
PC2 = c(pca$x[,2])

#Get axis bounds
min1 = min(PC1)
max1 = max(PC1)
min2 = min(PC2)
max2 = max(PC2)

mycolors = c(rep("#d1641b",12),rep("#8ebade",10),rep("#845d9c",12))
  
plot(PC1,PC2,xlim=c(-10,20),ylim=c(-10,10),main="SK-OV-3 Volatilomes",col=mycolors,pch=20,cex=3,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5)
legend(16,30,legend = c("SK-OV-3","RPMI1640 Media","System"),col = c("#d1641b","#8ebade","#845d9c"),pch=20,pt.cex = 1.5)


###### Barplot
par(mar=c(5,5,5,5))
VOCnumber = c(31,121,102)
bpcolors = c("#d1641b","#8ebade","#845d9c")
barplot(VOCnumber,col = bpcolors,ylim=c(0,250),main="Endogenous and Exogenous VOCs in the Biodome In Vitro System",ylab = "Number of VOCs",cex.axis = 1.5,cex.main=1.5,cex.lab=1.5)
abline(h=0,col="black")
axis(1,at=c(0.68,1.9,3.1),labels = c("Cellular Origin (Putative)","Media Origin","System Origin"),cex=1.5,cex.axis=1.25)
