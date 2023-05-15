#######################################  Biodome Automated Analysis  ####################################################################

#Written by: Jarrett Eshima
#Date: February, 2021
#For: Use by the Dr. Barbara Smith Lab at Arizona State University

############## USE

#Optimized to handle incorrectly named VOCs - frequent problem with isotopically labeled metabolites
#Optimized to handle poor alignment results in ChromaTOF (i.e. majority unknowns)
     #Note: ChromaTOF alignment algorithm is superior to this one. However, if ChromaTOF is having difficulties aligning different sample conditions
     #(i.e. Heavy vs Light), then this algorithm can supplement the alignment performed by the instrument to identify the same or similar VOCs based on
     #elution and mass spectral characteristics.

#######################################################################################################################################

#Load Libraries/Dependencies
library(ggplot2) #Plot Results
library(stringr) #for working with strings
library(qvalue) #Multiple Hypothesis Test Adjustment

#Load necessary functions
quickfun = function(x,y){
  return(x[[y]])
}

sumrow = function(x){
  nr = nrow(x)
  tmp = rep(NA,nr)
  for(i in 1:nr){
    tmp[i] = sum(x[i,],na.rm=T)
  }
  
  return(tmp)
}

sumcol = function(x){
  nc = ncol(x)
  tmp = rep(NA,nc)
  for(i in 1:nc){
    tmp[i] = sum(x[,i],na.rm=T)
  }
  
  return(tmp)
}


#######################################################################################################################################
#Pre-processing in ChromaTOF

#Step 1:
#When exporting the aligned compound table from ChromaTOF stat compare, selected the "Count" column, right click and select 'sort'
#Then highlight the compounds found in 4 or more samples (>66%), export selected analytes, delete column superscripts (sample names),

#Step 2:
#replace the column labels with the following:
#First two column headers
#VOC Mass
#Copy and paste the following headers for each sample:
#Class Area Dim1 Dim2

#Step 3:
#add excel column at the end with the following formula (add additional formula parts if more than 6 samples):
#=COUNTIF(C2,"Class1")+COUNTIF(G2,"Class1")+COUNTIF(K2,"Class1")+COUNTIF(O2,"Class1")+COUNTIF(S2,"Class1")+COUNTIF(W2,"Class1")
#Give the column the following header: "Number"

#Step 4: 
#Move all processed files to an empty folder, set that folder as the working directory below.
#All files must be saved as .CSV (comma separated values)

#######################################################################################################################################
#############################  AUTOMATED PEAK ALIGNMENT FOR ISOTOPICALLY LABELED METABOLITES  ##########################################
#######################################################################################################################################

###############################################   User Input   #################################################################

#Copy directory where files are stored/located
#wd = 'F:/Biodome Alignment Algorithm'
#wd = "C:/Users/jeshima/Desktop/Biodome Manuscript/Data/GCxGC/IsotopeLabeling"
wd = "C:/Users/jeshima/Desktop/Biodome Manuscript/Data/GCxGC/IsotopeLabeling/AlignedFiles2"
setwd(wd)

#Samples in each experimental condition (Positive Control, Negative Control, Group 1, Group 2, etc.)
#Depends on frequency of TDU change and VOC collection duration
nsamples = 4

#How many experimental conditions did you run? (blanks/negative controls/positive controls are all separate experimental conditions)
#Note: This algorithm cannot handle more than 5 experimental conditions currently
numberofconditions = 4

#Duplicate Threshold: This values determines the threshold for removing duplicated compound names.
#Note: some duplicate naming is unavoidable due to alignment parameters and the occurrence of "split peaks". A value of 8-10 is recommended to filter "only" the background contaminants
maxdups = 10

#Save the output?
writefile = TRUE
writewd = "C:/Users/jeshima/Desktop/Biodome Manuscript/Data/GCxGC/IsotopeLabeling/AlignedFiles2/Output"

#Give your experimental conditions names (this cannot be left blank, give unique names for each condition)
#IMPORTANT NOTE: EACH CONDITION NAME BELOW MUST BE THE FIRST SET OF CHARACTERS IN THE FILE NAME
#Example: Heavy_Biodome24hrtest_7-29-21 ; Empty_Biodome24hrtest_7-29-21
#NAMING IS CASE SENSITIVE - SO BE EXTRA CAREFUL WHEN NAMING
Condition1Name = "Heavy" #Experimental Condition
Condition2Name = "Light" #Positive Control
Condition3Name = "Media" #Negative Control ("optional")
Condition4Name = "Empty" #Negative Control 2 (optional)
Condition5Name = "" #Optional

## RETENTION TIME TOLERANCE
#User defined first dimension tolerance (in seconds) - typical values are 4-14 (smaller values = stronger restrictions & higher confidence in matching, but more likely to miss some appropriately aligned VOCs due to single outliers)
#Starting value of 10 recommended
Dim1Tol = 10
#User defined second dimension tolerance (in seconds) - typical values are 0.05-0.2
#Starting value of 0.12 recommended
Dim2Tol = 0.07
#Select the reference condition for alignment (for heavy labeling, the heavy condition is recommended as the reference)
ref = "Light"

  
#Note to user: Condition 1 is established as the reference level. Therefore, when interpreting SumRes (SummaryResults) rownames
#are set based on the VOC names assigned to Condition 1. These names may or may not be accurate depending on the quality of naming.
#Columns 2 and 3 show the index of the VOC, from the respective condition, that matched the VOC from Condition 1.

#Note to user: Filtering out compounds present in <50-80% samples is recommended in an effort to focus on high quality VOCs
#Filtering is can be performed using one of the following options:
#1) all files have VOCs removed that are missing in more than X% of samples (not recommended for negative controls, may remove compounds that are important)
#2) only the files corresponding to the experimental condition(s) are filtered
#3) no filtering is performed 
#4) Different filtering thresholds are used depending on the experimental conditions (BE VERY CAREFUL WITH THIS OPTION)

option = 2

############### For Options 1-3

### IF YOU SELECT OPTION 1 OR 2, PLEASE SPECIFY THE FILTERING THRESHOLD ###
#The greater the filtering threshold, the higher 'quality' the VOCs will be (i.e. present in more samples)
#The lower the filtering threshold, the more unknown or poorly resolved VOCs will be included
#A filterthresh value less than 0.5 is not recommended 
FilteringThresh = 0.6 #strictly less than this value

### PARAMETER FOR OPTION 2 ONLY - this can be left blank or not run at all if options 1 or 3 are selected
#Type the names of your experimental conditions
expconditions = c("Light","Heavy")

############### For Option 4 (currently only supports two different filtering thresholds)
FilteringThresh1 = 0.8
FilteringThresh2 = 0.6

expconditions1 = c("Heavy","Light") #This is one 'group' of experimental conditions which should be filtered in the same way
expconditions2 = c("Blank","Empty") #This is the second 'group' 
#Advanced Functionality: You can leave off experimental condition names, which will result in the condition not being filtered.


##Object Naming Scheme
#...Clean... = Background contaminants removed, but no further filtering
#...Final... = Background contaminants removed, and VOCs are filtered according to the option above
#...Output... = Data frame ready for user viewing 


### ONCE YOU HAVE SPECIFIED THE INPUT PARAMETERS ABOVE, HIT "SOURCE" ABOVE TO RUN THE ENTIRE SCRIPT

#AlignInfo files will be in the same order as the ConditionNames specified above

#Check Indices
# X = 44; Y = 2
# 
# cat("Dim 1")
# HighQualityLightDimOne[X,]
# HighQualityHeavyDimOne[Y,]
# 
# cat("Dim 2")
# HighQualityLightDimTwo[X,]
# HighQualityHeavyDimTwo[Y,]

###############################################   READ IN DATA   #################################################################

if(numberofconditions == 1){
  cat("Alignment algorithm is obsolete for a single experimental condition. Use ChromaTOF to perform alignment within an experimental condition.")
}else if(numberofconditions == 2){
  ncharc1 = as.numeric(nchar(Condition1Name))
  ncharc2 = as.numeric(nchar(Condition2Name))
  allfiles = list.files(path=wd,pattern = ".csv")
  C1 = as.numeric(which(substr(allfiles,1,ncharc1) == Condition1Name))
  C2 = as.numeric(which(substr(allfiles,1,ncharc2) == Condition2Name))
  AlignInfo1 = read.csv(allfiles[C1])
  AlignInfo2 = read.csv(allfiles[C2])
  conditions = c(Condition1Name,Condition2Name)
}else if(numberofconditions == 3){
  ncharc1 = as.numeric(nchar(Condition1Name))
  ncharc2 = as.numeric(nchar(Condition2Name))
  ncharc3 = as.numeric(nchar(Condition3Name))
  allfiles = list.files(path=wd,pattern = ".csv")
  C1 = as.numeric(which(substr(allfiles,1,ncharc1) == Condition1Name))
  C2 = as.numeric(which(substr(allfiles,1,ncharc2) == Condition2Name))
  C3 = as.numeric(which(substr(allfiles,1,ncharc3) == Condition3Name))
  AlignInfo1 = read.csv(allfiles[C1])
  AlignInfo2 = read.csv(allfiles[C2])
  AlignInfo3 = read.csv(allfiles[C3])
  conditions = c(Condition1Name,Condition2Name,Condition3Name)
}else if(numberofconditions == 4){
  ncharc1 = as.numeric(nchar(Condition1Name))
  ncharc2 = as.numeric(nchar(Condition2Name))
  ncharc3 = as.numeric(nchar(Condition3Name))
  ncharc4 = as.numeric(nchar(Condition4Name))
  allfiles = list.files(path=wd,pattern = ".csv")
  C1 = as.numeric(which(substr(allfiles,1,ncharc1) == Condition1Name))
  C2 = as.numeric(which(substr(allfiles,1,ncharc2) == Condition2Name))
  C3 = as.numeric(which(substr(allfiles,1,ncharc3) == Condition3Name))
  C4 = as.numeric(which(substr(allfiles,1,ncharc4) == Condition4Name))
  AlignInfo1 = read.csv(allfiles[C1])
  AlignInfo2 = read.csv(allfiles[C2])
  AlignInfo3 = read.csv(allfiles[C3])
  AlignInfo4 = read.csv(allfiles[C4])
  conditions = c(Condition1Name,Condition2Name,Condition3Name,Condition4Name)
}else if(numberofconditions == 5){
  ncharc1 = as.numeric(nchar(Condition1Name))
  ncharc2 = as.numeric(nchar(Condition2Name))
  ncharc3 = as.numeric(nchar(Condition3Name))
  ncharc4 = as.numeric(nchar(Condition4Name))
  ncharc5 = as.numeric(nchar(Condition5Name))
  allfiles = list.files(path=wd,pattern = ".csv")
  C1 = as.numeric(which(substr(allfiles,1,ncharc1) == Condition1Name))
  C2 = as.numeric(which(substr(allfiles,1,ncharc2) == Condition2Name))
  C3 = as.numeric(which(substr(allfiles,1,ncharc3) == Condition3Name))
  C4 = as.numeric(which(substr(allfiles,1,ncharc4) == Condition4Name))
  C5 = as.numeric(which(substr(allfiles,1,ncharc5) == Condition5Name))
  AlignInfo1 = read.csv(allfiles[C1])
  AlignInfo2 = read.csv(allfiles[C2])
  AlignInfo3 = read.csv(allfiles[C3])
  AlignInfo4 = read.csv(allfiles[C4])
  AlignInfo5 = read.csv(allfiles[C5])
  conditions = c(Condition1Name,Condition2Name,Condition3Name,Condition4Name,Condition5Name)
}

#AlignInfo files will be in the same order as the ConditionNames specified above

#######################################################################################################################################

############################################  ALIGNMENT ALGORITHM  ####################################################################

#######################################################################################################################################
## Reading in Conditions

for(i in 1:numberofconditions){
  
  #Grab VOC names from each condition
  file = paste("C",i,"VOC",sep="")
  dataname = paste("AlignInfo",i,sep="")
  tmp = get(dataname)
  
  #tmp[,1] = make.names(tmp[,1],unique = T)
  assign(file,tmp[,1])
  
  #Name Ref Column
  navec = rep(NA,nrow(tmp))
  navec = tmp[,1]
  
  #Grab the number of time points the VOC was observed in 
  file2 = paste("C",i,"Count",sep="")
  handle = "Number"
  assign(file2,quickfun(tmp,handle))
  #Clean up the matrix
  tmp2 = get(file2)
  tmp2 = data.frame(matrix(tmp2,ncol=1))
  tmp2 = cbind(navec,tmp2)
  colnames(tmp2) = c("VOC","Count")
  #rownames(tmp2) = tmp[,1]
  assign(file2,tmp2)
  
  #Grab the indices for each data type (Class, Area, Dimension1RT, Dimension2RT)
  I1 = paste("ClassIndex",conditions[i],sep="")
  I2 = paste("AreaIndex",conditions[i],sep="")
  I3 = paste("Dim1Index",conditions[i],sep="")
  I4 = paste("Dim2Index",conditions[i],sep="")
  tmp3 = get(dataname)
  ph = which(substr(colnames(tmp3),1,4) == "Clas")
  assign(I1,ph)
  ph = which(substr(colnames(tmp3),1,4) == "Area")
  assign(I2,ph)
  ph = which(substr(colnames(tmp3),1,4) == "Dim1")
  assign(I3,ph)
  ph = which(substr(colnames(tmp3),1,4) == "Dim2")
  assign(I4,ph)
  
  
  #Subset first dim info
  D1 = paste("DimOne",conditions[i],sep="")
  tmp4 = get(dataname)
  index = get(I3)
  tmp5 = tmp4[,index]
  vec = seq(1,nsamples,1)
  mycols = paste0(conditions[i],vec,sep="")
  tmp5 = cbind(navec,tmp5)
  #rownames(tmp5) = tmp[,1]
  colnames(tmp5) = c("VOC",mycols)
  assign(D1,tmp5)
  
  #Subset second dim info
  D2 = paste("DimTwo",conditions[i],sep="")
  tmp6 = get(dataname)
  index = get(I4)
  tmp7 = tmp6[,index]
  vec = seq(1,nsamples,1)
  mycols = paste0(conditions[i],vec,sep="")
  tmp7 = cbind(navec,tmp7)
  #rownames(tmp7) = tmp[,1]
  colnames(tmp7) = c("VOC",mycols)
  assign(D2,tmp7)
  
  #Subset abundance info
  D3 = paste(conditions[i],"Abundance",sep="")
  tmp8 = get(dataname)
  index = get(I2)
  tmp9 = tmp8[,index]
  vec = seq(1,nsamples,1)
  mycols = paste0(conditions[i],vec,sep="")
  tmp9 = cbind(navec,tmp9)
  #rownames(tmp9) = tmp[,1]
  colnames(tmp9) = c("VOC",mycols)
  assign(D3,tmp9)
}

rm(tmp,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9) #to help with system memory



#Remove VOCs that show up more than 'maxdups' times (mainly siloxanes, ethanol, and other contaminants)
for(i in 1:numberofconditions){
  dataname = paste("AlignInfo",i,sep="")
  tmp = get(dataname)
  #tmp[,1] = make.names(tmp[,1],unique = T)
  
  VOClist = table(tmp[,1])
  count = 1
  RemoveNames = rep(NA,length(VOClist))
  
  for(j in 1:length(VOClist)){
    if(VOClist[[j]] > maxdups){
      RemoveNames[count] = names(VOClist[j])
      count = count+1
    }
  }
  RemoveNames = RemoveNames[!is.na(RemoveNames)]
  
  dataname = paste(conditions[i],"Abundance",sep="")
  Atmp = get(dataname)
  
  dataname = paste("DimTwo",conditions[i],sep="")
  D2tmp = get(dataname)
  
  dataname = paste("DimOne",conditions[i],sep="")
  D1tmp = get(dataname)
  
  dataname = paste("C",i,"Count",sep="")
  Ctmp = get(dataname)
  
  Atmp2 = Atmp[! Atmp$VOC %in% RemoveNames,]
  D2tmp2 = D2tmp[ ! D2tmp$VOC %in% RemoveNames,]
  D1tmp2 = D1tmp[ ! D1tmp$VOC %in% RemoveNames,]
  Ctmp2 = Ctmp[! Ctmp$VOC %in% RemoveNames,]
  
  filt1 = paste("Clean",conditions[i],"Abundance",sep="")
  filt2 = paste("Clean",conditions[i], "DimTwo",sep="")
  filt3 = paste("Clean",conditions[i],"DimOne",sep="")
  filt4 = paste("Clean",conditions[i],"Count",sep="")
  
  assign(filt1,Atmp2)
  assign(filt2,D2tmp2)
  assign(filt3,D1tmp2)
  assign(filt4,Ctmp2)
}

rm(Atmp,Atmp2,D2tmp,D2tmp2,D1tmp,D1tmp2,Ctmp,Ctmp2)

##############################################################################################################################
## FILTERING LOW "COUNT" VOCs

#Optional: Remove VOCs not present in 80% heavy and light samples - for better looking plots (less missing data)
#This is not necessary if experienced ChromaTOF user can preselect the VOCs present in >60-80% of samples
#Code should be included so that unprocessed abundance data can be used

if(option == 1){
  thresh = FilteringThresh
  cat("Option 1 selected, filtering is being applied to all files/conditions.")
  
  for(i in 1:numberofconditions){
    dataname = paste("Clean",conditions[i],"Count",sep="")
    tmp = get(dataname)
    
    tmp$Count = tmp$Count/nsamples
    
    index = which(tmp$Count < thresh)
    
    dataname = paste("Clean",conditions[i],"Abundance",sep="")
    fdata1 = get(dataname)
    
    dataname = paste("Clean",conditions[i], "DimTwo",sep="")
    fdata2 = get(dataname)
    
    dataname = paste("Clean",conditions[i],"DimOne",sep="")
    fdata3 = get(dataname)
    
    dataname = paste("Clean",conditions[i],"Count",sep="")
    fdata4 = get(dataname)
    
    fdata1 = fdata1[-index,]
    fdata2 = fdata2[-index,]
    fdata3 = fdata3[-index,]
    fdata4 = fdata4[-index,]
    
    fin1 = paste("Final",conditions[i],"Abundance",sep="")
    fin2 = paste("Final",conditions[i], "DimTwo",sep="")
    fin3 = paste("Final",conditions[i], "DimOne",sep="")
    fin4 = paste("Final",conditions[i] , "Count",sep="")
    
    assign(fin1,fdata1)
    assign(fin2,fdata2)
    assign(fin3,fdata3)
    assign(fin4,fdata4)
  }
  
  
}else if(option == 2){
  thresh = FilteringThresh
  cat("Option 2 selected, filtering is being applied exclusively to the experimental conditions.")
  
  for(i in 1:length(expconditions)){
    
    experiment = expconditions[i]
    
    e = as.numeric(which(conditions == experiment))
    
    
    dataname = paste("Clean",conditions[e],"Count",sep="")
    tmp = get(dataname)
    
    tmp$Count = tmp$Count/nsamples
    
    index = which(tmp$Count < thresh)
    
    dataname = paste("Clean",conditions[e],"Abundance",sep="")
    fdata1 = get(dataname)
    
    dataname = paste("Clean",conditions[e], "DimTwo",sep="")
    fdata2 = get(dataname)
    
    dataname = paste("Clean",conditions[e],"DimOne",sep="")
    fdata3 = get(dataname)
    
    dataname = paste("Clean",conditions[e],"Count",sep="")
    fdata4 = get(dataname)
    
    fdata1 = fdata1[-index,]
    fdata2 = fdata2[-index,]
    fdata3 = fdata3[-index,]
    fdata4 = fdata4[-index,]
    
    fin1 = paste("Final",conditions[i],"Abundance",sep="")
    fin2 = paste("Final",conditions[i], "DimTwo",sep="")
    fin3 = paste("Final",conditions[i], "DimOne",sep="")
    fin4 = paste("Final",conditions[i] , "Count",sep="")
    
    assign(fin1,fdata1)
    assign(fin2,fdata2)
    assign(fin3,fdata3)
    assign(fin4,fdata4)
    
    
  }
  
}else if(option == 3){
  cat("Option 3 selected, no filtering is applied. Proceeding to the next step...")
}else if(option == 4){
  thresh1 = FilteringThresh1
  cat("Option 4 selected, filtering is being applied to the experimental conditions using different thresholds.")
  
  for(i in 1:length(expconditions1)){
    
    experiment = expconditions1[i]
    
    e = as.numeric(which(conditions == experiment))
    
    
    dataname = paste("Clean",conditions[e],"Count",sep="")
    tmp = get(dataname)
    
    tmp$Count = tmp$Count/nsamples
    
    index = which(tmp$Count < thresh1)
    
    dataname = paste("Clean",conditions[e],"Abundance",sep="")
    fdata1 = get(dataname)
    
    dataname = paste("Clean",conditions[e], "DimTwo",sep="")
    fdata2 = get(dataname)
    
    dataname = paste("Clean",conditions[e],"DimOne",sep="")
    fdata3 = get(dataname)
    
    dataname = paste("Clean",conditions[e],"Count",sep="")
    fdata4 = get(dataname)
    
    fdata1 = fdata1[-index,]
    fdata2 = fdata2[-index,]
    fdata3 = fdata3[-index,]
    fdata4 = fdata4[-index,]
    
    fin1 = paste("Final",conditions[i],"Abundance",sep="")
    fin2 = paste("Final",conditions[i], "DimTwo",sep="")
    fin3 = paste("Final",conditions[i], "DimOne",sep="")
    fin4 = paste("Final",conditions[i] , "Count",sep="")
    
    assign(fin1,fdata1)
    assign(fin2,fdata2)
    assign(fin3,fdata3)
    assign(fin4,fdata4)
  }
  
  
  thresh2 = FilteringThresh2
  cat("Option 2 selected, filtering is being applied to only the experimental conditions.")
  
  for(i in 1:length(expconditions2)){
    
    experiment = expconditions2[i]
    
    e = as.numeric(which(conditions == experiment))
    
    
    dataname = paste("Clean",conditions[e],"Count",sep="")
    tmp = get(dataname)
    
    tmp$Count = tmp$Count/nsamples
    
    index = which(tmp$Count < thresh2)
    
    dataname = paste("Clean",conditions[e],"Abundance",sep="")
    fdata1 = get(dataname)
    
    dataname = paste("Clean",conditions[e], "DimTwo",sep="")
    fdata2 = get(dataname)
    
    dataname = paste("Clean",conditions[e],"DimOne",sep="")
    fdata3 = get(dataname)
    
    dataname = paste("Clean",conditions[e],"Count",sep="")
    fdata4 = get(dataname)
    
    fdata1 = fdata1[-index,]
    fdata2 = fdata2[-index,]
    fdata3 = fdata3[-index,]
    fdata4 = fdata4[-index,]
    
    fin1 = paste("Final",conditions[i],"Abundance",sep="")
    fin2 = paste("Final",conditions[i], "DimTwo",sep="")
    fin3 = paste("Final",conditions[i], "DimOne",sep="")
    fin4 = paste("Final",conditions[i] , "Count",sep="")
    
    assign(fin1,fdata1)
    assign(fin2,fdata2)
    assign(fin3,fdata3)
    assign(fin4,fdata4)
  }
}else{
  cat("Please select a valid filtering option (1-4). See description at the start of the code.")
  break
}

rm(fdata1,fdata2,fdata3,fdata4,fin1,fin2,fin3,fin4)

###############################################################################################################################
# Subset into high confidence matches (for automated alignment) and low confidence matches requiring manual alignment
###############################################################################################################################
cat("Establishing parameters for downstream parsing...")
#Obtains first and second dimension averages for all files, according to the selected option
if(option == 3){
  
  for(i in 1:numberofconditions){
    dataname = paste("Clean",conditions[i],"DimOne",sep="")
    D1RT = get(dataname)
    D1RT = D1RT[,-1]
    
    dataname = paste("Clean",conditions[i],"DimTwo",sep="")
    D2RT = get(dataname)
    D2RT = D2RT[,-1]
    
    avgD1 = apply(D1RT,1,mean,na.rm=T)
    avgD2 = apply(D2RT,1,mean,na.rm=T)
    
    avgname1 = paste("AvgD1RT",conditions[i],sep="")
    avgname2 = paste("AvgD2RT",conditions[i],sep="")
    
    assign(avgname1,avgD1)
    assign(avgname2,avgD2)
  }
  
  
}else if(option == 1){

  for(i in 1:numberofconditions){
    dataname = paste("Final",conditions[e], "DimOne",sep="")
    D1RT = get(dataname)
    D1RT = D1RT[,-1]
    
    dataname = paste("Final",conditions[e], "DimTwo",sep="")
    D2RT = get(dataname)
    D2RT = D2RT[,-1]
  
    avgD1 = apply(D1RT,1,mean,na.rm=T)
    avgD2 = apply(D2RT,1,mean,na.rm=T)
  
    avgname1 = paste("AvgD1RT",conditions[i],sep="")
    avgname2 = paste("AvgD2RT",conditions[i],sep="")
  
    assign(avgname1,avgD1)
    assign(avgname2,avgD2)
  }
}else if(option == 2){
  
  for(i in 1:length(expconditions)){
    
    experiment = expconditions[i]
    
    e = as.numeric(which(conditions == experiment))
    
    dataname = paste("Final",conditions[e], "DimOne",sep="")
    D1RT = get(dataname)
    D1RT = D1RT[,-1]
    
    dataname = paste("Final",conditions[e], "DimTwo",sep="")
    D2RT = get(dataname)
    D2RT = D2RT[,-1]
    
    avgD1 = apply(D1RT,1,mean,na.rm=T)
    avgD2 = apply(D2RT,1,mean,na.rm=T)
    
    avgname1 = paste("AvgD1RT",conditions[e],sep="")
    avgname2 = paste("AvgD2RT",conditions[e],sep="")
    
    assign(avgname1,avgD1)
    assign(avgname2,avgD2)
  }
  
}else if(option == 4){
  for(i in 1:length(expconditions1)){
    
    experiment = expconditions1[i]
    
    e = as.numeric(which(conditions == experiment))
    
    dataname = paste("Final",conditions[e], "DimOne",sep="")
    D1RT = get(dataname)
    D1RT = D1RT[,-1]
    
    dataname = paste("Final",conditions[e], "DimTwo",sep="")
    D2RT = get(dataname)
    D2RT = D2RT[,-1]
    
    avgD1 = apply(D1RT,1,mean,na.rm=T)
    avgD2 = apply(D2RT,1,mean,na.rm=T)
    
    avgname1 = paste("AvgD1RT",conditions[e],sep="")
    avgname2 = paste("AvgD2RT",conditions[e],sep="")
    
    assign(avgname1,avgD1)
    assign(avgname2,avgD2)
  }
  
  for(i in 1:length(expconditions2)){
    
    experiment = expconditions2[i]
    
    e = as.numeric(which(conditions == experiment))
    
    dataname = paste("Final",conditions[e], "DimOne",sep="")
    D1RT = get(dataname)
    D1RT = D1RT[,-1]
    
    dataname = paste("Final",conditions[e], "DimTwo",sep="")
    D2RT = get(dataname)
    D2RT = D2RT[,-1]
    
    avgD1 = apply(D1RT,1,mean,na.rm=T)
    avgD2 = apply(D2RT,1,mean,na.rm=T)
    
    avgname1 = paste("AvgD1RT",conditions[e],sep="")
    avgname2 = paste("AvgD2RT",conditions[e],sep="")
    
    assign(avgname1,avgD1)
    assign(avgname2,avgD2)
  }
}

rm(avgD1,avgD2,D1RT,D2RT)

###############################################################################################################################
# Parse 
###############################################################################################################################
cat("Parsing low and high quality VOCs...")
if(option == 1){
  for(i in 1:numberofconditions){
    
    dataname = paste("Final",conditions[i], "DimOne",sep="")
    D1RT = get(dataname)
    dataname = paste("AvgD1RT",conditions[i],sep="")
    AvgD1RT = get(dataname)
    logicalindexD1 = data.frame(matrix(NA,nrow(D1RT),ncol=nsamples+1))
    colnames(logicalindexD1) = colnames(D1RT)
    logicalindexD1$VOC = D1RT$VOC
    
    dataname = paste("Final",conditions[i], "DimTwo",sep="")
    D2RT = get(dataname)
    dataname = paste("AvgD2RT",conditions[i],sep="")
    AvgD2RT = get(dataname)
    logicalindexD2 = data.frame(matrix(NA,nrow(D2RT),ncol=nsamples+1))
    colnames(logicalindexD2) = colnames(D2RT)
    logicalindexD2$VOC = D2RT$VOC
    
    for(j in 1:length(AvgD1RT)){
      for(k in 2:(nsamples+1)){
        if(is.na(D1RT[j,k])){
          logicalindexD1[j,k] = 0 #'High' Quality
        }else{
          if(abs(D1RT[j,k] - AvgD1RT[j]) > Dim1Tol){
            logicalindexD1[j,k] = 0 #Low Quality
          }else{
            logicalindexD1[j,k] = 1 #High Quality
          }
        }
      }
    }
    
    for(j in 1:length(AvgD2RT)){
      for(k in 2:(nsamples+1)){
        if(is.na(D2RT[j,k])){
          logicalindexD2[j,k] = 0 #'High' Quality
        }else{
          if(abs(D2RT[j,k] - AvgD2RT[j]) > Dim2Tol){
            logicalindexD2[j,k] = 0 #Low Quality
          }else{
            logicalindexD2[j,k] = 1 #High Quality
          }
        }
      }
    }
    
    tmp1 = logicalindexD1[,-1]
    tmp2 = logicalindexD2[,-1]
    LID1 = apply(tmp1,1,sum,na.rm=T)
    LID2 = apply(tmp2,1,sum,na.rm=T)
    
    HQ1 = rep(NA,length(LID1))
    HQ2 = rep(NA,length(LID2))
    
    dataname = paste("Final",conditions[i] , "Count",sep="")
    tmpcount = get(dataname)
    
    for(m in 1:length(LID1)){
      if(LID1[m] == tmpcount$Count[m]){
        HQ1[m] = m
      }
    }
    
    for(m in 1:length(LID2)){
      if(LID2[m] == tmpcount$Count[m]){
        HQ2[m] = m
      }
    }
    
    HQ = HQ1[HQ1 %in% HQ2] #This is the final list of indices associated with high quality VOCs
    #Also removes NAs
    
    dataname = paste("Final",conditions[i], "DimOne",sep="")
    FD1 = get(dataname)
    FD1.HQ = FD1[HQ,]
    FD1.LQ = FD1[-HQ,]
    HQnamed1 = paste("HighQuality",conditions[i],"DimOne",sep="")
    LQnamed1 = paste("LowQuality",conditions[i],"DimOne",sep="")
    assign(HQnamed1,FD1.HQ)
    assign(LQnamed1,FD1.LQ)
    
    dataname = paste("Final",conditions[i], "DimTwo",sep="")
    FD2 = get(dataname)
    FD2.HQ = FD2[HQ,]
    FD2.LQ = FD2[-HQ,]
    HQnamed2 = paste("HighQuality",conditions[i],"DimTwo",sep="")
    LQnamed2 = paste("LowQuality",conditions[i],"DimTwo",sep="")
    assign(HQnamed2,FD2.HQ)
    assign(LQnamed2,FD2.LQ)
    
    dataname = paste("Final",conditions[i], "Abundance",sep="")
    FA = get(dataname)
    FA.HQ = FA[HQ,]
    FA.LQ = FA[-HQ,]
    HQnameab = paste("HighQuality",conditions[i],"Abundance",sep="")
    LQnameab = paste("LowQuality",conditions[i],"Abundance",sep="")
    assign(HQnameab,FA.HQ)
    assign(LQnameab,FA.LQ)
    
    dataname = paste("Final",conditions[i], "Count",sep="")
    FC = get(dataname)
    FC.HQ = FC[HQ,]
    FC.LQ = FC[-HQ,]
    HQnameco = paste("HighQuality",conditions[i],"Count",sep="")
    LQnameco = paste("LowQuality",conditions[i],"Count",sep="")
    assign(HQnameco,FC.HQ)
    assign(LQnameco,FC.LQ)
    
  }
  
}else if(option == 2){
  
  for(i in 1:length(expconditions)){
    
    experiment = expconditions[i]
    e = as.numeric(which(conditions == experiment))
    
    dataname = paste("Final",conditions[e], "DimOne",sep="")
    D1RT = get(dataname)
    dataname = paste("AvgD1RT",conditions[e],sep="")
    AvgD1RT = get(dataname)
    logicalindexD1 = data.frame(matrix(NA,nrow(D1RT),ncol=nsamples+1))
    colnames(logicalindexD1) = colnames(D1RT)
    logicalindexD1$VOC = D1RT$VOC
    
    dataname = paste("Final",conditions[e], "DimTwo",sep="")
    D2RT = get(dataname)
    dataname = paste("AvgD2RT",conditions[e],sep="")
    AvgD2RT = get(dataname)
    logicalindexD2 = data.frame(matrix(NA,nrow(D2RT),ncol=nsamples+1))
    colnames(logicalindexD2) = colnames(D2RT)
    logicalindexD2$VOC = D2RT$VOC
    
    for(j in 1:length(AvgD1RT)){
      for(k in 2:(nsamples+1)){
        if(is.na(D1RT[j,k])){
          logicalindexD1[j,k] = 0#'High' Quality
        }else{
          if(abs(D1RT[j,k] - AvgD1RT[j]) > Dim1Tol){
            logicalindexD1[j,k] = 0 #Low Quality
          }else{
            logicalindexD1[j,k] = 1 #High Quality
          }
        }
      }
    }
    
    for(j in 1:length(AvgD2RT)){
      for(k in 2:(nsamples+1)){
        if(is.na(D2RT[j,k])){
          logicalindexD2[j,k] = 0 #'High' Quality
        }else{
          if(abs(D2RT[j,k] - AvgD2RT[j]) > Dim2Tol){
            logicalindexD2[j,k] = 0 #Low Quality
          }else{
            logicalindexD2[j,k] = 1 #High Quality
          }
        }
      }
    }
    
    tmp1 = logicalindexD1[,-1]
    tmp2 = logicalindexD2[,-1]
    LID1 = apply(tmp1,1,sum,na.rm=T)
    LID2 = apply(tmp2,1,sum,na.rm=T)
    
    HQ1 = rep(NA,length(LID1))
    HQ2 = rep(NA,length(LID2))
    
    dataname = paste("Final",conditions[e] , "Count",sep="")
    tmpcount = get(dataname)
    
    for(m in 1:length(LID1)){
      if(LID1[m] == tmpcount$Count[m]){
        HQ1[m] = m
      }
    }
    
    for(m in 1:length(LID2)){
      if(LID2[m] == tmpcount$Count[m]){
        HQ2[m] = m
      }
    }
    
    HQ = HQ1[HQ1 %in% HQ2] #This is the final list of indices associated with high quality VOCs
    
    dataname = paste("Final",conditions[e], "DimOne",sep="")
    FD1 = get(dataname)
    FD1.HQ = FD1[HQ,]
    FD1.LQ = FD1[-HQ,]
    HQnamed1 = paste("HighQuality",conditions[e],"DimOne",sep="")
    LQnamed1 = paste("LowQuality",conditions[e],"DimOne",sep="")
    assign(HQnamed1,FD1.HQ)
    assign(LQnamed1,FD1.LQ)
    
    dataname = paste("Final",conditions[e], "DimTwo",sep="")
    FD2 = get(dataname)
    FD2.HQ = FD2[HQ,]
    FD2.LQ = FD2[-HQ,]
    HQnamed2 = paste("HighQuality",conditions[e],"DimTwo",sep="")
    LQnamed2 = paste("LowQuality",conditions[e],"DimTwo",sep="")
    assign(HQnamed2,FD2.HQ)
    assign(LQnamed2,FD2.LQ)
    
    dataname = paste("Final",conditions[e], "Abundance",sep="")
    FA = get(dataname)
    FA.HQ = FA[HQ,]
    FA.LQ = FA[-HQ,]
    HQnameab = paste("HighQuality",conditions[e],"Abundance",sep="")
    LQnameab = paste("LowQuality",conditions[e],"Abundance",sep="")
    assign(HQnameab,FA.HQ)
    assign(LQnameab,FA.LQ)
    
    dataname = paste("Final",conditions[e], "Count",sep="")
    FC = get(dataname)
    FC.HQ = FC[HQ,]
    FC.LQ = FC[-HQ,]
    HQnameco = paste("HighQuality",conditions[e],"Count",sep="")
    LQnameco = paste("LowQuality",conditions[e],"Count",sep="")
    assign(HQnameco,FC.HQ)
    assign(LQnameco,FC.LQ)
    
  }
  
  
}else if(option == 3){
  
  for(i in 1:numberofconditions){
    
    dataname = paste("Clean",conditions[i],"DimOne",sep="")
    D1RT = get(dataname)
    dataname = paste("AvgD1RT",conditions[i],sep="")
    AvgD1RT = get(dataname)
    logicalindexD1 = data.frame(matrix(NA,nrow(D1RT),ncol=nsamples+1))
    colnames(logicalindexD1) = colnames(D1RT)
    logicalindexD1$VOC = D1RT$VOC
    
    dataname = paste("Clean",conditions[i],"DimTwo",sep="")
    D2RT = get(dataname)
    dataname = paste("AvgD2RT",conditions[i],sep="")
    AvgD2RT = get(dataname)
    logicalindexD2 = data.frame(matrix(NA,nrow(D2RT),ncol=nsamples+1))
    colnames(logicalindexD2) = colnames(D2RT)
    logicalindexD2$VOC = D2RT$VOC
    
    for(j in 1:length(AvgD1RT)){
      for(k in 2:(nsamples+1)){
        if(is.na(D1RT[j,k])){
          logicalindexD1[j,k] = 0 #'High' Quality
        }else{
          if(abs(D1RT[j,k] - AvgD1RT[j]) > Dim1Tol){
            logicalindexD1[j,k] = 0 #Low Quality
          }else{
            logicalindexD1[j,k] = 1 #High Quality
          }
        }
      }
    }
    
    for(j in 1:length(AvgD2RT)){
      for(k in 2:(nsamples+1)){
        if(is.na(D2RT[j,k])){
          logicalindexD2[j,k] = 0 #'High' Quality
        }else{
          if(abs(D2RT[j,k] - AvgD2RT[j]) > Dim2Tol){
            logicalindexD2[j,k] = 0 #Low Quality
          }else{
            logicalindexD2[j,k] = 1 #High Quality
          }
        }
      }
    }
    
    tmp1 = logicalindexD1[,-1]
    tmp2 = logicalindexD2[,-1]
    LID1 = apply(tmp1,1,sum,na.rm=T)
    LID2 = apply(tmp2,1,sum,na.rm=T)
    
    HQ1 = rep(NA,length(LID1))
    HQ2 = rep(NA,length(LID2))
    
    dataname = paste("Final",conditions[i] , "Count",sep="")
    tmpcount = get(dataname)
    
    for(m in 1:length(LID1)){
      if(LID1[m] == tmpcount$Count[m]){
        HQ1[m] = m
      }
    }
    
    for(m in 1:length(LID2)){
      if(LID2[m] == tmpcount$Count[m]){
        HQ2[m] = m
      }
    }
    
    HQ = HQ1[HQ1 %in% HQ2] #This is the final list of indices associated with high quality VOCs
    
    dataname = paste("Final",conditions[i], "DimOne",sep="")
    FD1 = get(dataname)
    FD1.HQ = FD1[HQ,]
    FD1.LQ = FD1[-HQ,]
    HQnamed1 = paste("HighQuality",conditions[i],"DimOne",sep="")
    LQnamed1 = paste("LowQuality",conditions[i],"DimOne",sep="")
    assign(HQnamed1,FD1.HQ)
    assign(LQnamed1,FD1.LQ)
    
    dataname = paste("Final",conditions[i], "DimTwo",sep="")
    FD2 = get(dataname)
    FD2.HQ = FD2[HQ,]
    FD2.LQ = FD2[-HQ,]
    HQnamed2 = paste("HighQuality",conditions[i],"DimTwo",sep="")
    LQnamed2 = paste("LowQuality",conditions[i],"DimTwo",sep="")
    assign(HQnamed2,FD2.HQ)
    assign(LQnamed2,FD2.LQ)
    
    dataname = paste("Final",conditions[i], "Abundance",sep="")
    FA = get(dataname)
    FA.HQ = FA[HQ,]
    FA.LQ = FA[-HQ,]
    HQnameab = paste("HighQuality",conditions[i],"Abundance",sep="")
    LQnameab = paste("LowQuality",conditions[i],"Abundance",sep="")
    assign(HQnameab,FA.HQ)
    assign(LQnameab,FA.LQ)
    
    dataname = paste("Final",conditions[i], "Count",sep="")
    FC = get(dataname)
    FC.HQ = FC[HQ,]
    FC.LQ = FC[-HQ,]
    HQnameco = paste("HighQuality",conditions[i],"Count",sep="")
    LQnameco = paste("LowQuality",conditions[i],"Count",sep="")
    assign(HQnameco,FC.HQ)
    assign(LQnameco,FC.LQ)
    
  }
  
}else if(option == 4){
  
  for(i in 1:length(expconditions1)){
    
    experiment = expconditions1[i]
    e = as.numeric(which(conditions == experiment))
    
    dataname = paste("Final",conditions[e], "DimOne",sep="")
    D1RT = get(dataname)
    dataname = paste("AvgD1RT",conditions[e],sep="")
    AvgD1RT = get(dataname)
    logicalindexD1 = data.frame(matrix(NA,nrow(D1RT),ncol=nsamples+1))
    colnames(logicalindexD1) = colnames(D1RT)
    logicalindexD1$VOC = D1RT$VOC
    
    dataname = paste("Final",conditions[e], "DimTwo",sep="")
    D2RT = get(dataname)
    dataname = paste("AvgD2RT",conditions[e],sep="")
    AvgD2RT = get(dataname)
    logicalindexD2 = data.frame(matrix(NA,nrow(D2RT),ncol=nsamples+1))
    colnames(logicalindexD2) = colnames(D2RT)
    logicalindexD2$VOC = D2RT$VOC
    
    for(j in 1:length(AvgD1RT)){
      for(k in 2:(nsamples+1)){
        if(is.na(D1RT[j,k])){
          logicalindexD1[j,k] = 0 #'High' Quality
        }else{
          if(abs(D1RT[j,k] - AvgD1RT[j]) > Dim1Tol){
            logicalindexD1[j,k] = 0 #Low Quality
          }else{
            logicalindexD1[j,k] = 1 #High Quality
          }
        }
      }
    }
    
    for(j in 1:length(AvgD2RT)){
      for(k in 2:(nsamples+1)){
        if(is.na(D2RT[j,k])){
          logicalindexD2[j,k] = 0 #'High' Quality
        }else{
          if(abs(D2RT[j,k] - AvgD2RT[j]) > Dim2Tol){
            logicalindexD2[j,k] = 0 #Low Quality
          }else{
            logicalindexD2[j,k] = 1 #High Quality
          }
        }
      }
    }
    
    tmp1 = logicalindexD1[,-1]
    tmp2 = logicalindexD2[,-1]
    LID1 = apply(tmp1,1,sum,na.rm=T)
    LID2 = apply(tmp2,1,sum,na.rm=T)
    
    HQ1 = rep(NA,length(LID1))
    HQ2 = rep(NA,length(LID2))
    
    dataname = paste("Final",conditions[e] , "Count",sep="")
    tmpcount = get(dataname)
    
    for(m in 1:length(LID1)){
      if(LID1[m] == tmpcount$Count[m]){
        HQ1[m] = m
      }
    }
    
    for(m in 1:length(LID2)){
      if(LID2[m] == tmpcount$Count[m]){
        HQ2[m] = m
      }
    }
    
    HQ = HQ1[HQ1 %in% HQ2] #This is the final list of indices associated with high quality VOCs
    
    dataname = paste("Final",conditions[e], "DimOne",sep="")
    FD1 = get(dataname)
    FD1.HQ = FD1[HQ,]
    FD1.LQ = FD1[-HQ,]
    HQnamed1 = paste("HighQuality",conditions[e],"DimOne",sep="")
    LQnamed1 = paste("LowQuality",conditions[e],"DimOne",sep="")
    assign(HQnamed1,FD1.HQ)
    assign(LQnamed1,FD1.LQ)
    
    dataname = paste("Final",conditions[e], "DimTwo",sep="")
    FD2 = get(dataname)
    FD2.HQ = FD2[HQ,]
    FD2.LQ = FD2[-HQ,]
    HQnamed2 = paste("HighQuality",conditions[e],"DimTwo",sep="")
    LQnamed2 = paste("LowQuality",conditions[e],"DimTwo",sep="")
    assign(HQnamed2,FD2.HQ)
    assign(LQnamed2,FD2.LQ)
    
    dataname = paste("Final",conditions[e], "Abundance",sep="")
    FA = get(dataname)
    FA.HQ = FA[HQ,]
    FA.LQ = FA[-HQ,]
    HQnameab = paste("HighQuality",conditions[e],"Abundance",sep="")
    LQnameab = paste("LowQuality",conditions[e],"Abundance",sep="")
    assign(HQnameab,FA.HQ)
    assign(LQnameab,FA.LQ)
    
    dataname = paste("Final",conditions[e], "Count",sep="")
    FC = get(dataname)
    FC.HQ = FC[HQ,]
    FC.LQ = FC[-HQ,]
    HQnameco = paste("HighQuality",conditions[e],"Count",sep="")
    LQnameco = paste("LowQuality",conditions[e],"Count",sep="")
    assign(HQnameco,FC.HQ)
    assign(LQnameco,FC.LQ)
    
  }
  
  for(i in 1:length(expconditions2)){
    
    experiment = expconditions2[i]
    e = as.numeric(which(conditions == experiment))
    
    dataname = paste("Final",conditions[e], "DimOne",sep="")
    D1RT = get(dataname)
    dataname = paste("AvgD1RT",conditions[e],sep="")
    AvgD1RT = get(dataname)
    logicalindexD1 = data.frame(matrix(NA,nrow(D1RT),ncol=nsamples+1))
    colnames(logicalindexD1) = colnames(D1RT)
    logicalindexD1$VOC = D1RT$VOC
    
    dataname = paste("Final",conditions[e], "DimTwo",sep="")
    D2RT = get(dataname)
    dataname = paste("AvgD2RT",conditions[e],sep="")
    AvgD2RT = get(dataname)
    logicalindexD2 = data.frame(matrix(NA,nrow(D2RT),ncol=nsamples+1))
    colnames(logicalindexD2) = colnames(D2RT)
    logicalindexD2$VOC = D2RT$VOC
    
    for(j in 1:length(AvgD1RT)){
      for(k in 2:(nsamples+1)){
        if(is.na(D1RT[j,k])){
          logicalindexD1[j,k] = 0 #'High' Quality
        }else{
          if(abs(D1RT[j,k] - AvgD1RT[j]) > Dim1Tol){
            logicalindexD1[j,k] = 0 #Low Quality
          }else{
            logicalindexD1[j,k] = 1 #High Quality
          }
        }
      }
    }
    
    for(j in 1:length(AvgD2RT)){
      for(k in 2:(nsamples+1)){
        if(is.na(D2RT[j,k])){
          logicalindexD2[j,k] = 0 #'High' Quality
        }else{
          if(abs(D2RT[j,k] - AvgD2RT[j]) > Dim2Tol){
            logicalindexD2[j,k] = 0 #Low Quality
          }else{
            logicalindexD2[j,k] = 1 #High Quality
          }
        }
      }
    }
    
    tmp1 = logicalindexD1[,-1]
    tmp2 = logicalindexD2[,-1]
    LID1 = apply(tmp1,1,sum,na.rm=T)
    LID2 = apply(tmp2,1,sum,na.rm=T)
    
    HQ1 = rep(NA,length(LID1))
    HQ2 = rep(NA,length(LID2))
    
    dataname = paste("Final",conditions[e] , "Count",sep="")
    tmpcount = get(dataname)
    
    for(m in 1:length(LID1)){
      if(LID1[m] == tmpcount$Count[m]){
        HQ1[m] = m
      }
    }
    
    for(m in 1:length(LID2)){
      if(LID2[m] == tmpcount$Count[m]){
        HQ2[m] = m
      }
    }
    
    HQ = HQ1[HQ1 %in% HQ2] #This is the final list of indices associated with high quality VOCs
    
    dataname = paste("Final",conditions[e], "DimOne",sep="")
    FD1 = get(dataname)
    FD1.HQ = FD1[HQ,]
    FD1.LQ = FD1[-HQ,]
    HQnamed1 = paste("HighQuality",conditions[e],"DimOne",sep="")
    LQnamed1 = paste("LowQuality",conditions[e],"DimOne",sep="")
    assign(HQnamed1,FD1.HQ)
    assign(LQnamed1,FD1.LQ)
    
    dataname = paste("Final",conditions[e], "DimTwo",sep="")
    FD2 = get(dataname)
    FD2.HQ = FD2[HQ,]
    FD2.LQ = FD2[-HQ,]
    HQnamed2 = paste("HighQuality",conditions[e],"DimTwo",sep="")
    LQnamed2 = paste("LowQuality",conditions[e],"DimTwo",sep="")
    assign(HQnamed2,FD2.HQ)
    assign(LQnamed2,FD2.LQ)
    
    dataname = paste("Final",conditions[e], "Abundance",sep="")
    FA = get(dataname)
    FA.HQ = FA[HQ,]
    FA.LQ = FA[-HQ,]
    HQnameab = paste("HighQuality",conditions[e],"Abundance",sep="")
    LQnameab = paste("LowQuality",conditions[e],"Abundance",sep="")
    assign(HQnameab,FA.HQ)
    assign(LQnameab,FA.LQ)
    
    dataname = paste("Final",conditions[e], "Count",sep="")
    FC = get(dataname)
    FC.HQ = FC[HQ,]
    FC.LQ = FC[-HQ,]
    HQnameco = paste("HighQuality",conditions[e],"Count",sep="")
    LQnameco = paste("LowQuality",conditions[e],"Count",sep="")
    assign(HQnameco,FC.HQ)
    assign(LQnameco,FC.LQ)
    
  }
  
}else{
  cat("Please select a valid filtering option (1-4). See description at the start of the code.")
  break
}

######## Code is fully automated to this point 8/18/21


###############################################################################################################################
# Remove siloxanes and TMS (trimethylsily)
###############################################################################################################################

cat("Filtering contaminants...")

#High Quality
for(i in 1:length(expconditions)){
  experiment = expconditions[i]
  e = as.numeric(which(conditions == experiment))
  
  dataname = paste("HighQuality",conditions[e],"DimOne",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "sil"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  
  dataname = paste("HighQuality",conditions[e], "DimTwo",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "sil"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  
  dataname = paste("HighQuality",conditions[e],"Abundance",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "sil"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  
  dataname = paste("HighQuality",conditions[e],"Count",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "sil"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
}

for(i in 1:length(expconditions)){
  experiment = expconditions[i]
  e = as.numeric(which(conditions == experiment))
  
  dataname = paste("HighQuality",conditions[e],"DimOne",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "TMS"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  
  dataname = paste("HighQuality",conditions[e], "DimTwo",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "TMS"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  
  dataname = paste("HighQuality",conditions[e],"Abundance",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "TMS"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  
  dataname = paste("HighQuality",conditions[e],"Count",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "TMS"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
}


#Low Quality
for(i in 1:length(expconditions)){
  experiment = expconditions[i]
  e = as.numeric(which(conditions == experiment))
  
  dataname = paste("LowQuality",conditions[e],"DimOne",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  if(nrow(tmp) != 0){
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "sil"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  }
  
  dataname = paste("LowQuality",conditions[e], "DimTwo",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  if(nrow(tmp) != 0){
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "sil"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  }
  
  dataname = paste("LowQuality",conditions[e],"Abundance",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  if(nrow(tmp) != 0){
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "sil"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  }
  
  dataname = paste("LowQuality",conditions[e],"Count",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  if(nrow(tmp) != 0){
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "sil"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  }
}

for(i in 1:length(expconditions)){
  experiment = expconditions[i]
  e = as.numeric(which(conditions == experiment))
  
  dataname = paste("LowQuality",conditions[e],"DimOne",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  if(nrow(tmp) != 0){
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "TMS"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  }
  
  dataname = paste("LowQuality",conditions[e], "DimTwo",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  if(nrow(tmp) != 0){
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "TMS"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  }
  
  dataname = paste("LowQuality",conditions[e],"Abundance",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  if(nrow(tmp) != 0){
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "TMS"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  }
  
  dataname = paste("LowQuality",conditions[e],"Count",sep="")
  tmp = get(dataname)
  removeindex = rep(NA,nrow(tmp))
  if(nrow(tmp) != 0){
  for(j in 1:nrow(tmp)){
    current = tmp$VOC[j]
    b = nchar(current)
    for(k in 1:b-2){
      a = 0+k
      if(str_sub(current,a,a+2) == "TMS"){
        removeindex[j] = j
      }
    }
  }
  removeindex = removeindex[!is.na(removeindex)]
  if(length(removeindex)>0){
    tmp = tmp[-removeindex,]
    assign(dataname,tmp)
  }
  }
  
}



###############################################################################################################################
# Retention Time Driven Alignment
###############################################################################################################################

#It is recommended to focus on 'High Quality' VOCs, unless you have a reason to look through the low quality ones
#This code currently only supports alignment of high quality VOCs - 5-10-22

#Allow for many matches between heavy and light - You still must manually pull up the mass spectrums and validate similarity
#Unique Mass needs to be integrated into this code

# CURRENTLY USES: DIM1, DIM2 - 5-10-22

### Automated Alignment - Conditions 1 and 2 ONLY 

#Define size of the logical index container
n1 = paste("HighQuality",ref,"Abundance",sep="") #Dependent on reference level specified by the user
c1 = get(n1)
nr = nrow(c1)
if(length(expconditions) == 2){
  n2 = paste("HighQuality",expconditions[! expconditions %in% ref],"Abundance",sep="")
  c2 = get(n2)
  nc = nrow(c2)
}else if(length(expconditions) > 2){
  n2 = paste("HighQuality",Condition2Name,"Abundance",sep="") #The code will automatically select the second listed condition if more than 2 are available
  c2 = get(n2)
  nc = nrow(c2)
}

#Get reference condition for alignment
refnam = paste("HighQuality",ref,"DimOne",sep="")
DIMONE = get(refnam)
refnam = paste("HighQuality",ref,"DimTwo",sep="")
DIMTWO = get(refnam)

#Condition 2 retention times

if(length(expconditions) == 2){
  matchname = paste("HighQuality",expconditions[! expconditions %in% ref],"DimOne",sep="")
  MatchOne = get(matchname)
  matchname = paste("HighQuality",expconditions[! expconditions %in% ref],"DimTwo",sep="")
  MatchTwo = get(matchname)
}else if(length(expconditions) > 2){
  matchname = paste("HighQuality",Condition2Name,"DimOne",sep="") #The code will automatically select the second listed condition if more than 2 are available
  MatchOne = get(matchname)
  matchname = paste("HighQuality",Condition2Name,"DimTwo",sep="") #The code will automatically select the second listed condition if more than 2 are available
  MatchTwo = get(matchname)
}

#Necessary cleanup
rownames(DIMONE) = make.unique(DIMONE$VOC)
DIMONE = DIMONE[,-1]
rownames(DIMTWO) = make.unique(DIMTWO$VOC)
DIMTWO = DIMTWO[,-1]
rownames(MatchOne) = make.unique(MatchOne$VOC)
MatchOne = MatchOne[,-1]
rownames(MatchTwo) = make.unique(MatchTwo$VOC)
MatchTwo = MatchTwo[,-1]

cat("Aligning conditions 1 and 2...")
logindex = matrix(NA,nrow = nr,ncol=nc) #logical index container
for(i in 1:nr){
  for(j in 1:nc){ 
    
    #Set first and second dimension tolerance values for matching
    maxtol1 = max(DIMONE[i,],na.rm = T)+Dim1Tol
    mintol1 = min(DIMONE[i,],na.rm = T)-Dim1Tol
    maxtol2 = max(DIMTWO[i,],na.rm = T)+Dim2Tol
    mintol2 = min(DIMTWO[i,],na.rm = T)-Dim2Tol
    
    tmp1 = MatchOne[j,]
    tmp2 = MatchTwo[j,]
    
    #Check that VOC first and second dim retention times fall within the user defined tolerance
    if(min(tmp1,na.rm = T) > mintol1 && max(tmp1,na.rm = T) < maxtol1 && min(tmp2,na.rm = T) > mintol2 && max(tmp2,na.rm = T) < maxtol2){
      
      logindex[i,j] = 1
    }else{
      logindex[i,j] = 0
    }
  }
}


matchmatrixsize = max(c(sumrow(logindex),sumcol(logindex)))

matchindex = matrix(NA,nrow=nrow(logindex),ncol=matchmatrixsize) #number of rows corresponds to the number of VOCs in the reference condition
rownames(matchindex) = rownames(DIMONE)
colnames(matchindex) = paste("Condition2Match",seq(1,matchmatrixsize,1))
for(i in 1:nrow(logindex)){
  ti = length(which(logindex[i,] == 1))
  if(ti > 0){
    matchindex[i,1:ti] = which(logindex[i,]==1)
  }else{
    matchindex[i,1:ti+1] = NA
  }
}

matchindex = cbind(rep(NA,nrow(matchindex)),matchindex)
colnames(matchindex)[1] = paste(ref,"Index",sep="")
for(i in 1:nr){
  for(j in 1:nr){
    if(rownames(matchindex)[i] == rownames(DIMONE)[j]){
      matchindex[i,1] = j
    }
  }
}


#This code has only been verified to work using option 2

if(numberofconditions == 3){
  if(option == 2){
    
    cat("Aligning condition 3...")
    
    #Define size of the logical index container
    n1 = paste("HighQuality",ref,"Abundance",sep="") #Dependent on reference level specified by the user
    c1 = get(n1)
    nr = nrow(c1)
    n2 = paste("Clean",Condition3Name,"Abundance",sep="") #The code will automatically select the second listed condition if more than 2 are available
    c2 = get(n2)
    nc = nrow(c2)
    
    #Get reference condition for alignment
    refnam = paste("HighQuality",ref,"DimOne",sep="")
    DIMONE = get(refnam)
    refnam = paste("HighQuality",ref,"DimTwo",sep="")
    DIMTWO = get(refnam)
    
    #Condition 3 retention times
    matchname = paste("Clean",Condition3Name,"DimOne",sep="")
    MatchOne = get(matchname)
    matchname = paste("Clean",Condition3Name,"DimTwo",sep="")
    MatchTwo = get(matchname)
    
    #Necessary cleanup
    rownames(DIMONE) = make.unique(DIMONE$VOC)
    DIMONE = DIMONE[,-1]
    rownames(DIMTWO) = make.unique(DIMTWO$VOC)
    DIMTWO = DIMTWO[,-1]
    rownames(MatchOne) = make.unique(MatchOne$VOC)
    MatchOne = MatchOne[,-1]
    rownames(MatchTwo) = make.unique(MatchTwo$VOC)
    MatchTwo = MatchTwo[,-1]
    
    
    logindex = matrix(NA,nrow = nr,ncol=nc) #logical index container
    for(i in 1:nr){
      for(j in 1:nc){ 
        
        #Set first and second dimension tolerance values for matching
        maxtol1 = max(DIMONE[i,],na.rm = T)+Dim1Tol
        mintol1 = min(DIMONE[i,],na.rm = T)-Dim1Tol
        maxtol2 = max(DIMTWO[i,],na.rm = T)+Dim2Tol
        mintol2 = min(DIMTWO[i,],na.rm = T)-Dim2Tol
        
        tmp1 = MatchOne[j,]
        tmp2 = MatchTwo[j,]
        
        #Check that VOC first and second dim retention times fall within the user defined tolerance
        if(min(tmp1,na.rm = T) > mintol1 && max(tmp1,na.rm = T) < maxtol1 && min(tmp2,na.rm = T) > mintol2 && max(tmp2,na.rm = T) < maxtol2){
          
          logindex[i,j] = 1
        }else{
          logindex[i,j] = 0
        }
      }
    }
    
    
    matchmatrixsize = max(c(sumrow(logindex),sumcol(logindex)))
    
    matchindex2 = matrix(NA,nrow=nrow(logindex),ncol=matchmatrixsize) #number of rows corresponds to the number of VOCs in the reference condition
    rownames(matchindex2) = rownames(DIMONE)
    colnames(matchindex2) = paste("Condition3Match",seq(1,matchmatrixsize,1))
    for(i in 1:nrow(logindex)){
      ti = length(which(logindex[i,] == 1))
      if(ti > 0){
        matchindex2[i,1:ti] = which(logindex[i,]==1)
      }else{
        matchindex2[i,1:ti+1] = NA
      }
    }
    
    FullMatchMatrix = cbind(matchindex,matchindex2)
    View(FullMatchMatrix)
    
  }else{
    cat("Sorry, this code does not currently support this analysis")
  }
  
}else if(numberofconditions == 4){
  if(option == 2){
    
    cat("Aligning conditions 3 and 4...")
    
    #Define size of the logical index container
    n1 = paste("HighQuality",ref,"Abundance",sep="") #Dependent on reference level specified by the user
    c1 = get(n1)
    nr = nrow(c1)
    n2 = paste("Clean",Condition3Name,"Abundance",sep="") #The code will automatically select the second listed condition if more than 2 are available
    c2 = get(n2)
    nc = nrow(c2)
    
    #Get reference condition for alignment
    refnam = paste("HighQuality",ref,"DimOne",sep="")
    DIMONE = get(refnam)
    refnam = paste("HighQuality",ref,"DimTwo",sep="")
    DIMTWO = get(refnam)
    
    #Condition 3 retention times
    matchname = paste("Clean",Condition3Name,"DimOne",sep="")
    MatchOne = get(matchname)
    matchname = paste("Clean",Condition3Name,"DimTwo",sep="")
    MatchTwo = get(matchname)
    
    #Necessary cleanup
    rownames(DIMONE) = make.unique(DIMONE$VOC)
    DIMONE = DIMONE[,-1]
    rownames(DIMTWO) = make.unique(DIMTWO$VOC)
    DIMTWO = DIMTWO[,-1]
    rownames(MatchOne) = make.unique(MatchOne$VOC)
    MatchOne = MatchOne[,-1]
    rownames(MatchTwo) = make.unique(MatchTwo$VOC)
    MatchTwo = MatchTwo[,-1]
    
    
    logindex = matrix(NA,nrow = nr,ncol=nc) #logical index container
    for(i in 1:nr){
      for(j in 1:nc){ 
        
        #Set first and second dimension tolerance values for matching
        maxtol1 = max(DIMONE[i,],na.rm = T)+Dim1Tol
        mintol1 = min(DIMONE[i,],na.rm = T)-Dim1Tol
        maxtol2 = max(DIMTWO[i,],na.rm = T)+Dim2Tol
        mintol2 = min(DIMTWO[i,],na.rm = T)-Dim2Tol
        
        tmp1 = MatchOne[j,]
        tmp2 = MatchTwo[j,]
        
        #Check that VOC first and second dim retention times fall within the user defined tolerance
        if(min(tmp1,na.rm = T) > mintol1 && max(tmp1,na.rm = T) < maxtol1 && min(tmp2,na.rm = T) > mintol2 && max(tmp2,na.rm = T) < maxtol2){
          
          logindex[i,j] = 1
        }else{
          logindex[i,j] = 0
        }
      }
    }
    
    
    matchmatrixsize = max(c(sumrow(logindex),sumcol(logindex)))
    
    matchindex2 = matrix(NA,nrow=nrow(logindex),ncol=matchmatrixsize) #number of rows corresponds to the number of VOCs in the reference condition
    rownames(matchindex2) = rownames(DIMONE)
    colnames(matchindex2) = paste("Condition3Match",seq(1,matchmatrixsize,1))
    for(i in 1:nrow(logindex)){
      ti = length(which(logindex[i,] == 1))
      if(ti > 0){
        matchindex2[i,1:ti] = which(logindex[i,]==1)
      }else{
        matchindex2[i,1:ti+1] = NA
      }
    }
    
    PartialMatch = cbind(matchindex,matchindex2)
    
    #Define size of the logical index container
    n1 = paste("HighQuality",ref,"Abundance",sep="") #Dependent on reference level specified by the user
    c1 = get(n1)
    nr = nrow(c1)
    n2 = paste("Clean",Condition4Name,"Abundance",sep="") #The code will automatically select the second listed condition if more than 2 are available
    c2 = get(n2)
    nc = nrow(c2)
    
    #Get reference condition for alignment
    refnam = paste("HighQuality",ref,"DimOne",sep="")
    DIMONE = get(refnam)
    refnam = paste("HighQuality",ref,"DimTwo",sep="")
    DIMTWO = get(refnam)
    
    #Condition 3 retention times
    matchname = paste("Clean",Condition4Name,"DimOne",sep="")
    MatchOne = get(matchname)
    matchname = paste("Clean",Condition4Name,"DimTwo",sep="")
    MatchTwo = get(matchname)
    
    #Necessary cleanup
    rownames(DIMONE) = make.unique(DIMONE$VOC)
    DIMONE = DIMONE[,-1]
    rownames(DIMTWO) = make.unique(DIMTWO$VOC)
    DIMTWO = DIMTWO[,-1]
    rownames(MatchOne) = make.unique(MatchOne$VOC)
    MatchOne = MatchOne[,-1]
    rownames(MatchTwo) = make.unique(MatchTwo$VOC)
    MatchTwo = MatchTwo[,-1]
    
    
    logindex = matrix(NA,nrow = nr,ncol=nc) #logical index container
    for(i in 1:nr){
      for(j in 1:nc){ 
        
        #Set first and second dimension tolerance values for matching
        maxtol1 = max(DIMONE[i,],na.rm = T)+Dim1Tol
        mintol1 = min(DIMONE[i,],na.rm = T)-Dim1Tol
        maxtol2 = max(DIMTWO[i,],na.rm = T)+Dim2Tol
        mintol2 = min(DIMTWO[i,],na.rm = T)-Dim2Tol
        
        tmp1 = MatchOne[j,]
        tmp2 = MatchTwo[j,]
        
        #Check that VOC first and second dim retention times fall within the user defined tolerance
        if(min(tmp1,na.rm = T) > mintol1 && max(tmp1,na.rm = T) < maxtol1 && min(tmp2,na.rm = T) > mintol2 && max(tmp2,na.rm = T) < maxtol2){
          
          logindex[i,j] = 1
        }else{
          logindex[i,j] = 0
        }
      }
    }
    
    
    matchmatrixsize = max(c(sumrow(logindex),sumcol(logindex)))
    
    matchindex3 = matrix(NA,nrow=nrow(logindex),ncol=matchmatrixsize) #number of rows corresponds to the number of VOCs in the reference condition
    rownames(matchindex3) = rownames(DIMONE)
    colnames(matchindex3) = paste("Condition4Match",seq(1,matchmatrixsize,1))
    for(i in 1:nrow(logindex)){
      ti = length(which(logindex[i,] == 1))
      if(ti > 0){
        matchindex3[i,1:ti] = which(logindex[i,]==1)
      }else{
        matchindex3[i,1:ti+1] = NA
      }
    }
    
    FullMatchMatrix = cbind(PartialMatch,matchindex3)
    View(FullMatchMatrix)
    
  }else{
    cat("Sorry, this code does not currently support this analysis")
  }
}else if(numberofconditions == 5){
  
  cat("Aligning conditions 3-5...")
  
  #Define size of the logical index container
  n1 = paste("HighQuality",ref,"Abundance",sep="") #Dependent on reference level specified by the user
  c1 = get(n1)
  nr = nrow(c1)
  n2 = paste("Clean",Condition3Name,"Abundance",sep="") #The code will automatically select the second listed condition if more than 2 are available
  c2 = get(n2)
  nc = nrow(c2)
  
  #Get reference condition for alignment
  refnam = paste("HighQuality",ref,"DimOne",sep="")
  DIMONE = get(refnam)
  refnam = paste("HighQuality",ref,"DimTwo",sep="")
  DIMTWO = get(refnam)
  
  #Condition 3 retention times
  matchname = paste("Clean",Condition3Name,"DimOne",sep="")
  MatchOne = get(matchname)
  matchname = paste("Clean",Condition3Name,"DimTwo",sep="")
  MatchTwo = get(matchname)
  
  #Necessary cleanup
  rownames(DIMONE) = make.unique(DIMONE$VOC)
  DIMONE = DIMONE[,-1]
  rownames(DIMTWO) = make.unique(DIMTWO$VOC)
  DIMTWO = DIMTWO[,-1]
  rownames(MatchOne) = make.unique(MatchOne$VOC)
  MatchOne = MatchOne[,-1]
  rownames(MatchTwo) = make.unique(MatchTwo$VOC)
  MatchTwo = MatchTwo[,-1]
  
  
  logindex = matrix(NA,nrow = nr,ncol=nc) #logical index container
  for(i in 1:nr){
    for(j in 1:nc){ 
      
      #Set first and second dimension tolerance values for matching
      maxtol1 = max(DIMONE[i,],na.rm = T)+Dim1Tol
      mintol1 = min(DIMONE[i,],na.rm = T)-Dim1Tol
      maxtol2 = max(DIMTWO[i,],na.rm = T)+Dim2Tol
      mintol2 = min(DIMTWO[i,],na.rm = T)-Dim2Tol
      
      tmp1 = MatchOne[j,]
      tmp2 = MatchTwo[j,]
      
      #Check that VOC first and second dim retention times fall within the user defined tolerance
      if(min(tmp1,na.rm = T) > mintol1 && max(tmp1,na.rm = T) < maxtol1 && min(tmp2,na.rm = T) > mintol2 && max(tmp2,na.rm = T) < maxtol2){
        
        logindex[i,j] = 1
      }else{
        logindex[i,j] = 0
      }
    }
  }
  
  
  matchmatrixsize = max(c(sumrow(logindex),sumcol(logindex)))
  
  matchindex2 = matrix(NA,nrow=nrow(logindex),ncol=matchmatrixsize) #number of rows corresponds to the number of VOCs in the reference condition
  rownames(matchindex2) = rownames(DIMONE)
  colnames(matchindex2) = paste("Condition3Match",seq(1,matchmatrixsize,1))
  for(i in 1:nrow(logindex)){
    ti = length(which(logindex[i,] == 1))
    if(ti > 0){
      matchindex2[i,1:ti] = which(logindex[i,]==1)
    }else{
      matchindex2[i,1:ti+1] = NA
    }
  }
  
  PartialMatch = cbind(matchindex,matchindex2)
  
  #Define size of the logical index container
  n1 = paste("HighQuality",ref,"Abundance",sep="") #Dependent on reference level specified by the user
  c1 = get(n1)
  nr = nrow(c1)
  n2 = paste("Clean",Condition4Name,"Abundance",sep="") #The code will automatically select the second listed condition if more than 2 are available
  c2 = get(n2)
  nc = nrow(c2)
  
  #Get reference condition for alignment
  refnam = paste("HighQuality",ref,"DimOne",sep="")
  DIMONE = get(refnam)
  refnam = paste("HighQuality",ref,"DimTwo",sep="")
  DIMTWO = get(refnam)
  
  #Condition 3 retention times
  matchname = paste("Clean",Condition4Name,"DimOne",sep="")
  MatchOne = get(matchname)
  matchname = paste("Clean",Condition4Name,"DimTwo",sep="")
  MatchTwo = get(matchname)
  
  #Necessary cleanup
  rownames(DIMONE) = make.unique(DIMONE$VOC)
  DIMONE = DIMONE[,-1]
  rownames(DIMTWO) = make.unique(DIMTWO$VOC)
  DIMTWO = DIMTWO[,-1]
  rownames(MatchOne) = make.unique(MatchOne$VOC)
  MatchOne = MatchOne[,-1]
  rownames(MatchTwo) = make.unique(MatchTwo$VOC)
  MatchTwo = MatchTwo[,-1]
  
  
  logindex = matrix(NA,nrow = nr,ncol=nc) #logical index container
  for(i in 1:nr){
    for(j in 1:nc){ 
      
      #Set first and second dimension tolerance values for matching
      maxtol1 = max(DIMONE[i,],na.rm = T)+Dim1Tol
      mintol1 = min(DIMONE[i,],na.rm = T)-Dim1Tol
      maxtol2 = max(DIMTWO[i,],na.rm = T)+Dim2Tol
      mintol2 = min(DIMTWO[i,],na.rm = T)-Dim2Tol
      
      tmp1 = MatchOne[j,]
      tmp2 = MatchTwo[j,]
      
      #Check that VOC first and second dim retention times fall within the user defined tolerance
      if(min(tmp1,na.rm = T) > mintol1 && max(tmp1,na.rm = T) < maxtol1 && min(tmp2,na.rm = T) > mintol2 && max(tmp2,na.rm = T) < maxtol2){
        
        logindex[i,j] = 1
      }else{
        logindex[i,j] = 0
      }
    }
  }
  
  
  matchmatrixsize = max(c(sumrow(logindex),sumcol(logindex)))
  
  matchindex3 = matrix(NA,nrow=nrow(logindex),ncol=matchmatrixsize) #number of rows corresponds to the number of VOCs in the reference condition
  rownames(matchindex3) = rownames(DIMONE)
  colnames(matchindex3) = paste("Condition4Match",seq(1,matchmatrixsize,1))
  for(i in 1:nrow(logindex)){
    ti = length(which(logindex[i,] == 1))
    if(ti > 0){
      matchindex3[i,1:ti] = which(logindex[i,]==1)
    }else{
      matchindex3[i,1:ti+1] = NA
    }
  }
  
  PartialMatch2 = cbind(PartialMatch,matchindex3)
  
  #Define size of the logical index container
  n1 = paste("HighQuality",ref,"Abundance",sep="") #Dependent on reference level specified by the user
  c1 = get(n1)
  nr = nrow(c1)
  n2 = paste("Clean",Condition4Name,"Abundance",sep="") #The code will automatically select the second listed condition if more than 2 are available
  c2 = get(n2)
  nc = nrow(c2)
  
  #Get reference condition for alignment
  refnam = paste("HighQuality",ref,"DimOne",sep="")
  DIMONE = get(refnam)
  refnam = paste("HighQuality",ref,"DimTwo",sep="")
  DIMTWO = get(refnam)
  
  #Condition 3 retention times
  matchname = paste("Clean",Condition5Name,"DimOne",sep="")
  MatchOne = get(matchname)
  matchname = paste("Clean",Condition5Name,"DimTwo",sep="")
  MatchTwo = get(matchname)
  
  #Necessary cleanup
  rownames(DIMONE) = make.unique(DIMONE$VOC)
  DIMONE = DIMONE[,-1]
  rownames(DIMTWO) = make.unique(DIMTWO$VOC)
  DIMTWO = DIMTWO[,-1]
  rownames(MatchOne) = make.unique(MatchOne$VOC)
  MatchOne = MatchOne[,-1]
  rownames(MatchTwo) = make.unique(MatchTwo$VOC)
  MatchTwo = MatchTwo[,-1]
  
  
  logindex = matrix(NA,nrow = nr,ncol=nc) #logical index container
  for(i in 1:nr){
    for(j in 1:nc){ 
      
      #Set first and second dimension tolerance values for matching
      maxtol1 = max(DIMONE[i,],na.rm = T)+Dim1Tol
      mintol1 = min(DIMONE[i,],na.rm = T)-Dim1Tol
      maxtol2 = max(DIMTWO[i,],na.rm = T)+Dim2Tol
      mintol2 = min(DIMTWO[i,],na.rm = T)-Dim2Tol
      
      tmp1 = MatchOne[j,]
      tmp2 = MatchTwo[j,]
      
      #Check that VOC first and second dim retention times fall within the user defined tolerance
      if(min(tmp1,na.rm = T) > mintol1 && max(tmp1,na.rm = T) < maxtol1 && min(tmp2,na.rm = T) > mintol2 && max(tmp2,na.rm = T) < maxtol2){
        
        logindex[i,j] = 1
      }else{
        logindex[i,j] = 0
      }
    }
  }
  
  
  matchmatrixsize = max(c(sumrow(logindex),sumcol(logindex)))
  
  matchindex4 = matrix(NA,nrow=nrow(logindex),ncol=matchmatrixsize) #number of rows corresponds to the number of VOCs in the reference condition
  rownames(matchindex4) = rownames(DIMONE)
  colnames(matchindex4) = paste("Condition5Match",seq(1,matchmatrixsize,1))
  for(i in 1:nrow(logindex)){
    ti = length(which(logindex[i,] == 1))
    if(ti > 0){
      matchindex4[i,1:ti] = which(logindex[i,]==1)
    }else{
      matchindex4[i,1:ti+1] = NA
    }
  }
  
  FullMatchMatrix = cbind(PartialMatch2,matchindex4)
  View(FullMatchMatrix)
}else{
  cat("This code does not support more than 5 experimental conditions")
}

###############################################################################################################################
# MS Comparison (performed in ChromaTOF)
###############################################################################################################################


a = paste("HighQuality",ref,"DimOne",sep="")
a = get(a)
View(a)
b = paste("HighQuality",ref,"DimTwo",sep="")
b = get(b)
View(b)
if(length(expconditions) == 2){
  c = paste("HighQuality",expconditions[! expconditions %in% ref],"DimOne",sep="")
  c = get(c)
  c = View(c)
  d = paste("HighQuality",expconditions[! expconditions %in% ref],"DimTwo",sep="")
  d = get(d)
  d = View(d)
}else if(length(expconditions) > 2){
  c = paste("HighQuality",Condition2Name,"DimOne",sep="")
  c = get(c)
  c = View(c)
  d = paste("HighQuality",Condition2Name,"DimTwo",sep="")
  d = get(d)
  d = View(d)
}

if(writefile == T){
  setwd(writewd)
  write.csv(FullMatchMatrix,"IsotopicallyLabeled_AlignmentMatrix.csv")
  write.csv(HighQualityLightDimOne,"HighQuality_Light_Dim1.csv")
  write.csv(HighQualityLightDimTwo,"HighQuality_Light_Dim2.csv")
  write.csv(HighQualityHeavyDimOne,"HighQuality_Heavy_Dim1.csv")
  write.csv(HighQualityHeavyDimTwo,"HighQuality_Heavy_Dim2.csv")
}
