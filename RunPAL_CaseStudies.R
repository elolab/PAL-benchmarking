library(robustlmm)
library(lme4)
library(lmerTest)
library(MASS)
library(PASI)

# Define where to save results from all runs
respath = "path/to/result/dir"

# Read in info tables
setwd("path/to/sample/info")
load("CollectedInfo.RData")

######### DAISY

# Read in data
setwd("path/to/daisy/data")
data = read.table("PXD007884_Entrez.txt", sep="\t", quote="", header=1)
info_DAISY = info_DAISY[colnames(data),]

# Keep only genes that are expressed in some data point by at least 2 case and 2 control donors
dropind = which(apply(data, 1, function(r){
  expresses_case = which(!is.na(r[info_DAISY$Labels != 0]))
  expresses_control = which(!is.na(r[info_DAISY$Labels == 0]))
  if((length(expresses_case)==0) | (length(expresses_control)==0)) return(T)
  expressingdonors_case = unique(gsub("_.*","",colnames(data)[info_DAISY$Labels != 0][expresses_case]))
  expressingdonors_control = unique(gsub("_.*","",colnames(data)[info_DAISY$Labels == 0][expresses_control]))
  if((length(expressingdonors_case) < 2) | (length(expressingdonors_control) < 2)) return(T)
  return(F)
}))
if(length(dropind) > 0) data = data[-dropind,]

setwd(respath)

# Run PAL with neutralizing for age and using time from seroconversion as the feature of interest
pal_time_DAISY = PAL(data=data, info=info_DAISY, grouplabels="Labels", neutralise="Age", mainfeature="TimeFromSeroconversion", userandom="Donor", 
                     neutralisationmodel="lmer", pathwaymodel="lmer", score="activity", nodemin=1, seed=1234)
write.table(pal_time_DAISY$pathwayscores, file="PAL_DAISY_NoAge_TimeFromSero.txt", sep="\t", quote=F)
write.table(pal_time_DAISY$significance, file="PAL_DAISY_NoAge_TimeFromSero_fdr.txt", sep="\t", quote=F)

# Run PAL with neutralizing for age and using sample groups as the feature of interest
pal_group_DAISY = PAL(data=data, info=info_DAISY, grouplabels="Labels", neutralise="Age", mainfeature="SampleGroup", userandom="Donor", 
                      neutralisationmodel="lmer", pathwaymodel="lmer", score="activity", nodemin=1, seed=1234)
write.table(pal_group_DAISY$pathwayscores, file="PAL_DAISY_NoAge_SampleGroup.txt", sep="\t", quote=F)
write.table(pal_group_DAISY$significance, file="PAL_DAISY_NoAge_SampleGroup_fdr.txt", sep="\t", quote=F)

# Write significant detections into a file
daisyfindings_time = rownames(pal_time_DAISY$significance)[which(as.numeric(pal_time_DAISY$significance[,"FDR_TimeFromSeroconversion"]) < 0.05)]
daisyfindings_group = rownames(pal_group_DAISY$significance)[which(as.numeric(pal_group_DAISY$significance[,"FDR_SampleGroupcontrol"]) < 0.05)]
linestowrite = c("TimeFromSero_DAISY:", daisyfindings_time, "", "SampleGroup_DAISY:", daisyfindings_group)
write(linestowrite, file="PAL_SignificantDetections_DAISY.txt", sep="\n")

######### Diabimmune

# Read in data
setwd("")
data_pbmc = read.table("Diabimmune_Entrez_PBMC.txt", sep="\t", quote="", header=1, check.names=F)
data_cd4 = read.table("Diabimmune_Entrez_CD4.txt", sep="\t", quote="", header=1, check.names=F)
data_cd8 = read.table("Diabimmune_Entrez_CD8.txt", sep="\t", quote="", header=1, check.names=F)

# Keep only genes with median expression above 1 in either sample group
caseind_pbmc = which(info_Diabimmune[colnames(data_pbmc),"SampleGroup"] == "Case")
controlind_pbmc = which(info_Diabimmune[colnames(data_pbmc),"SampleGroup"] == "Control")
data_pbmc = data_pbmc[pmax(apply(data_pbmc[,caseind_pbmc],1,median),apply(data_pbmc[,controlind_pbmc]==0,1,median)) > 1,]
caseind_cd4 = which(info_Diabimmune[colnames(data_cd4),"SampleGroup"] == "Case")
controlind_cd4 = which(info_Diabimmune[colnames(data_cd4),"SampleGroup"] == "Control")
data_cd4 = data_cd4[pmax(apply(data_cd4[,caseind_cd4],1,median),apply(data_cd4[,controlind_cd4]==0,1,median)) > 1,]
caseind_cd8 = which(info_Diabimmune[colnames(data_cd8),"SampleGroup"] == "Case")
controlind_cd8 = which(info_Diabimmune[colnames(data_cd8),"SampleGroup"] == "Control")
data_cd8 = data_cd8[pmax(apply(data_cd8[,caseind_cd8],1,median),apply(data_cd8[,controlind_cd8]==0,1,median)) > 1,]

setwd(respath)

# Run PAL with neutralizing for age and using time from seroconversion as the feature of interest
pal_time_pbmc = PAL(data=log2(data_pbmc+1), info=info_Diabimmune[colnames(data_pbmc),,drop=F], grouplabels="Label", neutralise="Age", mainfeature="TimeFromSero", 
                    userandom="Donor", neutralisationmodel="lmer", pathwaymodel="lmer", score="activity", seed=1234)
write.table(pal_time_pbmc$pathwayscores, file="PAL_Diabimmune_PBMC_NoAge_TimeFromSero.txt", sep="\t", quote=F)
write.table(pal_time_pbmc$significance, file="PAL_Diabimmune_PBMC_NoAge_TimeFromSero_fdr.txt", sep="\t", quote=F)

pal_time_cd4 = PAL(log2(data_cd4+1), info=info_Diabimmune[colnames(data_cd4),,drop=F], grouplabels="Label", neutralise="Age", mainfeature="TimeFromSero", 
                   userandom="Donor", neutralisationmodel="lmer", pathwaymodel="lmer", score="activity", seed=1234)
write.table(pal_time_cd4$pathwayscores, file="PAL_Diabimmune_CD4_NoAge_TimeFromSero.txt", sep="\t", quote=F)
write.table(pal_time_cd4$significance, file="PAL_Diabimmune_CD4_NoAge_TimeFromSero_fdr.txt", sep="\t", quote=F)

pal_time_cd8 = PAL(log2(data_cd8+1), info=info_Diabimmune[colnames(data_cd8),,drop=F], grouplabels="Label", neutralise="Age", mainfeature="TimeFromSero", 
                   userandom="Donor", neutralisationmodel="lmer", pathwaymodel="lmer", score="activity", seed=1234)
write.table(pal_time_cd8$pathwayscores, file="PAL_Diabimmune_CD8_NoAge_TimeFromSero.txt", sep="\t", quote=F)
write.table(pal_time_cd8$significance, file="PAL_Diabimmune_CD8_NoAge_TimeFromSero_fdr.txt", sep="\t", quote=F)

# Run PAL with neutralizing for age and using samplegroup as the feature of interest
pal_group_pbmc = PAL(log2(data_pbmc+1), info=info_Diabimmune[colnames(data_pbmc),,drop=F], grouplabels="Label", neutralise="Age", mainfeature="SampleGroup", 
                     userandom="Donor", neutralisationmodel="lmer", pathwaymodel="lmer", score="activity", seed=1234)
write.table(pal_group_pbmc$pathwayscores, file="PAL_Diabimmune_PBMC_NoAge_SampleGroup.txt", sep="\t", quote=F)
write.table(pal_group_pbmc$significance, file="PAL_Diabimmune_PBMC_NoAge_SampleGroup_fdr.txt", sep="\t", quote=F)

pal_group_cd4 = PAL(log2(data_cd4+1), info=info_Diabimmune[colnames(data_cd4),,drop=F], grouplabels="Label", neutralise="Age", mainfeature="SampleGroup", 
                    userandom="Donor", neutralisationmodel="lmer", pathwaymodel="lmer", score="activity", seed=1234)
write.table(pal_group_cd4$pathwayscores, file="PAL_Diabimmune_CD4_NoAge_SampleGroup.txt", sep="\t", quote=F)
write.table(pal_group_cd4$significance, file="PAL_Diabimmune_CD4_NoAge_SampleGroup_fdr.txt", sep="\t", quote=F)

pal_group_cd8 = PAL(log2(data_cd8+1), info=info_Diabimmune[colnames(data_cd8),,drop=F], grouplabels="Label", neutralise="Age", mainfeature="SampleGroup", 
                    userandom="Donor", neutralisationmodel="lmer", pathwaymodel="lmer", score="activity", seed=1234)
write.table(pal_group_cd8$pathwayscores, file="PAL_Diabimmune_CD8_NoAge_SampleGroup.txt", sep="\t", quote=F)
write.table(pal_group_cd8$significance, file="PAL_Diabimmune_CD8_NoAge_SampleGroup_fdr.txt", sep="\t", quote=F)

# Write significant detections into a file
pbmcfindings_time = rownames(pal_time_pbmc$significance)[which(as.numeric(pal_time_pbmc$significance[,"FDR_TimeFromSero"]) < 0.05)]
cd4findings_time = rownames(pal_time_cd4$significance)[which(as.numeric(pal_time_cd4$significance[,"FDR_TimeFromSero"]) < 0.05)]
cd8findings_time = rownames(pal_time_cd8$significance)[which(as.numeric(pal_time_cd8$significance[,"FDR_TimeFromSero"]) < 0.05)]
pbmcfindings_group = rownames(pal_group_pbmc$significance)[which(as.numeric(pal_group_pbmc$significance[,"FDR_SampleGroupControl"]) < 0.05)]
cd4findings_group = rownames(pal_group_cd4$significance)[which(as.numeric(pal_group_cd4$significance[,"FDR_SampleGroupControl"]) < 0.05)]
cd8findings_group = rownames(pal_group_cd8$significance)[which(as.numeric(pal_group_cd8$significance[,"FDR_SampleGroupControl"]) < 0.05)]
linestowrite = c("TimeFromSero_PBMC:",pbmcfindings_time,"","TimeFromSero_CD4:",cd4findings_time,"","TimeFromSero_CD8:", cd8findings_time, "",
                 "SampleGroup_PBMC:",pbmcfindings_group,"","SampleGroup_CD4:",cd4findings_group,"","SampleGroup_CD8:", cd8findings_group, "")
write(linestowrite, file="PAL_SignificantDetections_Diabimmune.txt", sep="\n")

######### BabyDiet

# Read in data
setwd("")
data = read.table("MTAB1724_BabyDiet_Entrez.txt", sep="\t", quote="", header=1)
keepsamples = intersect(colnames(data), rownames(info_BabyDiet))
info_BabyDiet = info_BabyDiet[keepsamples,] # Note: info_BabyDiet does not include age outliers, so they are excluded here
data = data[,keepsamples]

# Remove donors that gave only one sample
donorfreqs = table(info_BabyDiet$Donor)
dropdonors = names(donorfreqs)[which(donorfreqs == 1)]
dropindex = which(info_BabyDiet$Donor %in% dropdonors)
info_BabyDiet = info_BabyDiet[-dropindex,]
data = data[,-dropindex]

setwd(respath)

# Run PAL with neutralizing for age and using time from serovonversion as the feature of interest
pal_time_BabyDiet = PAL(data=data, info=info_BabyDiet, grouplabels="Labels", neutralise="Age", mainfeature="TimeFromSero", userandom="Donor", 
                        neutralisationmodel="lmer", pathwaymodel="lmer", score="activity", seed=1234)
write.table(pal_time_BabyDiet$pathwayscores, file="PAL_BabyDiet_NoAge_TimeFromSero.txt", sep="\t", quote=F)
write.table(pal_time_BabyDiet$significance, file="PAL_BabyDiet_NoAge_TimeFromSero_fdr.txt", sep="\t", quote=F)

# Run PAL with neutralizing for age and using sample group as the feature of interest
pal_group_BabyDiet = PAL(data=data, info=info_BabyDiet, grouplabels="Labels", neutralise="Age", mainfeature="SampleGroup", userandom="Donor", 
                         neutralisationmodel="lmer", pathwaymodel="lmer", score="activity", seed=1234)
write.table(pal_group_BabyDiet$pathwayscores, file="PAL_BabyDiet_NoAge_SampleGroup.txt", sep="\t", quote=F)
write.table(pal_group_BabyDiet$significance, file="PAL_BabyDiet_NoAge_SampleGroup_fdr.txt", sep="\t", quote=F)

# Write significant detections into a file
babydietfindings_time = rownames(pal_time_BabyDiet$significance)[which(as.numeric(pal_time_BabyDiet$significance[,"FDR_TimeFromSero"]) < 0.05)]
babydietfindings_group = rownames(pal_group_BabyDiet$significance)[which(as.numeric(pal_group_BabyDiet$significance[,"FDR_SampleGroupControl"]) < 0.05)]
linestowrite = c("TimeFromSero_BabyDiet:", babydietfindings_time, "", "SampleGroup_BabyDiet:", babydietfindings_group)
write(linestowrite, file="PAL_SignificantDetections_BabyDiet.txt", sep="\n")

