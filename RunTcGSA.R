library(TcGSA)
library(GSA)

# Construct gene sets (file loaded from Broad institutes web page)
setwd("/path/to/genesets")
load("NoiseLevelData.RData")
sampleinfo = read.table("SampleInfo_overlap10_200donors.txt", sep="\t", quote="")

# Load in data
setwd("path/to/data")
load("DataAndInfo.RData")

# DAISY
cases = grep("TD", colnames(data_DAISY), value=T)
tcgsa_DAISY = TcGSA.LR(expr=data_DAISY[,cases], gmt=genesets, design=info_DAISY[cases,c("Donor","Age","TimeFromSeroconversion")], 
                      subject_name="Donor", time_name="TimeFromSeroconversion", covariates_fixed="Age", minGSsize=1)
significant_tcgsa_DAISY = signifLRT.TcGSA(tcgsa_DAISY)$mixedLRTadjRes 
res_tcgsa_DAISY = multtest.TcGSA(tcgsa_DAISY)

# Diabimmune_PBMC
cases = grep("Case", colnames(data_pbmc), value=T)
tcgsa_pbmc = TcGSA.LR(expr=data_pbmc[,cases], gmt=genesets, design=info_Diabimmune[cases,c("Donor","Age","TimeFromSero")], 
                      subject_name="Donor", time_name="TimeFromSero", covariates_fixed="Age", minGSsize=5)
significant_tcgsa_pbmc = signifLRT.TcGSA(tcgsa_pbmc)$mixedLRTadjRes 
res_tcgsa_pbmc = multtest.TcGSA(tcgsa_pbmc)

# Diabimmune_CD4
cases = grep("Case", colnames(data_cd4), value=T)
tcgsa_cd4 = TcGSA.LR(expr=data_cd4[,cases], gmt=genesets, design=info_Diabimmune[cases,c("Donor","Age","TimeFromSero")], 
                      subject_name="Donor", time_name="TimeFromSero", covariates_fixed="Age", minGSsize=5)
significant_tcgsa_cd4 = signifLRT.TcGSA(tcgsa_cd4)$mixedLRTadjRes 
res_tcgsa_cd4 = multtest.TcGSA(tcgsa_cd4)

# Diabimmune_CD8
cases = grep("Case", colnames(data_cd8), value=T)
tcgsa_cd8 = TcGSA.LR(expr=data_cd8[,cases], gmt=genesets, design=info_Diabimmune[cases,c("Donor","Age","TimeFromSero")], 
                     subject_name="Donor", time_name="TimeFromSero", covariates_fixed="Age", minGSsize=5)
significant_tcgsa_cd8 = signifLRT.TcGSA(tcgsa_cd8)$mixedLRTadjRes 
res_tcgsa_cd8 = multtest.TcGSA(tcgsa_cd8)

# BabyDiet
cases = rownames(info_BabyDiet_MultiDonor)[info_BabyDiet_MultiDonor$SampleGroup == "Case"]
tcgsa_BabyDiet_MultiDonor = TcGSA.LR(expr=data_BabyDiet_MultiDonor[,cases], gmt=genesets, design=info_BabyDiet_MultiDonor[cases,c("Donor","Age","TimeFromSero")], 
                          subject_name="Donor", time_name="TimeFromSero", covariates_fixed="Age", minGSsize=5)
significant_tcgsa_BabyDiet_MultiDonor = signifLRT.TcGSA(tcgsa_BabyDiet_MultiDonor)$mixedLRTadjRes 
res_tcgsa_BabyDiet_MultiDonor = multtest.TcGSA(tcgsa_BabyDiet_MultiDonor)

# Save results
setwd("C:/Users/makija/Documents/Research/Project 14 PASI for longitudinal data/Results/TcGSA results")
save(file="Results_TcGSA.RData", list=c("tcgsa_DAISY", "res_tcgsa_DAISY", "tcgsa_pbmc", "res_tcgsa_pbmc", "tcgsa_cd4", "res_tcgsa_cd4",
                                        "tcgsa_cd8", "res_tcgsa_cd8", "tcgsa_BabyDiet", "res_tcgsa_BabyDiet", "tcgsa_BabyDiet_MultiDonor", "res_tcgsa_BabyDiet_MultiDonor"))
write.table(significant_tcgsa_DAISY, file="TcGSA_significant_DAISY.txt", sep="\t", quote=F, row.names=F)
write.table(significant_tcgsa_pbmc, file="TcGSA_significant_PBMC.txt", sep="\t", quote=F, row.names=F)
write.table(significant_tcgsa_cd4, file="TcGSA_significant_CD4.txt", sep="\t", quote=F, row.names=F)
write.table(significant_tcgsa_cd8, file="TcGSA_significant_CD8.txt", sep="\t", quote=F, row.names=F)
write.table(significant_tcgsa_BabyDiet, file="TcGSA_significant_BabyDiet.txt", sep="\t", quote=F, row.names=F)
write.table(significant_tcgsa_BabyDiet_MultiDonor, file="TcGSA_significant_BabyDiet_MultiDonor.txt", sep="\t", quote=F, row.names=F)

##### Simulated data

# Read in gene sets
setwd("/path/to/simulated/genesets")
genesets = GSA.read.gmt("GeneSets.gmt")

# Read in data
setwd("/path/to/simulated/data")
load("NoiseLevelData.RData")
sampleinfo = read.table("SampleInfo_overlap10_200donors.txt", sep="\t", quote="")
data = list(data1, data5, data10, data15, data20, data30, data40)
names(data) = c("data1", "data5", "data10", "data15", "data20", "data30", "data40")

# Run TcGSA
setwd("/path/to/result/dir")
res_TcGSA = list()
for(i in 1:length(data)){
  
  # Run TcGSA
  tempdata = data[[i]]
  cases = grep("Case", colnames(tempdata), value=T)
  tcgsa_sim = TcGSA.LR(expr=tempdata[,cases], gmt=genesets, design=sampleinfo[cases,c("Donor","Age","TimeFromEvent")], 
                       subject_name="Donor", time_name="TimeFromEvent", covariates_fixed="Age", minGSsize=1)
  significant_tcgsa_sim = signifLRT.TcGSA(tcgsa_sim)$mixedLRTadjRes 
  res_tcgsa_sim = multtest.TcGSA(tcgsa_sim)
  
  # Save results
  filename = paste("ResultsSimulated_TcGSA_", ".RData", sep=names(data)[i])
  save(file=filename, list=c("tcgsa_sim", "res_tcgsa_sim"))
  filename2 = paste("TcGSA_significant_Simulated_", ".txt", sep=names(data)[i])
  write.table(significant_tcgsa_sim, file=filename2, sep="\t", quote=F, row.names=F)
  
  res_TcGSA[[i]] = tcgsa_sim
}


