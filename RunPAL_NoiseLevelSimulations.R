library(lme4)
library(robustlmm)
library(PAL)

# Read in data
setwd("/path/to/data")
load("NoiseLevelData.RData")
sampleinfo = read.table("SampleInfo_overlap10_200donors.txt", sep="\t", quote="")

# Prepare input
sampleinfo[,"GroupLabel"] = rep(0, nrow(sampleinfo))
sampleinfo[which(sampleinfo$SampleGroup == "Case"),"GroupLabel"] = 1
path = "/path/to/pathways"
datas = list(data1=as.data.frame(data1), data5=as.data.frame(data5), 
             data10=as.data.frame(data10), data15=as.data.frame(data15),
             data20=as.data.frame(data20), data30=as.data.frame(data30),
             data40=as.data.frame(data40))

# Run PAL for different noise levels
reslist = list()
for(i in 1:length(datas)){
  tempdata = datas[[i]]
  genesd = apply(tempdata[,grep("Control",colnames(tempdata),value=T)], 1, sd) 
  tempdata = tempdata[genesd>0,]
  filename = paste("ResultsPAL_Noise", ".RData", sep=gsub("data","",names(datas[i])))
  res = PAL(tempdata, sampleinfo, "GroupLabel", adjust="Age", mainfeature="TimeFromEvent", 
             userandom="Donor", adjustmentformula=NULL, pathwayformula=NULL, 
             adjustmentmodel="lmer", pathwaymodel="lmer",  pathwayadress=path, 
             useKEGG=FALSE, score="activity", nodemin=5, seed=1234)
  save(file=filename, list=c("res"))
  reslist[[i]] = res
}
names(reslist) = paste("Noise", c(1,5,10,15,20,30,40), sep="")
save(file="ResultsPAL_NoiseLevels.RData", list=c("reslist"))
  
