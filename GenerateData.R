library(graphsim) # version 1.0.3
library(igraph)   # version 1.3.4

# Generate pathways and save them in PAL format
pathways = GeneratePathways(100, minsize=5, maxsize=60, overlapprop=0.1, seed=252525)
setwd("/path/to/wanted/dir")
WritePathways(pathways)
WritePathwayGmt(pathways, "GeneSets.gmt")

# Generate base expression level and variation in it for all pathway genes
allgenes = unique(unlist(lapply(pathways, function(p){V(p)$name})))
set.seed(252525)
baselevel = rnorm(length(allgenes), mean=10, sd=4)
genesd = rnorm(length(allgenes), mean=0.1, sd=0.02)
baselevel[baselevel < 0] = 0
names(baselevel) = allgenes
genesd[genesd < 0] = 0.01
genesd = genesd * baselevel
names(genesd) = allgenes

casenumber = 100
controlnumber = 100
donornumber = casenumber + controlnumber

# Generate ages for samples from donors (no duplicated time points (two decimals) within one donor)
seed = 252525
continue = T
while(continue){
  set.seed(seed) 
  timepoints = sample(2:8, size=donornumber, replace=T)
  ages = lapply(timepoints, function(x){sort(runif(x, min=1, max=10))})
  samplenumber = sum(timepoints)
  samplenames = paste(rep(paste(rep(c("Case","Control"),times=c(casenumber,controlnumber)), c(1:casenumber, 1:controlnumber), sep=""),times=timepoints), round(unlist(ages),digits=2), sep="_")
  if(all(!duplicated(samplenames))){
    continue = F
  } else seed = seed + 1
}

# Generate event times for 20 samples and calculate time from it
set.seed(252525)
event = runif(casenumber, min=0, max=10)
timefromevent = mapply(function(age, eventtime){
  return(age - eventtime)
}, age=ages[1:casenumber], eventtime=event)

# Collect sample info
casesamples = grep("Case", samplenames, value=T)
sampleinfo = as.data.frame(matrix(NA, nrow=samplenumber, ncol=4))
rownames(sampleinfo) = samplenames
colnames(sampleinfo) = c("Age","TimeFromEvent","Donor","SampleGroup")
sampleinfo[,"Age"] = unlist(ages)
sampleinfo[casesamples,"TimeFromEvent"] = unlist(timefromevent)
sampleinfo[,"Donor"] = as.character(gsub("_.*", "", samplenames))
sampleinfo[,"SampleGroup"] = gsub("\\d*$", "", sampleinfo[,"Donor"])

# Randomly pick 20 pathways with disease progression effect and 60 with age effect (check that these partly overlap)
set.seed(252525)
hitpathways = names(pathways)[sample(1:length(pathways), size=20)]
agepathways = names(pathways)[sample(1:length(pathways), size=60)]
length(intersect(hitpathways, agepathways)) # overlap: 9

# Define which nodes in those pathways are without the effect (realistic situation: not the whole pathway is altered neatly)
nonessentialnodes = IdentifyUninterestingNodes(pathways)

# Generate and fill expression related to pathways without age or disease effect
basedata = GenerateBaseExpression(means=baselevel, deviations=genesd, sampleinfo, 
                                  pathways, allgenes, agepathways, hitpathways, 
                                  seed=252525)
baseexpression = basedata$data

# Generate expression changes related to pathways affected by age
ageeffect = GenerateAgeExpression(agepathways, mineffect=0.01, maxeffect=0.05, 
                                  baseexpression, nonessentialnodes, pathways, 
                                  sampleinfo, seed=252525)
ageexpression = ageeffect$ageexpression

# Generate expression changes related to pathways that vary based on disease progress
diseaseeffect = GenerateDiseaseExpression(hitpathways, mineffect=0.01, maxeffect=0.05, 
                                          baseexpression, nonessentialnodes, pathways, 
                                          sampleinfo, seed=12345)
diseaseexpression = diseaseeffect$diseaseexpression

# Add donor effect
donoreffect = GenerateDonorEffect(baselevel, rownames(sampleinfo), mean=0, sdprop=0.05, seed=252525)

# Add noise
noise1 = GenerateNoise(baselevel, rownames(sampleinfo), mean=0, sdprop=0.01, seed=12345)
noise5 = GenerateNoise(baselevel, rownames(sampleinfo), mean=0, sdprop=0.05, seed=12345)
noise10 = GenerateNoise(baselevel, rownames(sampleinfo), mean=0, sdprop=0.1, seed=12345)
noise15 = GenerateNoise(baselevel, rownames(sampleinfo), mean=0, sdprop=0.15, seed=12345)
noise20 = GenerateNoise(baselevel, rownames(sampleinfo), mean=0, sdprop=0.2, seed=12345)
noise30 = GenerateNoise(baselevel, rownames(sampleinfo), mean=0, sdprop=0.3, seed=12345)
noise40 = GenerateNoise(baselevel, rownames(sampleinfo), mean=0, sdprop=0.4, seed=12345)

# Construct data
data = baseexpression + ageexpression + diseaseexpression + donoreffect
data1 = data + noise1
data5 = data + noise5
data10 = data + noise10
data15 = data + noise15
data20 = data + noise20
data30 = data + noise30
data40 = data + noise40

# Truncate values (only data40 had any (7) sch values)
data1[data1 > 50] = 50
data5[data5 > 50] = 50
data10[data10 > 50] = 50
data15[data15 > 50] = 50
data20[data20 > 50] = 50
data30[data30 > 50] = 50
data40[data40 > 50] = 50

# Collect info of simulated pathways
pathwayinfo = matrix(0, length(pathways), 5)
rownames(pathwayinfo) = names(pathways)
colnames(pathwayinfo) = c("DiseaseEffect","AgeEffect","Nodes","Relations","MaxCor")
pathwayinfo[hitpathways,"DiseaseEffect"] = diseaseeffect$diseasecoef
pathwayinfo[agepathways,"AgeEffect"] = ageeffect$agecoef
pathwayinfo[,"Nodes"] = unlist(lapply(pathways, function(p){length(V(p))}))
pathwayinfo[,"Relations"] = as.numeric(unlist(lapply(pathways, function(p){length(E(p))})))
pathwayinfo[,"MaxCor"] = basedata$maxcors

# Save
setwd("/where/to/save/stuff")
save(file="Pathways.RData", list=c("pathways"))
save(file="NoiseLevelData.RData", list=c("data1","data5","data10","data15","data20","data30","data40"))
write.table(pathwayinfo, file="PathwayInfo_overlap10.txt", sep="\t", quote=F)
write.table(sampleinfo, file="SampleInfo_overlap10_200donors.txt", sep="\t", quote=F)





#
#
#
GeneratePathways = function(n, minsize=5, maxsize=60, overlapprop=0.2, seed=252525){
  
  # Randomly set  number of nodes and probability parameter for edge selection
  set.seed(seed)
  nodenumbers = sample(minsize:maxsize, size=n, replace=T)
  
  # Generate graphs using the randomly selected parameters above
  graphs = list()
  for(i in 1:n){
    continue = T
    j = 1
    while(continue){
      j = j + 1
      set.seed(j)
      currentgraph = sample_gnp(nodenumbers[i], 0.04, directed=T, loops=F) # set edge type?
      if(is_dag(currentgraph)) continue = F
    }
    graphs[[i]] = currentgraph
    print(i)
  }
  
  # Set overlapping genes appearing in multiple pathways
  set.seed(seed)
  genenames = sample(as.character(1:sum(nodenumbers)))
  if(overlapprop > 0){
    dupgenes = sample(genenames, size=round(overlapprop*sum(nodenumbers)))
    dupindex = which(genenames %in% dupgenes)
    genenames[dupindex] = sample(dupgenes, size=length(dupindex), replace=T)
  }
  
  # Set gene ids for the nodes so that there is realistic overlap between the graphs/pathway
  pathways = list()
  for(k in 1:n){
    nodenames = genenames[1:nodenumbers[k]]
    genenames = genenames[-c(1:nodenumbers[k])]
    pathways[[k]] = set.vertex.attribute(graphs[[k]], "name", value=nodenames)
    pathways[[k]] = set.vertex.attribute(pathways[[k]], "label", value=nodenames)
  }
  names(pathways) = paste("Pathway", 1:n, sep="")
  
  return(pathways)
}


#
#
#
GenerateNoise = function(baselevels, samplenames, mean=0, sdprop=0.05, seed=12345){
  
  #
  set.seed(seed)
  noise = t(as.data.frame(lapply(baselevels, function(b){rnorm(length(samplenames), mean=mean, sd=sdprop*b)})))
  
  return(as.matrix(noise))
}



#
#
#
GenerateBaseExpression = function(means, deviations, sampleinfo, pathways, allgenes, agepathways, diseasepathways, seed=252525){
  
  # Initialise expression matrix
  data = matrix(0, nrow=length(allgenes), ncol=nrow(sampleinfo))
  rownames(data) = allgenes
  colnames(data) = rownames(sampleinfo)
  
  # Fill in expression
  notrendpathways = setdiff(names(pathways), unique(c(diseasepathways, agepathways)))
  set.seed(seed)
  maxcors = rep(0.9, length(pathways))
  notrendindex = match(notrendpathways, names(pathways), 0)
  maxcors[notrendindex] = runif(length(notrendpathways), min=0.1, max=0.9)
  maxcors_remember = maxcors
  pathwayorder = names(pathways)[order(maxcors, decreasing=TRUE)]
  for(p in pathways){ 
    # p = pathways[[z]]
    pathwaygenes = as.character(V(p)$name)
    generatedexpression = generate_expression(samplenumber, p, cor=maxcors[1], mean=baselevel[pathwaygenes], sd=deviations[pathwaygenes])
    maxcors = maxcors[-1]
    data[rownames(generatedexpression),] = generatedexpression
  }
  
  return(list(data=data, maxcors=maxcors_remember))
}

#
#
#
GenerateDonorEffect = function(baselevels, samplenames, mean=0, sdprop=0.05, seed=252525){
  
  #
  donors = gsub("_.*", "", samplenames)
  donorsamples = rle(donors)$lengths
  donornumber = length(unique(donors))
  
  #
  set.seed(seed)
  donoreffect = t(as.data.frame(lapply(baselevels, function(b){rep(rnorm(donornumber, mean=mean, sd=sdprop*b), times=donorsamples)})))
  
  return(as.matrix(donoreffect))
}


# Summary: Writes pathways into separate files in the format used by PAL.
# Input:   'plist' is a list of pathways in igraph format.
# Output:  Writes each given pathway into a .txt file (named as the list elements) in PAL format, where
#          the first line is pathway name, then all nodes (Entrez) are listed one per line, and finally all 
#          relations (one per line) are listed in format "StartEntrez EndEntrez InteractionType", where
#          interaction type is either + or - for activation and inhibition, respectively.
WritePathways = function(plist){
  for(i in 1:length(plist)){
    
    # extract name, nodes, and edges from the pathway
    p = plist[[i]]
    pathname = names(plist)[i]
    nodes = V(p)$name
    rawedges = E(p)
    edgenumber = length(rawedges)
    if(edgenumber > 0){
      rawedges = t(as.matrix(as.data.frame(lapply(1:edgenumber, function(e){V(p)[inc(e)]$name}))))
      edges = apply(cbind(rawedges,rep("+",nrow(rawedges))), 1, paste, collapse=" ") 
    } else edges = NULL
    
    # Write into a file
    filename = paste(pathname, ".txt", sep="")
    write(c(pathname, as.character(nodes), as.character(edges)), file=filename, sep="\n")
  }
}


# Summary: Saves pathways as one gene set file in the format used by TcGSA
# Input:   'pathways' is a list of pathways in igraph format, 
#          'filename' is a name of the gmt file where the gene sets are written
# Output:  Writes the given pathways in gmt format used by TcGSA into the file, which name is given as an argument.
WritePathwayGmt = function(pathways, filename){
  
  # Extract first two cols (name and in this case artificial link)
  towrite = names(pathways)
  towrite = paste(towrite, rep("http://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_RIBOFLAVIN_METABOLISM", length(pathways)), sep="\t")
  
  # Add genes 
  for(i in 1:length(pathways)){
    towrite[i] = paste(towrite[i], paste(V(pathways[[i]])$name, collapse="\t"), sep="\t")
  }
  
  # Write
  write(towrite, file=filename, sep="\n")
}


# Summary: Identifies the least important (islands and leafs without many arriving relations) nodes for given pathways
# Input:   'paths' is a named list of pathways in igraph format
# Output:  A named list (pathways) of uninteresting nodes
IdentifyUninterestingNodes = function(paths){
  
  # Extract not-so-important nodes from each pathway at time
  leastimportantnodes = lapply(paths, function(p){
    
    # Extract the number of leaving/arriving relations for each node
    neighbournumber = degree(p)
    childnumber = degree(p, mode="out")
    pickindex = which((neighbournumber < 2) & (childnumber == 0))
    
    # If there are not-so-important nodes, return them. Otherwise return NULL.
    if(length(pickindex) > 0){
      toreturn = V(p)$name[pickindex]
    } else toreturn = NULL
    
    return(toreturn)
  })
  
  names(leastimportantnodes) = names(paths)
  return(leastimportantnodes)
}


#
#
#
# Output: a named list of two elements: ageexpression is a matrix of expression changes
#         caused by age (can be negative) and agecoef is a numeric vector telling the
#         mean age coefficient of each affected pathway
GenerateAgeExpression = function(agepathways, mineffect, maxeffect, baseexpression, 
                                 nonessentialnodes, pathways, sampleinfo, seed){
  
  # Generate mean (over genes) age effect for each affected pathway
  # NOTE: the first affected pathway gets the smallest coefficient and so on so that 
  # the biggest effect dominates when the gene-level values are generated
  set.seed(seed) 
  meaneffect = runif(length(agepathways), min=mineffect, max=maxeffect)
  meaneffect = sample(c(-1,1), size=length(agepathways), replace=T) * meaneffect
  meaneffect = meaneffect[order(abs(meaneffect), decreasing=F)]
  
  # Identify which unimportant nodes do not have the time trend within a trending pathway (if many unimportant genes, choose a random subset)
  set.seed(seed)
  agepathways_missnodes = mapply(function(lowimportance, path){
    pathnodes = V(path)$name
    noeffect = lowimportance
    if(length(lowimportance) > (0.25*length(pathnodes))){ 
      noeffect = sample(lowimportance, size=round(0.25*length(pathnodes)), replace=F)
    } 
    return(noeffect)
  }, lowimportance=nonessentialnodes[agepathways], path=pathways[agepathways])
  names(agepathways_missnodes) = agepathways
  
  # Initialise age-based expression
  ageexpression = matrix(0, nrow=nrow(baseexpression), ncol=ncol(baseexpression))
  colnames(ageexpression) = colnames(baseexpression)
  rownames(ageexpression) = rownames(baseexpression)
  
  # Extract donor and age info from sample info
  # For convenience, age is in matrix format, but all rows are the same
  donors = gsub("_.*", "", rownames(sampleinfo))
  donornmber = length(unique(donors))
  ages = matrix(rep(sampleinfo$Age, nrow(baseexpression)), byrow=T, ncol=ncol(data))
  samplenumber = nrow(sampleinfo)
  
  # Fill in age effect by pathways (least affected first)
  for(i in 1:length(agepathways)){
    
    # Extract affected genes from the pathway
    p = agepathways[i]
    pgenes = setdiff(V(pathways[[p]])$name, agepathways_missnodes[[p]])
    
    set.seed(i)
    donorfactors = rep(rnorm(donornumber, mean=0, sd=0.01), times=rle(donors)[[1]])
    donorfactors = matrix(rep(donorfactors, length(pgenes)), byrow=T, ncol=samplenumber)
    effectsizes = rnorm(length(pgenes), mean=meaneffect[i], sd=abs(0.2*meaneffect[i])) 
    effectsizes = matrix(rep(effectsizes, samplenumber), byrow=F, ncol=samplenumber)
    ageeffect = effectsizes + donorfactors
    ageexpression[pgenes,] = ages[1:length(pgenes),] * ageeffect * baseexpression[pgenes,]
  }
  
  return(list(ageexpression=ageexpression, agecoef=meaneffect))
}


#
#
#
# Output: a named list of two elements: diseaseexpression is a matrix of expression changes
#         caused by disease (can be negative) and diseasecoef is a numeric vector telling the
#         mean disease coefficient of each affected pathway
GenerateDiseaseExpression = function(diseasepathways, mineffect, maxeffect, baseexpression, 
                                     nonessentialnodes, pathways, sampleinfo, seed){
  
  # Generate mean (over genes) age effect for each affected pathway
  # NOTE: the first affected pathway gets the smallest coefficient and so on so that 
  # the biggest effect dominates when the gene-level values are generated
  set.seed(seed) 
  meaneffect = runif(length(diseasepathways), min=mineffect, max=maxeffect)
  meaneffect = sample(c(-1,1), size=length(diseasepathways), replace=T) * meaneffect
  meaneffect = meaneffect[order(abs(meaneffect), decreasing=F)]
  
  # Identify which unimportant nodes do not have the time trend within a trending pathway (if many unimportant genes, choose a random subset)
  set.seed(seed)
  diseasepathways_missnodes = mapply(function(lowimportance, path){
    pathnodes = V(path)$name
    noeffect = lowimportance
    if(length(lowimportance) > (0.25*length(pathnodes))){ 
      noeffect = sample(lowimportance, size=round(0.25*length(pathnodes)), replace=F)
    } 
    return(noeffect)
  }, lowimportance=nonessentialnodes[diseasepathways], path=pathways[diseasepathways])
  names(diseasepathways_missnodes) = diseasepathways
  
  # Initialise age-based expression
  diseaseexpression = matrix(0, nrow=nrow(baseexpression), ncol=ncol(baseexpression))
  colnames(diseaseexpression) = colnames(baseexpression)
  rownames(diseaseexpression) = rownames(baseexpression)
  
  # Extract donor and age info from sample info
  # For convenience, age is in matrix format, but all rows are the same
  donors = gsub("_.*", "", rownames(sampleinfo))
  donornmber = length(unique(donors))
  diseasestate = matrix(rep(sampleinfo$TimeFromEvent, nrow(baseexpression)), byrow=T, ncol=ncol(data))
  diseasestate[is.na(diseasestate)] = 0 # Control samples will have effect 0
  samplenumber = nrow(sampleinfo)
  
  # Fill in disease effect by pathways (least affected first)
  for(i in 1:length(diseasepathways)){
    
    # Extract affected genes from the pathway
    p = diseasepathways[i]
    pgenes = setdiff(V(pathways[[p]])$name, diseasepathways_missnodes[[p]])
    
    set.seed(i)
    donorfactors = rep(rnorm(donornumber, mean=0, sd=0.01), times=rle(donors)[[1]])
    donorfactors = matrix(rep(donorfactors, length(pgenes)), byrow=T, ncol=samplenumber)
    effectsizes = rnorm(length(pgenes), mean=meaneffect[i], sd=abs(0.2*meaneffect[i])) 
    effectsizes = matrix(rep(effectsizes, samplenumber), byrow=F, ncol=samplenumber)
    diseaseeffect = effectsizes + donorfactors
    diseaseexpression[pgenes,] = diseasestate[1:length(pgenes),] * diseaseeffect * baseexpression[pgenes,]
  }
  
  return(list(diseaseexpression=diseaseexpression, diseasecoef=meaneffect))
}

