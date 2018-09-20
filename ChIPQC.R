#Specify genome
#Runs ChIPQC on specified working directory ($1) and genome ($2). Working directory should include sampleSheet.
library(ChIPQC)

args <- commandArgs(TRUE)

wd <- as.double(args[1])
setwd(wd)

#Specify genome
gnm <- as.double(args[2])

#Specify blacklist file
bl <- "/mnt/iscsi_speed/blelloch/deniz/genomes/blacklist/mm10.blacklist.bed"

#Load samplesheet (containing columns SampleID  Tissue  Factor
#Condition  Treatment Replicate bamReads  ControlID bamControl  Peaks PeakCaller)
samples <- dir(wd,"_sampleSheet.csv")
samples <- read.csv(file.path(wd,samples[1]))

#Construct ChIPQCexperiment object

exp <- ChIPQC(samples, consensus = T, bCount = T, summits = 250,
               annotation = gnm, blacklist = b)

ChIPQCreport(exp,facetBy=c("Tissue","Treatment"))
