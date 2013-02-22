# testScript3.R
setwd('/home/iain/Documents/Work/Programming/R_scripts/wCCSPackage')
rm(list=ls())

# load the package, data and set output directory
# library(wCCS)
# load the functions
source('wCCS/R/upDown.R')
source('wCCS/R/interactions.R')
source('wCCS/R/geneCounts.R')
source('wCCS/R/makeMetric.R')
source('wCCS/R/getOverlaps.R')
source('wCCS/R/prioritiseTargets.R')

# load the data
load('wCCS/data/TgtScanData.rda')
dir.create('dataOut') # create a data for output

# load genes and regMirs
genes <- read.table('/home/iain/Documents/Work/Results/forOtherPeople/ditteAndersen/rescriptsetc/ratMuscleGenes.txt', header=TRUE, sep='\t')
mirs <- read.table('/home/iain/Documents/Work/Results/forOtherPeople/ditteAndersen/rescriptsetc/mirFCData.txt', header=TRUE, sep='\t')

# 1. Separate the regulated miRNA into those up or down regulated
deMirs <- upDown(mirs)

# 2. Get the interactions (mRNA/miRNA pairs) for the regulated miRNAs from the TgtScan data
interactionTest <- interactions(deMirs, TgtScanData, genes)

# 3. Get the counts of targets and targetting mirs
countsTest <- geneCounts(interactionTest)
write.table(countsTest$geneCounts, 'dataOut/hitsPerGene.txt', sep='\t', quote=F)
write.table(countsTest$mirCounts, 'dataOut/tgtsPerMir.txt', sep='\t', quote=F)

# 4. Get the cumulative metric for each targeted gene
upCumulMetric <- makeMetric(interactionTest$upTgts, deMirs$up) 
downCumulMetric <- makeMetric(interactionTest$downTgts, deMirs$down) 
 
# 5. Get mRNA targeted by up and down miRs
cumulMetric <- getOverlaps(upCumulMetric, downCumulMetric)

# 6. Output a figure showing the distribution of the wCCS
pdf('dataOut/wCCSDists.pdf', height=10, width=12, pointsize=12, paper='special')
par(mfrow = c(1,2))
hist(cumulMetric$upTgts[,3], xlab='wCCS', ylab='Freq', main='up wCCS Distribution', lwd=2, col='gray70')
hist(cumulMetric$downTgts[,3], xlab='wCCS', ylab='Freq', main='down wCCS Distribution', lwd=2, col='gray70')
dev.off()

png('dataOut/wCCSDists.png', height=700, width=800, pointsize=12)
par(mfrow = c(1,2), cex=1.5)
hist(upCumulMetric[,3], xlab='wCCS', ylab='Freq', main='up wCCS Distribution', lwd=2, col='gray70')
hist(downCumulMetric[,3], xlab='wCCS', ylab='Freq', main='down wCCS Distribution', lwd=2, col='gray70')
dev.off()

# 7. Prioritise the targets
upTgts <- cumulMetric$upTgts
downTgts <- cumulMetric$downTgts

priUps <- prioritiseTargets(upTgts)
priDowns <- prioritiseTargets(downTgts)

# 8. Write out data

write.table(priUps$all, 'dataOut/allTargetsUP.txt', sep='\t', quote=F, row.names=F)
write.table(priUps$top, 'dataOut/topTargetsUP.txt', sep='\t', quote=F, row.names=F)
write.table(priUps$bottom, 'dataOut/bottomTargetsUP.txt', sep='\t', quote=F, row.names=F)


write.table(priDowns$all, 'dataOut/allTargetsDOWN.txt', sep='\t', quote=F, row.names=F)
write.table(priDowns$top, 'dataOut/topTargetsDOWN.txt', sep='\t', quote=F, row.names=F)
write.table(priDowns$bottom, 'dataOut/bottomTargetsDOWN.txt', sep='\t', quote=F, row.names=F)

### TEST WITH ONLY UPREG MIRS ###
library(wCCS)
dir.create('dataOut', showWarnings = FALSE) # create a data for output
data(TgtScanData)
data(genes)
data(mirs)

# only up mirs
mirs <- mirs[-which(mirs[,2] < 0),]

# 1. Separate the regulated miRNA into those up or down regulated
deMirs <- upDown(mirs)

# 2. Get the interactions (mRNA/miRNA pairs) for the regulated miRNAs from the TgtScan data
interactionTest <- interactions(deMirs, TgtScanData, genes)

# 3. Get the counts of targets and targetting mirs
countsTest <- geneCounts(interactionTest)

# 4. Get the cumulative metric for each targeted gene
upCumulMetric <- makeMetric(interactionTest$upTgts, deMirs$up) 

# 5. Get the overlaps
# there will be no overlaps so silly to test

# 6. Output a figure showing the distribution of the wCCS
pdf('dataOut/wCCSDists.pdf', height=10, width=12, pointsize=12, paper='special')
par(cex=1.5)
hist(upCumulMetric[,3], xlab='wCCS', ylab='Freq', main='up wCCS Distribution', lwd=2, col='gray70')
dev.off()

# 7. Prioritise the targets
priUps <- prioritise(upCumulMetric)

# 8. Write out data
write.table(priUps$all, 'dataOut/allTargetsUP.txt', sep='\t', quote=F, row.names=F)
write.table(priUps$top, 'dataOut/topTargetsUP.txt', sep='\t', quote=F, row.names=F)
write.table(priUps$bottom, 'dataOut/bottomTargetsUP.txt', sep='\t', quote=F, row.names=F)

# make a GeneSetCollection for further analysis
upGeneSet <- geneSetFunc(interactionTest$up)
downGeneSet <- geneSetFunc(interactionTest$down)

# 9. Are the genes targted by e.g. miR-495 enriched for any genesets

