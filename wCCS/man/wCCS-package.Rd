\name{wCCS-package}
\alias{wCCS-package}

\alias{wCCS}

\docType{package}

\title{
Prioritise Gene Targets of Several MicroRNA
}

\description{
For predicted mRNA/miRNA interactions the \href{www.targetscan.org}{TargetScan} database provides a context score indicating the confidence in any particular interaction. Using this context score and  information about miRNA abundance from e.g. a microarray or qPCR survey the wCCS package prioritises the predicted targets of all selected miRNAs.
}

\details{
\tabular{ll}{
Package: \tab wCCS\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2012-07-27\cr
License: \tab GPL-3\cr
}
}

\author{
Iain J Gallagher

Maintainer: Iain J Gallagher <iaingallagher@gmail.com>
}

\references{Gallagher, I. J., Scheele, C., Keller, P., Nielsen, A. R., Remenyi, J., Fischer, C. P., Roder, K., Babraj, J., Wahlestedt, C., Hutvagner, G., Pedersen, P., Timmons, J.A. 2010 Integration of microRNA Changes in Vivo Identifies Novel Molecular Features of Muscle Insulin Resistance in Type 2 Diabetes. \emph{Genome Medicine} \bold{2}

Garcia, D. M., Baek, D., Shin, C., Bell, G. W., Grimson, A., Bartel, D. P. 2011 Weak Seed-pairing Stability and High Target-site Abundance Decrease the Proficiency of Lsy-6 and Other microRNAs. \emph{Nature Structural and Molecular Biology} \bold{18}, 1139--1146.
}

\keyword{ package }

\examples{
# load the required data
data(TgtScanData)
data(genes)
data(mirs)

# create a directory to write data out
dir.create('dataOut', showWarnings = FALSE)
# 1. Separate the regulated miRNA into those up or down regulated
deMirs <- upDown(mirs)
# 2. Get the interactions (mRNA/miRNA pairs) for the regulated miRNAs from the TgtScan data
interactionTest <- interactions(deMirs, TgtScanData, genes)
# 3. Get the counts of targets and targeting mirs
countsTest <- geneCounts(interactionTest)
# 4. Get the cumulative metric for each targeted gene
upCumulMetric <- makeMetric(interactionTest$upTgts, deMirs$up) 
downCumulMetric <- makeMetric(interactionTest$downTgts, deMirs$down) 
# 5. Re-score targeted by up and down miRs
olaps <- getOverlaps(upCumulMetric, downCumulMetric)
# 7. Prioritise the targets
ups <- olaps$upTgts
downs <- olaps$downTgts
priUps <- prioritiseTargets(ups)
priDowns <- prioritiseTargets(downs)
# 8. Write out data
write.table(priUps$all, 'dataOut/allTargetsUP.txt', sep='\t', quote=FALSE, row.names=FALSE)
write.table(priUps$top, 'dataOut/topTargetsUP.txt', sep='\t', quote=FALSE, row.names=FALSE)
write.table(priUps$bottom, 'dataOut/bottomTargetsUP.txt', sep='\t', quote=FALSE, row.names=FALSE)
write.table(priDowns$all, 'dataOut/allTargetsDOWN.txt', sep='\t', quote=FALSE, row.names=FALSE)
write.table(priDowns$top, 'dataOut/topTargetsDOWN.txt', sep='\t', quote=FALSE, row.names=FALSE)
write.table(priDowns$bottom, 'dataOut/bottomTargetsDOWN.txt', sep='\t', quote=FALSE, row.names=FALSE)
}
