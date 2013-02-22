\name{getOverlaps}
\alias{getOverlaps}

\title{
Calculate a Revised Weighted Cumulative Context Score (wCCS) for Genes Targeted by Both Up and Downregulated miRNAs
}

\description{
Identifies those genes targeted by both up and down-regulated miRNA and assigns a new weighted cumulative context score (wCCS) based on whether the effect of the up or down-regulated miRNAs is predicted to be greater.
}

\usage{
getOverlaps(df1, df2)
}

\arguments{
 \item{df1}{Dataframe of genes targeted by up (or down) regulated miRNA, i.e. the  output of the makeMetric function.}
 \item{df2}{Dataframe of genes targeted by down (or up) regulated miRNA, i.e. the output of the makeMetric function.}
}

\details{
It is likely that some genes will be targeted by more than one miRNA of interest. In order to take this into account we examine the overlap between those genes targeted by up and down regulated miRNA and calculate a new cumulative metric for these based on whether they are more targeted by one or other set of regulated miRNAs. Overlaps (those genes targeted by up and down miRNA) are identified, a new wCCS is calculated by subtracting the higher (less negative) wCCS from the lower wCCS (more negative). Finally the gene is assigned as a predicted target of the up or down-regulated miRNAs based on the wCCS.
}

\value{
Returns a named list with elements upTgts and downTgts containing the Entrez gene ID, the gene symbol and the weighted cumulative context score (wCCS) for each predicted target of the supplied miRNA(s).
}

\author{
Iain J Gallagher
}

\examples{
data(genes)
data(mirs)
data(TgtScanData)

deMirs <- upDown(mirs)
intExample <- interactions(deMirs, TgtScanData, genes)

upCumulMetric <- makeMetric(intExample$upTgts, deMirs$up) 
downCumulMetric <- makeMetric(intExample$downTgts, deMirs$down)

cumulMetric <- getOverlaps(upCumulMetric, downCumulMetric)
}

\keyword{ manip }
