\name{makeMetric}
\alias{makeMetric}

\title{
Create Weighted Cumulative Context Score (wCCS) To Rank miRNA Targeted Genes
}

\description{
For a cohort of regulated or interesting miRNA calculates the weighted cumulative context score (wCCS) for each predicted gene target.  
}

\usage{
makeMetric(miRGenes, deMirs)
}

\arguments{
  \item{miRGenes}{
A dataframe, ideally one element from the named list created using interactions function.
}
  \item{deMirs}{
A dataframe, ideally on element from the named list created by the upDown function (see details).
}
}

\details{
The miRGenes argument is ideally a named element from list object created by the the interactions function. The deMirs argument to makeMetric should be one element of the list created by the upDown function i.e. identifiers for regulated miRNA and the log2 fold change (or absolute value of log2 fold change for down-regulated miRNAs).  
}

\value{
A dataframe containing the Entrez gene ID, gene symbol and weighted cumulative context score (wCCS) for all genes predicted to be targets of miRNA detailed the deMirs object passed to the function.  
}

\references{Gallagher, I. J., Scheele, C., Keller, P., Nielsen, A. R., Remenyi, J., Fischer, C. P., Roder, K., Babraj, J., Wahlestedt, C., Hutvagner, G., Pedersen, P., Timmons, J.A. 2010 Integration of microRNA Changes in Vivo Identifies Novel Molecular Features of Muscle Insulin Resistance in Type 2 Diabetes. \emph{Genome Medicine} \bold{2}

Garcia, D. M., Baek, D., Shin, C., Bell, G. W., Grimson, A., Bartel, D. P. 2011 Weak Seed-pairing Stability and High Target-site Abundance Decrease the Proficiency of Lsy-6 and Other microRNAs. \emph{Nature Structural and Molecular Biology} \bold{18}, 1139--1146.
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
}

\keyword{ manip }
