\name{geneCounts}
\alias{geneCounts}

\title{
Count miRNA Predicted to Target Genes and Predicted Targets of Each miRNA.
}

\description{
Counts the number of miRNA predicted to target each gene and the number of predicted targets for each miRNA.
}
\usage{
geneCounts(lst)
}

\arguments{
  \item{lst}{
A named list containing the targets of up and down regulated miRNA i.e. the output of the interactions function.
}
}

\details{
This function simply counts the number of miRNA targeting each gene identified as a potential target from the supplied list and the number of mRNA each miRNA is predicted to target. If there is no data in TargetScan for a miRNA included in the list of regulated miRNA then no data is returned.
}

\value{
Returns a named list containing two elements. geneCounts - the gene symbol and the number of miRNA predicted to target that gene. mirCounts - the number of genes in the supplied list that each miRNA is predicted to target.
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
counts <- geneCounts(intExample)
}

\keyword{ manip }

