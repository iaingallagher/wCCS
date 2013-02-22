\name{interactions}
\alias{interactions}

\title{
Selects miRNA Targets From TargetScan Data
}

\description{
Using the miRNA and gene data provided parses the TargetScan data for specific miRNA/gene interactions.  
}
\usage{
interactions(deMirsLst, mirInfo, genes)
}

\arguments{
\item{deMirsLst}{
A named list of regulated miRNA i.e. the output of getUpDown (see details).
}
\item{mirInfo}{
A dataframe from TargetScan data containing details of mRNA/miRNA interactions.
}
\item{genes}{
A character vector of genes of interest (see details).
}
}

\details{
The identifier for the miRNAs should be in the same form as the miRNA labels in TargetScan e.g. hsa-miR-451. The case matching is important: hsa-mir-451 will not work. Data in the correct form is included in the package. Entrez gene IDs should be used as the gene identifiers. The package uses the conserved data from TargetScan and the Entrez gene identifiers included in the TargetScan data are human. For human miRNA studies the package uses Entrez gene identifiers to map the gene / miRNA interactions but for other species it uses gene symbols. 
}

\value{
A named list with elements upTgts and downTgts containing the genes from the TargetScan database interacting with the up-regulated or down-regulated miRNAs of interest respectively. Each element of the list contains the Entrez ID, the gene symbol, interacting miRNA and the TargetScan context score for each of the predicted miRNA targets.  
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
}

\keyword{ manip }

