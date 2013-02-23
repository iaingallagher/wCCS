\name{createGeneSet}
\alias{createGeneSet}

\title{Create a GeneSet From miRNA Targeted Genes}

\description{Create R Geneset objects for analysis of Geneset enrichment in the genes predicted as miRNA targets.}
\usage{createGeneSet(df)}

\arguments{A dataframe of genes, ideally one element of the list produced by the interactions function.}

\details{This function uses the GSEABase package to create a GeneSetCollection object for use by the plotOverlaps function in this package or in other bioconductor packages which require this input. This GeneSetCollection object is made up of GeneSets with genes grouped by miRNA predicted to target them.}

\value{A GeneSetCollection object.}

\author{Iain J Gallagher}

\examples{
data(genes)
data(mirs)
data(TgtScanData)
deMirs <- upDown(mirs)
intExample <- interactions(deMirs, TgtScanData, genes)

geneSetEx <- createGeneSet(intExample$upTgts)
}

\keyword{manip}