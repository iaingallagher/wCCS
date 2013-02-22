\name{upDown}
\alias{upDown}

\title{
Separate Up and Down-regulated miRNAs
}
\description{
Generates separate lists of up and down-regulated miRNA from an input dataframe.  
}
\usage{
upDown(x)
}

\arguments{
\item{x}{A dataframe of regulated miRNA (see details).}
}

\details{
The first column of the input data must be headed 'Probe' and is the miRNA identifer. The second column must be headed 'FC' and is fold change data. The fold change for each miRNA should be on the log2 scale. This ensures the calculation of the weighted cumulative context score (wCCS) is equivalent for up and down regulated miRNA and the results are comparable. The absolute value of the fold change for the downregulated miRNAs is returned in the resulting list and is used in the calculation of the wCCS.  
}

\value{A named list with the following elements:
\item{up}{The identifier and fold change for up-regulated miRNA}
\item{down}{The identifier and (absolute) fold change for down-regulated miRNA}
}


\author{
Iain J Gallagher
}

\examples{
data(mirs)
deMirs <- upDown(mirs)
}
\keyword{ manip }
