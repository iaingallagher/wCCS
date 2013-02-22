\name{prioritiseTargets}
\alias{prioritiseTargets}

\title{
Prioritise the Targets of miRNAs
}

\description{
The gene targets of miRNAs are prioritised based on the weighted cumulative context score (wCCS).
}

\usage{
prioritiseTargets(df, quant = 0.25)
}


\arguments{
  \item{df}{A dataframe of regulated miRNA (see details).}
  \item{quant}{The quantile cutoff to decide which mRNA are the best / worst predicted targets of the miRNA of interest. The default quant argument is 0.25 and it should be set > 0 and < 1.}
}

\details{
The df argument is a dataframe containing the information returned from either the getOverlaps function (in the case of up and down regulated miRNA) or the makeMetric function (if only miRNA regulated in one direction are being examined). The quant argument is used to select the top quant\% gene targets by wCCS. The default is to select the top and bottom 25\%. 
}

\value{
A named list is returned containing
  \item{all}{A dataframe containing the all the results from the analysis.}
  \item{top}{The genes which are the best predicited targets of the miRNA. This is decided by the quant argument and consists of those genes where the wCCS is in the top quant\% of all the wCCS.}
  \item{bottom}{The genes which are the worst predicted targets of the miRNA. This is decided by the quant argument and consists of those genes where the wCCS is in the bottom quant\% of all the wCCS.}
}

\references{
Gallagher, I. J., Scheele, C., Keller, P., Nielsen, A. R., Remenyi, J., Fischer, C. P., Roder, K., Babraj, J., Wahlestedt, C., Hutvagner, G., Pedersen, P., Timmons, J.A. 2010 Integration of microRNA Changes in Vivo Identifies Novel Molecular Features of Muscle Insulin Resistance in Type 2 Diabetes. \emph{Genome Medicine} \bold{2}
}

\author{
Iain J Gallagher
}

\examples{
# fake up some data
data(genes)
data(mirs)
data(TgtScanData)

deMirs <- upDown(mirs)
intExample <- interactions(deMirs, TgtScanData, genes)

upCumulMetric <- makeMetric(intExample$upTgts, deMirs$up) 
downCumulMetric <- makeMetric(intExample$downTgts, deMirs$down)

cumulMetric <- getOverlaps(upCumulMetric, downCumulMetric)

upTgts <- cumulMetric$upTgts
downTgts <- cumulMetric$downTgts

priUp <- prioritiseTargets(upTgts)
priDown <- prioritiseTargets(downTgts)
}

\keyword{ manip }
