\name{TgtScanData}
\docType{data}

\alias{TgtScanData}

\title{Conserved Site Data From the TargetScan Database (ver 6.2)}

\description{
This data includes the Entrez gene ID, the gene symbol, the interacting miRNA and the context score. These are the data required by the wCCS functions. This data represents the conserved site data only.
}

\usage{data(TgtScanData)}

\format{
A dataframe consisting of four columns, GeneID (the Entrez gene ID), Symbol (the gene symbol), miRNA (the miRNA identifier (e.g. hsa-miR-451)) and context.score (the TargetScan the context score for the particular mRNA/miRNA interaction).
}

\source{\href{www.targetscan.org}{TargetScan}}

\references{Garcia, D. M., Baek, D., Shin, C., Bell, G. W., Grimson, A., Bartel, D. P. 2011 Weak Seed-pairing Stability and High Target-site Abundance Decrease the Proficiency of Lsy-6 and Other microRNAs. \emph{Nature Structural and Molecular Biology} \bold{18}, 1139--1146.
}

\keyword{ datasets }
