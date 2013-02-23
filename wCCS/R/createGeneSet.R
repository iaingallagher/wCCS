# function to generate a GeneSetCollection for each list of miRNA targets
# input - df - input - df - select either upTgts or downTgts from output of interactions.R
# output - geneset collection for further analysis (e.g. GSEA)

createGeneSet <- function(listIn){

	geneMiRList <- split(listIn[,2], listIn[,3]) # make list - names are miRs; elements are gene symbols
	uniqueList <- lapply(geneMiRList, unique)# need unique values in list elements to make genesets

    sets <- Map(GeneSet, uniqueList, setName=names(uniqueList),
                MoreArgs=list(geneIdType = SymbolIdentifier() ))
    GeneSetCollection(sets)
}

# this doesn't work at the moment because of extra factors in TgtScanData - need to recreate TgtScanData
# and drop those