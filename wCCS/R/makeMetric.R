makeMetric <-
function(miRGenes, deMirs){
	
# need to match position of each miR in deMirPresGenes with its FC to form a vector of FC in correct order
	
	fcVector <- as.numeric(with (deMirs, FC[match(miRGenes[,3], Probe)] ) )

	# multiply fc by context score for each interaction	
    metric <- fcVector * as.numeric(as.character(miRGenes[,4]))

    geneMetric <- cbind(miRGenes[,1], metric)
    colnames(geneMetric) <- c('egID', 'wCCS')

    # make cumul by aggregate
    listMetric <- aggregate(as.numeric(geneMetric[,2]), list(geneMetric[,1]), sum)#returns a dataframe

    listMetric$symbol <- miRGenes [with(listMetric, match(listMetric[,1],miRGenes[,1])),2]
    listMetric <- listMetric[ ,c(1,3,2)] # reorder for tidyness

    colnames(listMetric) <- c('egID','symbol', 'wCCS')
    
	# return whole list
	return = (listMetric)# dataframe
}

