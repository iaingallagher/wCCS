geneCounts <- function(lst){

    df1 <- lst[[1]]
    df2 <- lst[[2]]
    dfAll <- rbind(df1, df2)

	mirCounts <- as.data.frame(sapply(unstack(dfAll, Gene.Symbol~miRNA), length))
    colnames(mirCounts) <- c('numberTargets')
    geneCounts <- as.data.frame(sapply(unstack(dfAll, miRNA~Gene.Symbol), length))
    colnames(geneCounts) <- c('numberHits')

	return = list(geneCounts = geneCounts, mirCounts = mirCounts)
}
