getOverlaps <-
function(df1, df2){
    allOverlaps <- intersect(df1[,1], df2[,1]) # id overlaps
        
    overlapsOne <- df1[ which(df1[,1] %in% allOverlaps), ] # get overlap data
    overlapsTwo <- df2[ which(df2[,1] %in% allOverlaps), ]

# remove overlaps from original lists
    uniqueOne <- df1[ -which(df1[,1] %in% allOverlaps), ]
    uniqueTwo <- df2[ -which(df2[,1] %in% allOverlaps), ]
    
# calculate the wCCS ratios
    ratios <- overlapsOne[,3] / overlapsTwo[,3] # >1 belongs to overlapsOne, < 1 belongs to overlapsTwo
    
# calculate the new wCCS for each overlapping gene
    newScore <- pmin(overlapsOne[,3], overlapsTwo[,3]) - pmax(overlapsOne[,3], overlapsTwo[,3])
    
# make the new score the wCCS    
    overlapsOne[,3] <- newScore
    overlapsOne$ratio <- ratios
    
# which gene belongs to which list
    oneOnly <- overlapsOne[which(overlapsOne$ratio > 1),]
    twoOnly <- overlapsOne[which(overlapsOne$ratio < 1),]

# add as appropriate
    df1 <- rbind(df1, oneOnly[,-4])
    df2 <- rbind(df2, twoOnly[,-4])

# return data
    return = list(upTgts = df1, downTgts = df2)
}
