prioritiseTargets <-
function (df, quant=0.25){

	topTargets <- df [ which(df[,3] < quantile(df[,3], quant)), ]
	bottomTargets <- df [ which(df[,3] > quantile(df[,3], 1-quant)), ]

    
	# return whole list
	return = list(all = df, top = topTargets, bottom = bottomTargets)# dataframe
}
