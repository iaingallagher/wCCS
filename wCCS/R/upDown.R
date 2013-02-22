upDown <- function(x){

	downInd <- which(x[, 2] < 0)
	down <- x[downInd, ]
    	down[,2] <- abs(down[,2])
	upInd <- which(x[, 2] > 0)
	up <- x[upInd,]
	
	return = list(up=up, down=down)
}

