interactions <- function(deMirsLst, mirInfo, genes){
  
  species <- substr(deMirsLst[[1]][1,1],1,3) # are these human mirs
  
  if (species =='hsa')
    mirInfo <- mirInfo[ which(mirInfo[,1] %in% genes[,1]), ] # if human get data for expressed eg ids
  else
    mirInfo <- mirInfo[ which(toupper(mirInfo[,2]) %in% toupper(genes[,1])), ] # if other get data for exp symbols
  
  upReg <- deMirsLst$up
  downReg <- deMirsLst$down
   
  deUpMirTgts <- mirInfo [which(mirInfo[,3] %in% upReg[,1]), ] # only regulated mir data

  deDownMirTgts <- mirInfo [which(mirInfo[,3] %in% downReg[,1]), ] # only regulated mir data

  return = list(upTgts = deUpMirTgts, downTgts = deDownMirTgts) 
  
}

