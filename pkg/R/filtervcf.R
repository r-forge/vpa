#Rscript for filter vcf data, only for one sample
filtervcf <- function(x, alter=NULL, alter.PL=20, QUAL=20, DP=10, GQ=NULL, FILTER=NULL, INDEL=NULL){
  if(class(x)!="vcf"){
    stop("x is not a vcf object")
  }
  if(!is.null(FILTER)){
    index1 <- x$FILTER %in% FILTER & as.numeric(x$QUAL)>=QUAL & as.numeric(x$INFO[,"DP"])>=DP
  }else{
    index1 <- as.numeric(x$QUAL)>=QUAL & as.numeric(x$INFO[,"DP"])>=DP
  }
  
  if(!is.null(GQ)){
    index2 <- as.numeric(x[[length(x)]][, "GQ"])>=GQ
  }else{
    index2 <- TRUE
  }
  if(!is.null(INDEL)){
    index3 <- x$INFO[, "INDEL"]==INDEL
  }else{
    index3 <- TRUE
  }
  if(!is.null(alter)){
    if(alter){
      # index4 <- x$REF != x$ALT
      PL <- x$SAMPLE1[,"PL"]
      #PL1 <- unlist(lapply(strsplit(PL, split=","), function(x)x[[1]]))
      PL <- lapply(strsplit(PL, split=","), function(x)x)
      index4 <- lapply(PL, function(x)as.numeric(x)[1]>=alter.PL | sum(as.numeric(x)[-1]<alter.PL)>0)
      index4[is.na(index4)] <- FALSE
      index4 <- unlist(index4) 
    }else{#not mutate
      PL <- x$SAMPLE1[,"PL"]
      #PL1 <- unlist(lapply(strsplit(PL, split=","), function(x)x[[1]]))
      #index4 <- as.numeric(PL1)<alter.PL
      PL <- lapply(strsplit(PL, split=","), function(x)x)
      index4 <- lapply(PL, function(x)as.numeric(x)[1]<alter.PL & sum(as.numeric(x)[-1]>=alter.PL)==length(x[-1]))
      index4[is.na(index4)] <- FALSE
      index4 <- unlist(index4)
    }
  }else{
    index4 <- TRUE
  }
  index <- index1 & index2 & index3 & index4
  index[is.na(index)] <- FALSE
  x$CHROM <- x$CHROM[index]
  x$POS <- x$POS[index]
  x$ID <- x$ID[index]
  x$REF <- x$REF[index]
  x$ALT <- x$ALT[index]
  x$QUAL <- x$QUAL[index]
  x$FILTER <- x$FILTER[index]
  x$FORMAT <- x$FORMAT[index]
  INFOID <- colnames(x$INFO)
  FORMATID <- colnames(x$SAMPLE1)
  x$INFO <- matrix(x$INFO[index,], nrow=sum(index))
  colnames(x$INFO) <- INFOID
  for(n in 1:(length(x)-10)){
    eval(parse(text=paste("x$SAMPLE",n," <- matrix(x$SAMPLE",n,"[index,], nrow=sum(index))", sep="")))
    eval(parse(text=paste("colnames(x$SAMPLE",n,") <- FORMATID", sep="")))
  }
  return(x)
}

subvcf <- function(x, CHR=NULL, POS=NULL, CHRPOS=NULL, samples=NULL, TF=NULL){
  if(class(x)!="vcf"){
    stop("x is not a vcf object")
  }
  if(!is.null(CHR)){
    index1 <- x$CHROM %in% CHR
  }else{
    index1 <- TRUE
  }
  if(!is.null(POS)){
    index2 <- x$POS %in% POS
  }else{
    index2 <- TRUE
  }
  if(!is.null(CHRPOS)){
    chrom <- paste(x$CHROM, x$POS, sep=":")
    index3 <- chrom %in% CHRPOS
  }else{
    index3 <- TRUE
  }
  if(!is.null(TF)){
    index4 <- TF
  }else{
    index4 <- TRUE
  }
  index <- index1 & index2 & index3 & index4
  index[is.na(index)] <- FALSE

  
  if(!is.null(samples)){
    s <- names(x)[11:length(x)] %in% samples
  }else{
    s <- rep(TRUE, length(x)-10)
  }
  s <- c(rep(TRUE, 10), s)
  
  x$CHROM <- x$CHROM[index]
  x$POS <- x$POS[index]
  x$ID <- x$ID[index]
  x$REF <- x$REF[index]
  x$ALT <- x$ALT[index]
  x$QUAL <- x$QUAL[index]
  x$FILTER <- x$FILTER[index]
  x$FORMAT <- x$FORMAT[index]
  INFOID <- colnames(x$INFO)
  FORMATID <- colnames(x$SAMPLE1)
#  x$INFO <- matrix(x$INFO[index,], nrow=sum(index))
#  colnames(x$INFO) <- INFOID
  x$INFO <- rbind(x$INFO[index,])
  for(n in 1:(length(x)-10)){
    #eval(parse(text=paste("x$SAMPLE",n," <- matrix(x$SAMPLE",n,"[index,], nrow=sum(index))", sep="")))
    eval(parse(text=paste("x$SAMPLE",n," <- rbind(x$SAMPLE",n,"[index,])", sep="")))
    eval(parse(text=paste("colnames(x$SAMPLE",n,") <- FORMATID", sep="")))
  }
  x <- x[s]
  class(x) <- "vcf"
  return(x)
}
