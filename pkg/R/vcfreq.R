vcfreq <- function(x){
  if(class(x)!="varlist")stop(paste(x, "is not a varlist object."))
  if(!is.null(x$VarVCF)){
    vcflist <- x$VarVCF
  }else if(!is.null(x$vcflist)){
    vcflist <- x$vcflist
  }else{
    stop("Input should be a varlist.")
  }
  pos <- lapply(vcflist, function(x)paste(x$CHROM, x$POS, sep=":"))
  Pos <- unique(unlist(pos))
  #frequency
  countm <- matrix(unlist(lapply(pos, function(x)Pos %in% x)), ncol=length(pos))
  colnames(countm) <- names(pos)
  groups <- x$Samples[,2]
  frequency <- t(apply(countm, 1, function(x)tapply(x, groups, function(y)sum(y)/length(y))))
  #REF
  ref <- lapply(vcflist, function(x)x$REF)
  ref1 <- unlist(ref)
  names(ref1) <- as.character(unlist(pos))
  gt <- lapply(vcflist, function(x)cbind(x$REF, x$ALT))
  REF <- ref1[match(Pos, names(ref1))]
  #ALT
  alt <- lapply(vcflist, function(x)x$ALT)

  altm <- matrix("-", nrow=length(Pos), ncol=length(pos))
  for(i in 1:length(pos)){
    altm[match(pos[[i]], Pos),i] <- alt[[i]]
  }
  colnames(altm) <- names(pos)
  
  freq <- cbind(REF, altm, frequency)
  rownames(freq) <- Pos
  return(freq)
}
