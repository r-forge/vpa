vcfreq <- function(vcf, method="fisher.test", ...){
  if(class(vcf)!="varlist")stop(paste(vcf, "is not a varlist object."))
  if(!is.null(vcf$VarVCF)){
    vcflist <- vcf$VarVCF
  }else if(!is.null(vcf$vcflist)){
    vcflist <- vcf$vcflist
  }else{
    stop("Input should be a varlist.")
  }
  pos <- lapply(vcflist, function(x)paste(x$CHROM, x$POS, sep=":"))
  Pos <- unique(unlist(pos))
  groups <- vcf$Samples[,2]
  ## #frequency
  ## countm <- matrix(unlist(lapply(pos, function(x)Pos %in% x)), ncol=length(pos))
  ## colnames(countm) <- names(pos)
  ## frequency <- t(apply(countm, 1, function(x)tapply(x, groups, function(y)sum(y)/length(y))))
  #REF
  ref <- lapply(vcflist, function(x)x$REF)
  ref1 <- unlist(ref)
  names(ref1) <- as.character(unlist(pos))
  gt <- lapply(vcflist, function(x)cbind(x$REF, x$ALT))
  REF <- ref1[match(Pos, names(ref1))]
  #ALT
  alt <- lapply(vcflist, function(x)x$ALT)
  GT <- lapply(vcflist, function(x)x$SAMPLE1[, "GT"])
  altm <- matrix(".", nrow=length(Pos), ncol=length(pos))
  for(i in 1:length(pos)){
    #altm[match(pos[[i]], Pos),i] <- alt[[i]]
    GTi <- GT[[i]]
    alti <- alt[[i]]
    REFi <- REF[match(pos[[i]], Pos)]
    Alti <- c()
    het <- GTi=="0/1"
    hom <- GTi=="1/1"
    Alti[het] <- paste(REFi[het], alti[het], sep="/")
    Alti[hom] <- paste(alti[hom], alti[hom], sep="/")    
    altm[match(pos[[i]], Pos), i] <- Alti
  }
  colnames(altm) <- names(pos)

  #alt allel frequecy
  gtm <- matrix("0/0", nrow=length(Pos), ncol=length(pos))
  for(j in 1:length(pos)){
    gtm[match(pos[[j]], Pos), j] <- GT[[j]]
  }
  acount <- apply(gtm, 1, function(x)tapply(x, groups, function(y){
    table(unlist(strsplit(y, split="/")))
  } ))
  frequency <- lapply(acount, function(x)lapply(x, function(y)
    ifelse("1" %in% names(y), y[names(y)=="1"]/sum(y), 0
           )))                    
  frequency <- matrix(unlist(frequency), byrow=TRUE, nrow=length(Pos))
  colnames(frequency) <- names(acount[[1]])

  #test
  p.value <- c()
  for(i in 1:length(Pos)){
    atable <- lapply(acount[[i]], function(x)x[match(c("1", "0"), names(x))])
    atable <- matrix(unlist(atable), ncol=2)
    atable[is.na(atable)] <- 0
    rownames(atable) <- c("1", "0")
    if(!is.na(pmatch(method, "fisher.test"))){
      p.value[i] <- fisher.test(atable, ...)$p.value
    }else if(!is.na(pmatch(method, "chisq.test"))){
      p.value[i] <- chisq.test(atable, ...)$p.value
    }else{
      stop("invalid method")
    }
  }
  
  freq <- cbind(REF, altm, frequency, p.value)
  rownames(freq) <- Pos
  freq <- freq[order(p.value),]
  return(freq)
}
