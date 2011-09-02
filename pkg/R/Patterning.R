Patterning <- function(x, pattern, paired=FALSE, not.covered=NULL, var.PL=NULL){
  if(sum(c("vcflist", "VarIndex", "Samples") %in% names(x))<3){
    stop("Please check input data!")
  }
  samples <- x$Samples
  if(sum(unique(samples[,2]) %in% colnames(pattern))<length(unique(samples[,2]))){
    stop("Pattern matrix don't match with sample groups!")
  }
  vcflist <- x$vcflist
  vari <- x$VarIndex
  if(sum(names(vcflist)==samples[,1] & colnames(vari)==samples[,1])!=length(samples[,1])){
    stop("Samples names should be concordant!")
  }  
  #variant frequency
  if(is.null(not.covered)){
    vari <- na.omit(vari)
  }else if(not.covered==TRUE){
    vari[is.na(vari)] <- TRUE
  }else if(not.covered==FALSE){
    vari[is.na(vari)] <- FALSE
  }
  #PL
  if(is.null(var.PL)){
    vari <- vari
  }else if(var.PL==TRUE){
    vari[vari=="PL"] <- TRUE
  }else if(var.PL==FALSE){
    vari[vari=="PL"] <- FALSE
  }
  vari <- vari==TRUE
  
  AltF <- t(apply(vari, 1, function(x)tapply(x, samples[,2], function(x)sum(mean(x)))))
  pattern <- pattern[,match(colnames(AltF), colnames(pattern))]
  
  #pattern
  if(paired){
    coni <- grep("con", samples[,2], ignore.case=TRUE)
    if(length(coni)==0){
      stop("Please specify control group in label column.")
    }
    control <- unique(samples[coni,2])
    casei <- c(1:nrow(samples))[-coni]
    
    for(n in coni){
      ci <- which(samples[,1] %in% samples[n,1])
      ci <- setdiff(ci, coni)
      vari[,ci] <- vari[,ci]==TRUE & vari[,n]==FALSE
    }
    AltF <- t(apply(vari, 1, function(x)tapply(x, samples[,2], function(x)sum(mean(x)))))
    fns <- colnames(AltF)
    AltF <- AltF[,colnames(AltF)!=control]
    AltF <- cbind(AltF)
    colnames(AltF) <- fns[fns!=control]
    pns <- colnames(pattern)
    pattern <- pattern[,colnames(pattern)!=control]
    pattern <- cbind(pattern)
    colnames(pattern) <- pns[pns!=control]
    AltP <- apply(AltF, 1, function(x)x >= pattern[1,] & x <= pattern[2,])
    AltP <- t(rbind(AltP))
    AltP <- rowSums(AltP)==ncol(AltP)
    VarI <- vari[,casei] & AltP
  }else{
    AltP <- t(apply(AltF, 1, function(x)x >= pattern[1,] & x <= pattern[2,]))
    AltP <- rowSums(AltP)==ncol(AltP)
    VarI <- vari & AltP
  }

  pos <- rownames(VarI)
  VarVCF <- list()
  if(paired){
    for(i in 1:ncol(VarI)){
      VarVCF[[i]] <- subvcf(vcflist[casei][[i]], CHRPOS=pos[VarI[,i]])
    }
    names(VarVCF) <- samples[casei, 1]
  }else{
    for(i in 1:ncol(VarI)){
      VarVCF[[i]] <- subvcf(vcflist[[i]], CHRPOS=pos[VarI[,i]])
    }
    names(VarVCF) <- samples[,1]
  }

  PVCF <- list(VarVCF=VarVCF, VarFrequency=AltF, Pattern=pattern, Samples=samples)
  class(PVCF) <- "varlist"
  return(PVCF)
}
