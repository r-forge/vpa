#Rscript for filter vcf data, only for one sample
filtervcf <- function(vcf, alter=NULL, alter.PL=20, alter.AD=3, alter.ADP=NULL, QUAL=20, DP=c(10,500), GQ=NULL, FILTER=NULL, INDEL=NULL){
  if(class(vcf)!="vcf"){
    stop("The input is not a vcf object")
  }
  if(!is.null(FILTER)){
    index1 <- vcf$FILTER %in% FILTER
  }else{
    index1 <- TRUE
  }
  if(!is.null(QUAL)){
    index2 <- as.numeric(vcf$QUAL)>=QUAL
  }else{
    index2 <- TRUE
  }
  if(!is.null(DP)){
    if(length(DP)==2){
      index3 <- as.numeric(vcf$INFO[,"DP"])>=DP[1] & as.numeric(vcf$INFO[,"DP"])<=DP[2]
    }else{
      stop("DP should be a two numeric vector including the minimum and maximum depths.")
    }
  }else{
    index3 <- TRUE
  }
  if(!is.null(GQ)){
    index4 <- as.numeric(vcf[[length(vcf)]][, "GQ"])>=GQ
  }else{
    index4 <- TRUE
  }
  if(!is.null(INDEL)){
    index5 <- vcf$INFO[, "INDEL"]==INDEL
  }else{
    index5 <- TRUE
  }
  if(!is.null(alter)){
    if(!is.null(alter.AD) | !is.null(alter.ADP)){ #depth of alt allele
      if("AD" %in% colnames(vcf$SAMPLE1)){
        AD <- vcf$SAMPLE1[, "AD"]
        AD <- unlist(lapply(strsplit(AD, split=","), function(x)x[2]))
      }else if("DP4" %in% colnames(vcf$INFO)){
        DP4 <- vcf$INFO[, "DP4"]
        AD <- unlist(lapply(strsplit(DP4, split=","), function(x)sum(as.integer(x)[3:4])))
      }else{
        warning("no depth of alter allele found!")
      }
    }
    if(alter){
      # index6 <- vcf$REF != vcf$ALT
      PL <- vcf$SAMPLE1[,"PL"]
      #PL1 <- unlist(lapply(strsplit(PL, split=","), function(x)x[[1]]))
      PL <- lapply(strsplit(PL, split=","), function(x)x)
      index6 <- lapply(PL, function(x)as.numeric(x)[1]>=alter.PL | sum(as.numeric(x)[-1]<alter.PL)>0)
      index6[is.na(index6)] <- FALSE
      index6 <- unlist(index6)
      if(!is.null(alter.AD)){
        index7 <- as.numeric(AD)>=alter.AD
        index7[is.na(index7)] <- FALSE
      }else{
        index7 <- TRUE
      }
      if(!is.null(alter.ADP)){
        index8 <- as.numeric(AD)/as.numeric(vcf$INFO[,"DP"])>=alter.ADP
        index8[is.na(index8)] <- FALSE
      }else{
        index8 <- TRUE
      }
      index6 <- index6 & index7 & index8
    }else{#not mutate
      PL <- vcf$SAMPLE1[,"PL"]
      #PL1 <- unlist(lapply(strsplit(PL, split=","), function(x)x[[1]]))
      #index6 <- as.numeric(PL1)<alter.PL
      PL <- lapply(strsplit(PL, split=","), function(x)x)
      index6 <- lapply(PL, function(x)as.numeric(x)[1]<alter.PL & sum(as.numeric(x)[-1]>=alter.PL)==length(x[-1]))
      index6[is.na(index6)] <- FALSE
      index6 <- unlist(index6)
      if(!is.null(alter.AD)){
        index7 <- as.numeric(AD)<=alter.AD
        index7[is.na(index7)] <- TRUE
      }else{
        index7 <- TRUE
      }
      if(!is.null(alter.ADP)){
        index8 <- as.numeric(AD)/as.numeric(vcf$INFO[,"DP"])<=alter.ADP
        index8[is.na(index8)] <- TRUE
      }else{
        index8 <- TRUE
      }
      index6 <- index6 & index7 & index8
    }
  }else{
    index6 <- TRUE
  }
  index <- index1 & index2 & index3 & index4 & index5 & index6
  index[is.na(index)] <- FALSE

  filteredvcf <- subvcf(vcf, TF=index)
  droppedvcf <- subvcf(vcf, TF=!index)
  return(list(filtered=filteredvcf, dropped=droppedvcf))
  
 ##  vcf$CHROM <- vcf$CHROM[index]
 ##  vcf$POS <- vcf$POS[index]
 ##  vcf$ID <- vcf$ID[index]
 ##  vcf$REF <- vcf$REF[index]
 ##  vcf$ALT <- vcf$ALT[index]
 ##  vcf$QUAL <- vcf$QUAL[index]
 ##  vcf$FILTER <- vcf$FILTER[index]
 ##  vcf$FORMAT <- vcf$FORMAT[index]
 ##  #INFOID <- colnames(vcf$INFO)
 ##  #FORMATID <- colnames(vcf$SAMPLE1)
 ##  vcf$INFO <- cbind(vcf$INFO[index,])
 ## # colnames(vcf$INFO) <- INFOID
 ##  for(n in 1:(length(vcf)-10)){
 ##    eval(parse(text=paste("vcf$SAMPLE",n," <- cbind(vcf$SAMPLE",n,"[index,])", sep="")))
 ##    #eval(parse(text=paste("colnames(vcf$SAMPLE",n,") <- FORMATID", sep="")))
 ##  }
 ##  return(vcf)
}

subvcf <- function(vcf, CHR=NULL, POS=NULL, CHRPOS=NULL, samples=NULL, TF=NULL){
  if(class(vcf)!="vcf"){
    stop("The input is not a vcf object")
  }
  if(!is.null(CHR)){
    index1 <- vcf$CHROM %in% CHR
  }else{
    index1 <- TRUE
  }
  if(!is.null(POS)){
    index2 <- vcf$POS %in% POS
  }else{
    index2 <- TRUE
  }
  if(!is.null(CHRPOS)){
    chrom <- paste(vcf$CHROM, vcf$POS, sep=":")
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
    s <- names(vcf)[11:length(vcf)] %in% samples
  }else{
    s <- rep(TRUE, length(vcf)-10)
  }
  s <- c(rep(TRUE, 10), s)
  
  vcf$CHROM <- vcf$CHROM[index]
  vcf$POS <- vcf$POS[index]
  vcf$ID <- vcf$ID[index]
  vcf$REF <- vcf$REF[index]
  vcf$ALT <- vcf$ALT[index]
  vcf$QUAL <- vcf$QUAL[index]
  vcf$FILTER <- vcf$FILTER[index]
  vcf$FORMAT <- vcf$FORMAT[index]
  INFOID <- colnames(vcf$INFO)
#  FORMATID <- colnames(vcf$SAMPLE1)
#  vcf$INFO <- matrix(vcf$INFO[index,], nrow=sum(index))
#  colnames(vcf$INFO) <- INFOID
  vcf$INFO <- cbind(vcf$INFO[index,])
  for(n in 1:(length(vcf)-10)){
    #eval(parse(text=paste("vcf$SAMPLE",n," <- matrix(vcf$SAMPLE",n,"[index,], nrow=sum(index))", sep="")))
    eval(parse(text=paste("vcf$SAMPLE",n," <- cbind(vcf$SAMPLE",n,"[index,])", sep="")))
 #   eval(parse(text=paste("colnames(vcf$SAMPLE",n,") <- FORMATID", sep="")))
  }
  vcf <- vcf[s]
  class(vcf) <- "vcf"
  return(vcf)
}
