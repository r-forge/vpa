#Rscript for writing to vcf file.
write.vcf <- function(x, file="", HEAD=TRUE, ...)UseMethod("write.vcf")
write.vcf.default <- function(x, file="", HEAD=TRUE, ...){
  #INFO
  INFO <- x$INFO
  INFOc <- c()
  for(i in 1:nrow(x$INFO)){
    #info <- c()
    info1 <- na.omit(x$INFO[i,])
    if("INDEL" %in% names(info1)){
      if(info1["INDEL"]==TRUE){
        info <- "INDEL"
      }else{
        info <- NA
      }
      info1 <- info1[-grep("INDEL", names(info1))]
      info1n <- names(info1)
      info <- na.omit(c(info, paste(info1n, info1, sep="=")))
    }else{
      info1n <- names(info1)
      if(!is.null(info1n)){
        info <- paste(info1n, info1, sep="=")
      }else{
        info <- info1
      }
    }
    info <- paste(info, collapse=";")
    INFOc <- c(INFOc, info)
  }
  #SAMPLE
  SAMPLEC <- c()
  for(n in 1:(length(x)-10)){
    eval(parse(text=paste("SAMPLE <- x$SAMPLE", n, sep="")))
    SAMPLE <- rbind(SAMPLE)
    SAMPLEc <- sapply(1:nrow(SAMPLE), function(y)paste(SAMPLE[y, ][match(unlist(strsplit(x$FORMAT[y], split=":")), names(SAMPLE[y,]))], collapse=":"))
    #SAMPLEc <- apply(SAMPLE, 1, function(x)paste(na.omit(x), collapse=":"))
    SAMPLEC <- cbind(SAMPLEC, SAMPLEc)
  }
  vcfdata <- cbind(x$CHROM, x$POS, x$ID, x$REF, x$ALT, x$QUAL, x$FILTER, INFOc, x$FORMAT, SAMPLEC)
  if(HEAD){
    colnames(vcfdata) <- strsplit(tail(x$HEAD, n=1), split="\t|#")[[1]][-1]
  }else{
    colnames(vcfdata) <- NULL
  }
  if (file == ""){
    return(vcfdata)
  }
  if(HEAD){
    write(x$HEAD, file, sep="\t")
    write.table(vcfdata, file=file, append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, ...)
  }else{
    write.table(vcfdata, file=file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, ...)
  }
}
