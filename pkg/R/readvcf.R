#Rscript for read vcf format
#file="/user/songliu/u2/group/Qiang/Seq/output/Exome/mpileup/1151HZ0001.flt.vcf"
read.vcf <- function(file, VCF=NULL, INFOID=NULL, FORMATID=NULL, ...)UseMethod("read.vcf")
read.vcf.default <- function(file, VCF=NULL, INFOID=NULL, FORMATID=NULL, ...){
  if(is.null(VCF)){
    cat("Read data...\n")
    vcf <- read.table(file, comment.char="#", sep="\t")
    vcfhead <- readLines(file, n=100)
    vcfhead <- vcfhead[grep("^#", vcfhead)]
  }else{
    vcf <- VCF$vcf
    vcfhead <- VCF$header
  }
  cat("Format to vcf object...\n")
  vcfdata <- list(
                  HEAD=vcfhead,
                  CHROM=as.character(vcf[,1]),
                  POS=as.integer(vcf[,2]),
                  ID=as.character(vcf[,3]),
                  REF=as.character(vcf[,4]),
                  ALT=as.character(vcf[,5]),
                  QUAL=as.numeric(vcf[,6]),
                  FILTER=as.character(vcf[,7]),
                  INFO=as.character(vcf[,8]),
                  FORMAT=as.character(vcf[,9])
                  )
  #titles <- strsplit(tail(vcfhead,1), split="\t")[[1]]
  #samples <- titles[10:length(titles)]
  for(n in 1:(ncol(vcf)-9)){
    assign(paste("SAMPLE", n, sep=""), as.character(vcf[,n+9]))
    eval(parse(text=paste("vcfdata$SAMPLE",n," <- SAMPLE",n, sep="")))
  }

  #INFO
  if(is.null(INFOID)){
    INFOh <- vcfhead[grep("##INFO", vcfhead)]
    INFOID <- unlist(lapply(strsplit(INFOh, split="=|,"), function(x)x[3]))
  }
  info1 <- strsplit(vcfdata$INFO, split=";")
  
  info1n <- lapply(info1, function(a)sapply(a, function(x)unlist(strsplit(x, split="="))[[1]]))
  info1v <- lapply(info1, function(a)sapply(a, function(x){inn <- unlist(strsplit(x, split="=")); ifelse(length(inn)>1, inn[2], TRUE)}))

  info2 <- sapply(1:length(info1), function(x)info1v[[x]][match(INFOID, info1n[[x]])])
  info2 <- t(rbind(info2))
  colnames(info2) <- INFOID
  
##   info2 <- lapply(info1, function(x)x[match(INFOID, x)+1])  
##   info2 <- lapply(info1, function(x)x[match(INFOID, x)+1])
##   INFO <- matrix(unlist(info2), ncol=length(INFOID), byrow=TRUE)
##   colnames(INFO) <- INFOID

  if("INDEL" %in% INFOID){
    info2[, "INDEL"] <- ifelse(is.na(info2[,"INDEL"]), FALSE, TRUE)
  }
  
  vcfdata$INFO <- info2
  #SAMPLE
  if(is.null(FORMATID)){
    Sh <- vcfhead[grep("##FORMAT", vcfhead)]
    FORMATID <- unlist(lapply(strsplit(Sh, split="=|,"), function(x)x[3]))
  }
  for(n in 1:(ncol(vcf)-9)){
    eval(parse(text=paste("SAMPLE <- SAMPLE",n, sep="")))
    score1 <- strsplit(vcfdata$SAMPLE, split=":")
    format1 <- strsplit(vcfdata$FORMAT, split=":")  
    score2 <- lapply(1:length(score1), function(x)score1[[x]][match(FORMATID, format1[[x]])])
    SAMPLE <- matrix(unlist(score2), ncol=length(FORMATID), byrow=TRUE)
    colnames(SAMPLE) <- FORMATID
    eval(parse(text=paste("vcfdata$SAMPLE",n," <- SAMPLE", sep="")))
  }
  class(vcfdata) <- "vcf"
  return(vcfdata)
}

print.vcf <- function(x, ...){
  cat("VCF data\n")
  cat(paste("Calls: ", length(x$POS), " postions, ", length(x)-10, " sample(s)\n", sep=""))
  cat(paste("Names: ", paste(names(x), collapse=" "), "\n", sep=""))
}

summary.vcf <- function(object, ...){
  cat(paste("CHROM: ", paste(head(unique(object$CHROM), n=3), collapse=","), "...", paste(tail(unique(object$CHROM), n=3), collapse=","), " (", length(unique(object$CHROM)), ")\n", sep=""))
  cat(paste("POS: ", length(object$POS), " calls\n", sep=""))
  cat(paste("ID: ", paste(head(object$ID, n=3), collapse=","), ",...\n", sep=""))
  cat(paste("REF: ", paste(head(object$REF, n=3), collapse=","), ",...\n", sep=""))
  cat(paste("ALT: ", paste(head(object$ALT, n=3), collapse=","), ",...\n", sep=""))
  cat(paste("QUAL:\n", sep=""))
  print(summary(object$QUAL))
  cat(paste("FILTER: ", paste(unique(object$FILTER), collapse=","), "\n", sep=""))
  cat(paste("INFO: ", paste(colnames(object$INFO), collapse=","), "\n", sep=""))
  cat(paste("SAMPLE: ", paste(colnames(object$SAMPLE), collapse=","), "\n", sep=""))  
}
