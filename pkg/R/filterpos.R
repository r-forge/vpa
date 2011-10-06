
#VCF, bed, tbi for dbSNP, 1KG, custom vcf
filterpos <- function(vcf, position=NULL, file="", type="vcf", tbi=FALSE, chr=TRUE, tabix="tabix", ...){
  if(chr==TRUE){
    if(is.na(pmatch("chr", vcf$CHROM))){
      CHR <- paste("chr", vcf$CHROM, sep="")
    }else{
      CHR <- vcf$CHROM
    }
  }else{
    if(!is.na(pmatch("chr", vcf$CHROM))){
      CHR <- sub("chr", "", vcf$CHROM)
    }else{
      CHR <- vcf$CHROM
    }
  }
  Pos <- cbind(CHR, vcf$POS)

  if(!is.null(position)){
    posf <- position
  }else{    
    if(tbi==TRUE){
      seqs <- pos2seq(Pos, Seqfile=file, tabix=tabix)$vcf
    }else{
      seqs <- try(read.table(file, ...), TRUE)
      if(inherits(seqs, "try-error")){
        stop("Wrong file format or arguments for read.table")
      }
    }
    if(type=="bed"){
      posf <- cbind(seqs[, 1], as.numeric(seqs[, 2])+1, seqs[, 3])
    }else if(type=="vcf"){
      posf <- seqs[, c(1,2,2)]
    }else if(type=="gff"){
      posf <- seqs[, c(1,4,5)]
    }else{
      stop("Invalid file type")
    }    
  }

  index <- c()
  for(i in 1:nrow(posf)){
    index1 <- Pos[,1]==posf[i,1]
    index2 <- as.numeric(Pos[,2])>= as.numeric(posf[i,2]) & as.numeric(Pos[,2])<= as.numeric(posf[i,3])
    index <- c(index, which(index1 & index2))
  }
  filter <- rep(FALSE, nrow(Pos))
  filter[index] <- TRUE
  
  #filter <- paste(Pos[,1], Pos[,2]) %in% paste(posf[,1], posf[,2])

  filteredvcf <- subvcf(vcf, TF=!filter)
  droppedvcf <- subvcf(vcf, TF=filter)
  return(list(filtered=filteredvcf, dropped=droppedvcf))
}
