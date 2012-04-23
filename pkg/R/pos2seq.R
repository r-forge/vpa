pos2seq <- function(Pos, Seqfile, file="", tabix="tabix", region=5000){  #Pos: 2 column sdata.frame
  if(is.null(tabix)){
    require(Rsamtools)
    tbx <- TabixFile(Seqfile)
    param <- GRanges(Pos[,1], IRanges(start=as.integer(Pos[,2])-1, end=as.integer(Pos[,2])))
    res <- scanTabix(tbx, param=param)
    seq1 <- matrix(unlist(strsplit(unlist(res), split="\t")), nrow=length(res), byrow=TRUE)
    header <- headerTabix(tbx)$header
  }else{
  #read header
    pos1 <- "1:0-0"
    con <- pipe(paste(tabix, "-h", Seqfile, pos1, sep=" "))
    header <- readLines(con)
    close(con)
    
    Pos2region <- function(pos){
      paste(as.character(pos[1]), ":", as.integer(pos[2]), "-", as.integer(pos[2]), sep="")
    }
    PosReg <- apply(Pos, 1, Pos2region)
  #PosReg <- sub("M", "MT", PosReg) #sub M to MT
  
  #region <- 5000
    regs <- seq(1, length(PosReg), region)
  
    P2Seq <- list()
    if(length(regs)==1){
      con <- pipe(paste(tabix, Seqfile, paste(PosReg, collapse=" "), sep=" "))
      if(length(scan(con, what="", n=1, quiet=TRUE))>0){
        P2Seq[[1]] <- read.table(con, sep="\t")
      }else{
        P2Seq[[1]] <- NULL
        close(con)
      }
    }else{
      for(i in 1:(length(regs)-1)){
        con <- pipe(paste(tabix, Seqfile, paste(PosReg[(regs[i]):(i*region)], collapse=" "), sep=" "))
        if(length(scan(con, what="", n=1, quiet=TRUE))>0){
          P2Seq[[i]] <- read.table(con, sep="\t")
        }else{
          P2Seq[[i]] <- NULL
        }
        #pos2seq <- read.table(con, sep="\t")
        #P2Seq[[i]] <- pos2seq
      #print(i)
      }
      con <- pipe(paste(tabix, Seqfile, paste(PosReg[max(regs):length(PosReg)], collapse=" "), sep=" "))
      if(length(scan(con, what="", n=1, quiet=TRUE))>0){
        pos2seq <- read.table(con, sep="\t")
      }else{
        pos2seq <- NULL
      }
      #pos2seq <- read.table(con, sep="\t")
      P2Seq[[length(regs)]] <- pos2seq
    }  
    seq1 <- c()
    if(length(P2Seq)>0){
    for(i in 1:length(P2Seq)){
        seq1 <- rbind(seq1, P2Seq[[i]])
      }
      #seq1 <- as.matrix(seq1)
    }else{
      seq1 <- matrix(nrow=0, ncol=10)
    }
  }
  seq1 <- unique(seq1)
  if(length(grep("^#CHROM", header))>0){
    colnames(seq1) <- strsplit(tail(header, n=1), split="\t|#")[[1]][-1]
  }
#  closeAllConnections()  
  if (file == ""){
    list(header=header, vcf=seq1)
  }else{
    write(header, file, sep="\t")
    write.table(seq1, file=file, append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
}
