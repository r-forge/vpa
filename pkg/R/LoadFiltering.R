
LoadFiltering <- function(file, datadir=NULL, filtering=TRUE, alter.PL=20, alter.AD=3, alter.ADP=NULL, QUAL=20, DP=c(10,500), GQ=20, FILTER=NULL, tabix="tabix", parallel=FALSE, pn=4, type=NULL, ...)UseMethod("LoadFiltering")
LoadFiltering.default <- function(file, datadir=NULL, filtering=TRUE, alter.PL=20, alter.AD=3, alter.ADP=NULL, QUAL=20, DP=c(10,500), GQ=20, FILTER=NULL, tabix="tabix", parallel=FALSE, pn=4, type=NULL, ...){
  samples <- as.matrix(read.table(file, head=FALSE))
  if(!is.null(datadir)){
    vfilepath <- file.path(datadir, samples[,3])
  }else{
    vfilepath <- samples[,3]
  }
  groups <- samples[,2]
  #vcf
  vcfdata <- list()
  pos <- c()
  for(i in 1:nrow(samples)){
    cat(paste("Load variants for ", samples[i, 1], ":\n", sep=""))
    m1 <- read.vcf(file=vfilepath[i]) #read vcf
    pos1 <- cbind(m1$CHROM, m1$POS)
    if(filtering){ 
      m1 <- filtervcf(m1, alter=TRUE, alter.PL=alter.PL, alter.AD=alter.AD, alter.ADP=alter.ADP, QUAL=QUAL, DP=DP, GQ=GQ, FILTER=FILTER, INDEL=NULL)$filtered
    }
    vcfdata[[i]] <- m1
    pos <- unique(rbind(pos, pos1))
  }
  names(vcfdata) <- samples[,1]
  #load sequence file
  if(ncol(samples)>=4){
    if(length(grep("gz", samples[,4], ignore.case=TRUE))!=nrow(samples)){
      stop("The fourth column should be indexed VCF files of all positions.")
    }
    if(!is.null(datadir)){
      sfilepath <- file.path(datadir, samples[,4])
    }else{
      sfilepath <- samples[,4]
    }
    
    loadseq <- function(n, pos, samples, sfilepath, tabix){
      cat(paste("retrieve calls from ", samples[n, 1], " with tabix...\n", sep=""))
      VCF <- pos2seq(pos, sfilepath[n], tabix=tabix) #require pos2seq, tabix
      seqvcf <- read.vcf(VCF=VCF, INFOID="DP", FORMATID="PL")
      return(seqvcf)
    }
  
    cat("Extract sequence level data ...\n")
    seqdata <- list()
    if(parallel){
      library(snowfall)
      for(name in c("samples", "pos", "sfilepath", "loadseq", "tabix")){
        eval(parse(text=(paste("assign('",name,"', ", name, ", envir=.GlobalEnv)", sep=""))))
      }
      sfInit(parallel=TRUE, cpus=pn, type=type, ...)
      sfLibrary(VPA)
      sfExport(list=c("samples", "pos", "sfilepath", "loadseq", "tabix"))
      seqdata <- sfLapply(1:nrow(samples), function(n)loadseq(n=n, pos=pos, samples=samples, sfilepath=sfilepath, tabix=tabix))
      sfStop()
    }else{
      for(n in 1:nrow(samples)){ #controls of case1
        seqdata[[n]] <- loadseq(n=n, pos=pos, samples=samples, sfilepath=sfilepath, tabix=tabix)
      }
    }
    names(seqdata) <- samples[,1]
  }

  #INDEX
  POS <- paste(pos[,1], pos[,2], sep=":")
  vartf <- lapply(vcfdata, function(x)POS %in% paste(x$CHROM, x$POS, sep=":"))
  if(exists("seqdata")){
   # covtf <- list()
    for(n in 1:nrow(samples)){
      p1 <- paste(seqdata[[n]]$CHROM, seqdata[[n]]$POS, sep=":")
      p1tf <- ifelse(as.numeric(seqdata[[n]]$INFO[, "DP"])>=DP[1], TRUE, NA)
      if(!is.null(alter.PL)){ #possible variant
        p1pl <- filtervcf(seqdata[[n]], alter=TRUE, alter.PL=alter.PL, alter.AD=NULL, alter.ADP=NULL, QUAL=0, DP=DP)$filtered
        plp <- setdiff(paste(p1pl$CHROM, p1pl$POS, sep=":"), paste(vcfdata[[n]]$CHROM, vcfdata[[n]]$POS, sep=":"))
        vartf[[n]][match(plp, POS)] <- "PL"
      }
      p1tf <- unique(cbind(p1, p1tf))
      cov1tf <- p1tf[match(POS, p1tf[,1]), 2]
      vartf[[n]][is.na(cov1tf)] <- NA
    }
  }
  VarIndex <- matrix(unlist(vartf), ncol=nrow(samples))
  colnames(VarIndex) <- samples[,1]
  rownames(VarIndex) <- POS
    
  varbatch <- list(vcflist=vcfdata, VarIndex=VarIndex, Samples=samples)
  class(varbatch) <- "varlist"
  return(varbatch)
}


print.varlist <- function(x, ...){
  for(i in 1:length(x)){
    if(is.null(dim(x[[i]]))){
      print(x[i])
    }else{
      print(eval(parse(text=paste("list(",names(x[i]),"=head(x[[",i,"]], n=3))", sep=""))))
      cat("... ...\n")
    }
  }
}
