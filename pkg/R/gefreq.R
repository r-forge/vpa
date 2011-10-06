gefreq <- function(vcf, method="fisher.test", level="gene", ref="hg19", ...){
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

  #Pos2Gene
  PosAnn <- c()
  for(i in 1:length(Pos)){
    chrpos <- unlist(strsplit(Pos[i], split=":"))
    PosAnn <- rbind(PosAnn, Pos2Gene(chr=chrpos[1], pos=chrpos[2], level=level, ref=ref))
  }
  PosAnn <- cbind(Pos, PosAnn)
  if(length(grep("/", PosAnn[,3]))>0){ #multiple genes for one variant
    gann <- strsplit(PosAnn[,3], split="/")
    len <- unlist(lapply(gann, length))
    PosAnn <- cbind(rep(PosAnn[,1], len), rep(PosAnn[,2], len), unlist(gann))
  }
  #counts
  GenesA <- paste(PosAnn[,2], PosAnn[,3], sep=":")
  Genes <- unique(GenesA)
  
  altm <- matrix("0", nrow=length(Genes), ncol=length(pos))
  for(i in 1:length(pos)){
    genei <- GenesA[PosAnn[,1] %in% pos[[i]]]
    genec <- table(genei)
    altm[match(names(genec), Genes), i] <- genec
  }
  colnames(altm) <- vcf$Samples[,1]
  rownames(altm) <- Genes
  #frequency
  frequency <- t(apply(altm, 1, function(x)tapply(x, groups, function(y)sum(y>0)/length(y))))
  #test
  p.value <- c()
  for(i in 1:length(Genes)){
    atable <- rbind(
                    tapply(altm[i,], groups, function(x)sum(x>0)),
                    tapply(altm[i,], groups, function(x)sum(x==0))
                    )
    if(!is.na(pmatch(method, "fisher.test"))){
      p.value[i] <- fisher.test(atable, ...)$p.value
    }else if(!is.na(pmatch(method, "chisq.test"))){
      p.value[i] <- chisq.test(atable, ...)$p.value
    }else{
      stop("invalid method")
    }
  }
  freq <- cbind(altm, frequency, p.value)
  freq <- freq[order(Genes,p.value),]
  #annotation
  colnames(PosAnn) <- c("position", "annotation", "genename")
  ann <- lapply(pos, function(x)PosAnn[match(x, PosAnn[,1]),])
  
  list(frequency=freq, annotation=ann)  
}
