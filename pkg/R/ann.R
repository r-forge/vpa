
getref <- function(ref="hg19"){
  tmpfile <- tempfile()
  download.file(paste("http://hgdownload.cse.ucsc.edu/goldenPath/", 
                      ref, "/database/refFlat.txt.gz", sep = ""), tmpfile, mode = "wb")
  refGene <- read.delim(conn <- gzfile(tmpfile, open = "rt"), header = FALSE, sep = "\t", stringsAsFactors=FALSE)
  close(conn)
  download.file(paste("http://hgdownload.cse.ucsc.edu/goldenPath/", 
                      ref, "/database/refFlat.sql", sep = ""), tmpfile, mode = "wb")
  refName <- readLines(tmpfile)
  ni <- grep("CREATE|KEY", refName)
  rfn <- c()
  for(i in (ni[1]+1):(ni[2]-1)){
    rfn <- c(rfn, strsplit(refName[i], split="`")[[1]][2])
  }
  colnames(refGene) <- rfn
  return(refGene)
}

Pos2Gene <- function(chr, pos, level="gene", show.dist=FALSE, ref="hg19"){
  if (!exists("refflat")) {
    reftmp <- list()
    reftmp[[ref]] <- getref(ref)
    assign("refflat", reftmp, .GlobalEnv)
  }
  else if (!(ref %in% names(refflat))) {
    refflat[[ref]] <- getref(ref)
  }
  ann <- refflat[[ref]]
  
  if(!is.na(pmatch("chr", chr))){
    chr <- chr
  }else{
    chr <- paste("chr", chr, sep="")
  }
  pos <- as.integer(pos)
  ann1 <- ann[ann[, "chrom"]==chr,]  
  gi <- which(pos>=as.integer(ann1[, "txStart"])+1 & pos<=as.integer(ann1[, "txEnd"]))
  if(length(gi)>0){
    if(!is.na(pmatch(level, "gene"))){
      gene <- paste(unique(ann1[gi, "geneName"]), collapse="/")    
      Gene <- c("gene", gene)      
    }else if(!is.na(pmatch(level, "exonic"))){
      estart <- ann1[gi, "exonStarts"]
      estart <- strsplit(estart, split=",")
      eend <- ann1[gi, "exonEnds"]
      eend <- strsplit(eend, split=",")
      ann1e <- c()
      for(i in 1:length(gi)){
        ann1e <- rbind(ann1e, cbind(ann1[gi[i], "geneName"], estart[[i]], eend[[i]]))        
      }
      ei <- which(pos>=as.integer(ann1e[, 2])+1 & pos<=as.integer(ann1e[, 3]))
      if(length(ei)>0){
        gene <- paste(unique(ann1e[ei, 1]), collapse="/")    
        Gene <- c("exon", gene)
      }else{
        gene <- paste(unique(ann1[gi, "geneName"]), collapse="/")
        Gene <- c("intron", gene)
      }
    }else{
      stop("Invalid level")
    }    
  }else{  
    s1 <- as.integer(pos)-as.integer(ann1[, "txEnd"])
    m1 <- min(s1[s1>=0])
    g1 <- ann1[match(m1, s1), "geneName"]
    e1 <- as.integer(ann1[, "txStart"])+1-as.integer(pos)
    m2 <- min(e1[e1>=0])
    g2 <- ann1[match(m2, e1), "geneName"]
    if(show.dist){
      gene <- paste(g1, "+", m1, ";", m2, "-", g2,  sep="")
    }else{
      gene <- paste(unique(c(g1, g2)), collapse=";")
    }
    Gene <- c("intergene", gene)
  }
  return(Gene)
}
