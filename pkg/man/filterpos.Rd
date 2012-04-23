\name{filterpos}
\alias{filterpos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Filter variants against known SNP dataset}
\description{
  The function is used to filter variants against known SNP dataset in VCF, bed,
  gff or user-specified position files. For example, variants in VCF format
  can be filtered against dbSNP, 1000 genome project dataset, customized VCF data and so on.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filterpos(vcf, position=NULL, file="", type="vcf", tbi=FALSE, chr=TRUE,
tabix="tabix", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf}{A VCF object for filtering.
    %%     ~~Describe \code{x} here~~
  }
  \item{position}{A data.frame or matrix with chromosome names in the
    first column, start positions in the second column and end positions
    in the third column (1-based). This can be used to filter against
    customized VCF data.
  }
  \item{file}{The file containing the known SNPs.}
  \item{type}{The date format of input file. It can be 'vcf', 'bed' or 'gff'.}
  \item{tbi}{Logical value. If TRUE, the input file should be indexed by
    tabix for efficient information retrieval.}
  \item{chr}{Logical value. If TRUE, the chromosome names of the input
    file should have the prefix of 'chr', e.g. 'chr1'. If FALSE, the
    chromosome names don't have the 'chr' prefix.}
  \item{tabix}{The path of tabix function. if NULL, scanTabix function
    from Rsamtools will be used instead.}
  \item{...}{More arguments for read.table when reading the input file.}
}
\details{
  Variants can be filtered against dbSNP and 1000 genome project dataset
  to eliminate common variants.
  
  For example, the dbSNP 132 can be download from UCSC
  (\url{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp132.txt.gz}). The
  2nd-5th columns of the dataset can be extracted easily using 'cut' or
  'awk' into a bed format file. The bed file can be indexed by
  'tabix' for efficient information retrieval. The filterpos function can be used to 
   eliminate the variants observed in the indexed dataset file, with
  arguments type="bed" and tbi=TRUE.
  
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The input vcf data will be filtered against known SNP database or user-specified
  position files. A list including filtered vcf data and dropped vcf data will be returned.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Qiang Hu
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{pos2seq}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# ivcffile1 <- system.file("extdata", "1151HZ0001.flt.vcf", package="VPA")
# vcfdata1 <- read.vcf(vcffile1)
# vcffile2 <- system.file("extdata", "1151HZ0006.flt.vcf", package="VPA")
# vcfdata2 <- read.vcf(vcffile2)
# vcf <- filterpos(vcfdata1, position=cbind(vcfdata2$CHROM, vcfdata2$POS,
#  vcfdata2$POS), chr=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ filter }
