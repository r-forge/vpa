\name{filterpos}
\alias{filterpos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Position filtering}
\description{
  The function is used to filter variants with positions from VCF, bed,
  gff or user-specified positions. For example, variants in VCF format
  can be filtered by this function with positions from dbSNP, 1000
  genome, customized VCF data and so on.
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
    in the third column (1-based).
  }
  \item{file}{The file is used to filtered the VCF object.}
  \item{type}{The date format of input file. It can be 'vcf', 'bed' or 'gff'.}
  \item{tbi}{Logical value. If TRUE, the input file should be indexed by
    tabix.}
  \item{chr}{Logical value. If TRUE, the chromosome names of the input
    file should have the prefix of 'chr', e.g. 'chr1'. If FALSE, the
    chromosome names don't have the 'chr' prefix.}
  \item{tabix}{The path of tabix function. if NULL, scanTabix function
    from Rsamtools will be used instead.}
  \item{...}{More arguments for read.table when reading the input file.}
}
\details{
  Variants can be filtered by dbSNP and 1000 genome to find novel
  mutations.
  
  For example, the dbSNP 132 can be download from UCSC
  (\url{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp132.txt.gz}). The
  2-5 columns of the data set can be extracted easily by 'cut' or
  'awk' to be a bed format file. The bed file can be indexed by
  'tabix'. Then the indexed file can be used to filter the variants
  using the function with arguments type="bed" and tbi=TRUE.
  
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The input vcf data will be filtered by SNP databases or user-specified
  positions. A list including filtered vcf data and dropped vcf data
  will return.
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
vcffile1 <- system.file("extdata", "1151HZ0001.flt.vcf", package="VPA")
vcfdata1 <- read.vcf(vcffile1)
vcffile2 <- system.file("extdata", "1151HZ0006.flt.vcf", package="VPA")
vcfdata2 <- read.vcf(vcffile2)
vcf <- filterpos(vcfdata1, position=cbind(vcfdata2$CHROM, vcfdata2$POS, vcfdata2$POS), chr=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ filter }
