\name{filtervcf}
\alias{filtervcf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variant filtering}
\description{
  VCF format file contains various score to assess the positional-level quality of
  variant and sequence call. Function filtervcf can be used to filter
  variants with user-specified quality criteria.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filtervcf(vcf, alter = NULL, alter.PL=20, alter.AD=3, alter.ADP=NULL,
QUAL = 20, DP = c(10,500), GQ = NULL, FILTER = NULL, INDEL = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf}{A VCF object for filtering.
%%     ~~Describe \code{x} here~~
}
  \item{alter}{Logical value. If TRUE, the variant positions are kept. If FALSE, the variant positions are discarded. If NULL, the option will be ignored.
%%     ~~Describe \code{alter} here~~
}
\item{alter.PL}{Phred-scaled genotype likelihoods of variant call to define a variant.}
\item{alter.AD}{The minimum depth of alter allele when alter is TRUE. The
  alt allel depths are read from AD (GATK) or DP4 (Samtools) in the VCF
  data.}
\item{alter.ADP}{The minimum depth pecentage of alter allele.}
  \item{QUAL}{phred-scaled variant likelihoods of variant call.
%%     ~~Describe \code{QUAL} here~~
}
  \item{DP}{position-level read depth.
%%     ~~Describe \code{DP} here~~
}
  \item{GQ}{Phred-scaled score for most likely genotype at position of interest.}
%%     ~~Describe \code{GQ} here~~
  \item{FILTER}{Quality label to keep. The VCF format of variant call by GATK
    VariantFiltration will label quality of each position.
%%     ~~Describe \code{FILTER} here~~
}
  \item{INDEL}{Logical value. If TRUE, only INDELs are evaluated. If FALSE, only point variants are evaluated.
%%     ~~Describe \code{INDEL} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The input vcf data will be filtered by user-specified quality
  criteria. A list including filtered vcf data and dropped vcf data will
  return.
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
  \code{\link{subvcf}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#Filter alignment artifacts to get promising SNP
vcffile <- system.file("extdata", "1151HZ0001.flt.vcf", package="VPA")
vcfdata <- read.vcf(vcffile)
vcflt <- filtervcf(vcfdata, alter=TRUE, alter.AD=3, QUAL=20, DP=c(10,500), GQ=20, INDEL=FALSE)$filtered
write.vcf(vcflt)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ filter }
