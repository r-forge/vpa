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
filtervcf(x, alter = TRUE, alter.PL=20, QUAL = 20, DP = 10, GQ = NULL,
FILTER = NULL, INDEL = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A VCF object for filtering.
%%     ~~Describe \code{x} here~~
}
  \item{alter}{Logical value. If TRUE, the positions of interest are kept. If FALSE, the positions of interest are discarded.
%%     ~~Describe \code{alter} here~~
}
\item{alter.PL}{Phred-scaled genotype likelihoods of variant call.}
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
  The input and output variant data are both VCF object.
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
  Qiang Hu%, Dan Wang, Li Yan, Song Liu
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
vcflt <- filtervcf(vcfdata, alter=TRUE, QUAL=20, DP=10, GQ=20, INDEL=FALSE)
write.vcf(vcflt)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ filter }
