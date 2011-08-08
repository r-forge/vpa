\name{subvcf}
\alias{subvcf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Subset of VCF data}
\description{
  To obtain a subset of a vcf object data.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subvcf(x, CHR = NULL, POS = NULL, CHRPOS=NULL, samples = NULL, TF =
NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vcf object data.
    %%     ~~Describe \code{x} here~~
  }
  \item{CHR}{Chromosome. To get the subset of input vcf data based on
    specified chromosome(s).
    %%     ~~Describe \code{CHR} here~~
  }
  \item{POS}{Position. To get the subset of input vcf data based on
    specified positions.
    %%     ~~Describe \code{POS} here~~
  }
  \item{CHRPOS}{Position within a chromosome seperated by colon.}
  \item{samples}{To specify sample(s) of interest.
    %%     ~~Describe \code{samples} here~~
  }
  \item{TF}{A vector of logical values. To define which of the
    corresponding positions will be extracted.
    %%     ~~Describe \code{TF} here~~
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A vcf object returned.
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
  Qiang Hu%, Dan Wang, Li Yan, Song Liu.
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{filtervcf}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

vcffile <- system.file("extdata", "1151HZ0001.flt.vcf", package="VPA")
vcfdata <- read.vcf(vcffile)

#extract calls in position 1:985999
subvcf(vcfdata, CHRPOS="1:985999")

#extract calls by TF
tf <- c(rep(TRUE, 10), rep(FALSE, length(vcfdata$POS)-10))
subvcf(vcfdata, TF=tf)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ subset }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
