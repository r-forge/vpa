\name{vcfreq}
\alias{vcfreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variant frequency of a varlist}
\description{
  To summarize the frequency of variants and estimate the statistical
  significance of frequency difference.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vcfreq(vcf, method="fisher.test", p=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf}{A varlist object.}
  \item{method}{The allele frequency test method to be used. This must
    be one of "fisher.test" and "chisq.test".}
  \item{p}{The maximum of the p values.}
  \item{...}{More arguments for the test method.}
%%     ~~Describe \code{file} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A data frame to list the position, REF, genotypes, variant allele frequencies and
  p.value of all variants across study groups.

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
  \code{\link{gefreq}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##read example vcf file in data directory.
#vcfreq(varRes1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ frequency }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
