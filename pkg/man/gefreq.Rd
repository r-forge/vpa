\name{gefreq}
\alias{gefreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gene frequency of a varlist}
\description{
  To summarize the frequency of variants in gene level and estimate the
  statistical significance of frequency difference.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gefreq(vcf, method="fisher.test", p=1, level="gene", ref="hg19", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf}{A varlist object.}
  \item{method}{The test method for mutated genes. This must
    be one of "fisher.test" and "chisq.test".}
  \item{p}{The maximum of the p values.}
  \item{level}{The annotation level for variants. It can be either
    "gene" (i.e., including intron region) or "exon" only (i.e., without
    including intron region).}
  \item{ref}{The version of reference genome, e.g. "hg19".}
  \item{...}{More arguments for the test method.}
  %%     ~~Describe \code{file} here~~
}

\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  A list contains a dataframe of frequencies between groups and an
  annotation list of each samples.
  \item{frequency }{A data frame to list the gene name, variation
    distribution, variation frequency and p.value of all genes with
    variants across study groups.}
  \item{otherfreq }{A data frame to list the frequency results of the
    genes that are not in the specified level.}
  \item{annotation }{A list including gene annotations for the variants of
    each sample.}
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
  \code{\link{vcfreq}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#gefreq(varRes1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ frequency }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
