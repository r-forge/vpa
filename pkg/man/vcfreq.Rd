\name{vcfreq}
\alias{vcfreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variant frequency of a varlist}
\description{
  To summarize the frequency of extracted variants.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vcfreq(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A varlist object.
%%     ~~Describe \code{file} here~~
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A data frame to list the REF, ALT and frequecies of all
  positions across study groups.

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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##read example vcf file in data directory.
#vcfreq(varpat1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ frequency }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
