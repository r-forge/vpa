\name{getref}
\alias{getref}
\title{Download reference gene table}
\description{To download reference gene table from UCSC golden path.
}
\usage{
getref(ref="hg19")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{The version of reference genome, such as hg18, hg19, etc.}
%%     ~~Describe \code{x} here~~
}
\details{
  List of reference genomes can be found in \url{http://hgdownload.cse.ucsc.edu/goldenPath/}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A data.frame from flat table of reference genome.
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
  \code{\link{Pos2Gene}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#gereq(varRes1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ reference genome }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
