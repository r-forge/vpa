\name{write.vcf}
\alias{write.vcf.default}
\alias{write.vcf}
\alias{print.vcf}
\alias{summary.vcf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write VCF object}
\description{To write a vcf object to VCF format file
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write.vcf(x, file = "", HEAD=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vcf object.
%%     ~~Describe \code{x} here~~
}
  \item{file}{Character. The file (including path) to which a vcf object will be written.
    %%     ~~Describe \code{file} here~~
  }
  \item{HEAD}{Logical value. If TRUE, head information will be output.}
  \item{...}{not used.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
  \code{\link{read.vcf}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
vcffile <- system.file("extdata", "1151HZ0001.flt.vcf", package="VPA")
vcfdata <- read.vcf(vcffile)
summary(vcfdata)
write.vcf(vcfdata)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ vcf }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
