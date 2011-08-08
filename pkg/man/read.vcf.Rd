\name{read.vcf}
\alias{read.vcf.default}
\alias{read.vcf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read VCF file}
\description{
  Load VCF format file of variant and/or sequence calls into a flexible VCF object in R environment.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.vcf(file, VCF=NULL, INFOID = NULL, FORMATID = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{VCF format file of sequence calls from tools such as
    samtools and GATK.
%%     ~~Describe \code{file} here~~
}
\item{VCF}{An object list from the function \code{pos2seq} output. }
\item{INFOID}{Characters. Only specified elements in INFO column of
    VCF file will be read.
%%     ~~Describe \code{INFOID} here~~
}
  \item{FORMATID}{Characters. Only specified elements in FORMAT column
    of VCF file will be read.
%%     ~~Describe \code{FORMATID} here~~
}
  \item{\dots}{not used.
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A vcf object. Also it is a list. Each element of the list is from a column of the
  VCF file. See HEAD for details.
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
  \code{\link{write.vcf}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##read example vcf file in data directory.
vcffile <- system.file("extdata", "1151HZ0001.flt.vcf", package="VPA")
vcfdata <- read.vcf(vcffile)
summary(vcfdata)
write.vcf(vcfdata)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ vcf }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
