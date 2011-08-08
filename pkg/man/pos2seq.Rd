\name{pos2seq}
\alias{pos2seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Positions to sequencing calls}
\description{Function to retrieve variant or sequence call position of interest
  from tabix indexed files.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pos2seq(Pos, Seqfile, file="", tabix="tabix", region = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pos}{A two columns data.frame or matrix includes chromosome and
    position for each variant (or sequence) call.
    %%     ~~Describe \code{Pos} here~~
  }
  \item{Seqfile}{A tabix indexed VCF file include all variant (or sequence) calls to
    retrieve.
    %%     ~~Describe \code{Seqfile} here~~
  }
  \item{file}{File path to write out the retrieved results as a plain file in VCF format.}
  \item{tabix}{The path of tabix function. if NULL, scanTabix function from Rsamtools will be used instead.}
  \item{region}{The number of positions for tabix to retrieve at the
    same time. Too big number will not work for tabix. The default is
    5000.
    %%     ~~Describe \code{region} here~~
  }
}
\details{
  \code{pos2seq} requires tabix function from SAMtools. The path of tabix could be specified in
  the optional argument of function if it is not in the PATH system environment. 
  More details: \url{http://samtools.sourceforge.net/tabix.shtml}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A list includes header and VCF data. The results can also be outputted as a plain text file in VCF format.
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
  \code{\link{LoadFiltering}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

vcffile <- system.file("extdata", "1151HZ0001.flt.vcf", package="VPA")
vcfdata <- read.vcf(vcffile)

##extract calls from tabix indexed data
Pos <- cbind(vcfdata$CHROM, vcfdata$POS)
gzfile <- system.file("extdata", "1151HZ0006.vcf.gz", package="VPA")
calls <- pos2seq(Pos, gzfile)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tabix }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
