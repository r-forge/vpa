\name{Pos2Gene}
\alias{Pos2Gene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Map position to gene}
\description{The function is used to annoate a variant to its
  gene. The annotation information is based on the refseq table downloaded from the UCSC genome browser.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Pos2Gene(chr, pos, level="gene", show.dist=FALSE, ref="hg19")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{Chromosome name of the variant, such as 'chr1' or '1'.
    %%     ~~Describe \code{x} here~~
  }
  \item{pos}{The location of variant in the Chromosome.}
  \item{level}{The annotation level for variants. It can be either
    "gene" (i.e., including intron region) or 
    "exon" (i.e., without including intron region).
}
  \item{show.dist}{Logical value. When a position is mapped as
    inter-gene, whether to show the distances from the two genes.
    %%     ~~Describe \code{file} here~~
  }
  \item{ref}{The version of reference genome.}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  A vector including gene annotation results. 
  %If the level is 'gene', a position will be annotated as 'gene' or 'intergene'. 
  %If it is 'exon', a position will be annotated as 'exon', 'intron' or 'intergene'.

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
  \code{\link{getref}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Pos2Gene("1", "1000000", level="exon", show.dist=TRUE, ref="hg19")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ annotation }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
