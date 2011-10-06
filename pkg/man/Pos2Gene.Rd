\name{Pos2Gene}
\alias{Pos2Gene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write VCF object}
\description{The function is used to map a location in the genome to its
  gene annotation in transcript or exon level. The annotation
  information is based on the refseq table in the UCSC genome.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Pos2Gene(chr, pos, level="gene", show.dist=FALSE, ref="hg19")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{Chromosome name, such as 'chr1' or '1'.
    %%     ~~Describe \code{x} here~~
  }
  \item{pos}{A location in the reference genome.}
  \item{level}{Annotation level for a position. It can be 'gene' or
    'exon' so far.}
  \item{show.dist}{Logical value. When a position is mapped as
    inter-gene, wether to show the distance between the two genes.
    %%     ~~Describe \code{file} here~~
  }
  \item{ref}{The version of reference genome.}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  A vector including annotation and gene name. If the level is 'gene', a
  position will be annotated as 'gene' or 'intergene'. If it is 'exon', a
  position will be annotated as 'exon', 'intron' or 'intergene'.

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
