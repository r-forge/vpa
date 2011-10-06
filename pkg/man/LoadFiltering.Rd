\name{LoadFiltering}
\alias{LoadFiltering.default}
\alias{LoadFiltering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{To load and filter variants in batch mode}
\description{
  To load data from study subjects and perform position-level quality filtering. 
  The index.txt file contains group status and VCF file location of each subject. 
  The function take index.txt   file as input to load variant and sequence call 
  files automatically. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LoadFiltering(file, datadir=NULL, filtering=TRUE, alter.PL=20,
alter.AD=3, alter.ADP=NULL, QUAL=20, DP=c(10,500), GQ=20, FILTER=NULL,
tabix="tabix", parallel=FALSE, pn=4, type=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
    Formatted input file including the annotation information of study
    subjects.
    %%     ~~Describe \code{x} here~~
  }
  \item{datadir}{The work directory of the index file and variants
    data. If it is NULL, the absolute path of variants files should be 
    provided in the annotation file.}
  \item{filtering}{Logical value. Whether to filter VCF data by
    specified quality criteria.}
  \item{alter.PL}{Phred-scaled genotype likelihoods to define a
    variant. If it's not NULL, PL scores from sequence-level calls are
    used to decide a variant besides from variant calls.}
  \item{alter.AD}{The minimum depth of alter allele.}
  \item{alter.ADP}{The minimum depth percentage of alter allele.}
  \item{QUAL}{Phred-like quality score cut-off to filter variants when
    \code{filtering} is TRUE.}
  \item{DP}{Sequencing depth cut-off to filter variants when
    \code{filtering} is TRUE.}
  \item{GQ}{Phred-like genotype score cut-off to filter variants when
    \code{filtering} is TRUE.}
  \item{FILTER}{Filter label when \code{filtering} is
    TRUE.}
  \item{tabix}{The file path of executable tabix.}
  \item{parallel}{If TRUE, the function will run in parallel model.}
  \item{pn}{The CPU numbers to be used if parallel is TRUE.}
  \item{type}{MPI type. See detail in \code{help(sfInit)}}
  \item{\dots}{
    Arguments to pass to the method \code{sfInit} of the snowfall package.
    %%     ~~Describe \code{\dots} here~~
  }
}
\details{
  \code{file} {The input file contains the annotation
    information of each sample. Each row is for one sample. The four
    columns are separated by tab, including sample name (required), group status (required),
    variant call file name (required) and sequence call file name
    (optional). Sample name column lists the sample name. Group status
    column lists the status (e.g., aggressive, benign or normal) of group 
    each sample belongs to. Variant call file name column lists the path
    of VCF formatted variant call file. Sequence call file name column
    lists the path of compressed VCF sequence call file. The high-volume
    data in tab-delimited VCF formats can be efficiently compressed by
    bgzip program and retrieved through tabix program from open-source
    SAMtool package. If the VCF format file is compressed by bgzip,
    tabix should be installed. The path of tabix should be specified in
    the function if it is not in the PATH system environment.}

  \code{Quality criteria} {The detail of quality scores in VCF data can
    be found at \url{http://www.1000genomes.org/node/101}.}
  
  \code{parallel} {This function will extract calls in sequential mode. If parallel is true, the function will extract calls in parallel mode. The package
    \code{Rmpi} and \code{snowfall} are required for parallel mode.}
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The value returned is a varlist, including \code{vcflist},
  \code{VarIndex} and \code{Samples}.
  \item{varlist}{A list of vcf objects, one for each sample. If the filtering is true, the variant data are filtered by
    specified quality criteria.}
  \item{VarIndex}{The indexes for all variant positions. TRUE
    denotes the presence of variant. FALSE denotes the absence of variant. NA denotes low coverage.}
  \item{Sample}{Samples annotation from the input index file.}
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

setwd(system.file("extdata", package="VPA"))
varflt <- LoadFiltering(file="index1.txt", filtering=TRUE, alter.PL=20, alter.AD=3)
pattern <- cbind(A=c(1/4,1), B=c(0,0))
varRes1 <- Patterning(varflt, pattern, var.PL=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ variant }% __ONLY ONE__ keyword per line
