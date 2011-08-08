\name{Patterning}
\alias{Patterning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract variants with user-specified variant pattern}
\description{
  To prioritize variants in user-specified frequency pattern. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Patterning(x, pattern, paired=FALSE, not.covered=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A \code{varlist} class data from the function \code{LoadFiltering}.
    %%     ~~Describe \code{x} here~~
  }
  \item{pattern}{The variant frequency matrix. Each column of the matrix
    is defined as the minimum and maximum value of variant frequency
    for each group of interest.}
  \item{paired}{Logical value. Whether cases and controls are paired. If
  paired is TRUE, control group label in index file should be marked as
  control. Sample names should be matched between case and its matched control.}
  \item{not.covered}{Logical value for the position with sequence coverage
    less than specified depth. If TURE, such low-coverage positions will be
    considered as variant. If FALSE, such low-coverage positions will be
    considered as reference. If it's NULL (default), such low-coverage positions
    will be filtered.}
}
\details{
   This function is used to extract variant with user-specified frequency pattern across study subjects. 
   The pattern matrix is specified by users in advance. The column names should be matched with sample group names.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The value returned is a \code{varlist}, including \code{VarVCF}, \code{VarFrequency}, 
  \code{Pattern} and \code{Samples}.
  \item{VarVCF}{A list variants with user-specified frequency pattern in each sample.}
  \item{VarFrequency}{Variant frequencies for input positions.}
  \item{Sample}{Samples annotation from input file.}
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

#setwd(system.file("extdata", package="VPA"))
#varflt <- LoadFiltering(file="index1.txt", filtering=TRUE)
#pattern <- cbind(A=c(1/4,1), B=c(0,0))
#varRes1 <- Patterning(varflt, pattern)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ variant }
\keyword{ pattern }% __ONLY ONE__ keyword per line
