\name{VPA-package}
\alias{VPA-package}
\alias{VPA}
\docType{package}
\title{Extract variants from VCF data with specified variant frequency pattern}
\description{
  VPA (Variant Pattern Analyzer) is a package for prioritizing variants with user-specified frequency
  pattern from multiple study subjects in next-generation sequencing study. The package starts from individual 
  files of sequence and   variant calls and the extract variants with user-specified frequency
  pattern across the study subjects of interest.  The frequency pattern can be analyzed at both variant level and gene level, 
  and functions are provided to assess the statistical significance of observed frequency difference. 
}
\details{
\tabular{ll}{
Package: \tab VPA\cr
Type: \tab Package\cr
Version: \tab 0.3.4\cr
Date: \tab 2011-11-04\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

}
\author{
Qiang Hu%, Dan Wang, Li Yan, Song Liu

Maintainer: Qiang Hu <huqmail@gmail.com>
}
\references{
}
\keyword{ package }
\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{


#setwd(system.file("extdata", package="VPA"))
#varflt <- LoadFiltering(file="index1.txt", filtering=TRUE)
#pattern <- cbind(A=c(1/4,1), B=c(0,0))
#varRes1 <- Patterning(varflt, pattern)


}
