\name{Export2rqtl}
\alias{Export2rqtl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function transform the output to Rqtl
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Export2rqtl(genoFile, phenoFile, output.name = paste0("/Users/yanjunzan/Documents/impute/results/F2_833.within.fam._Tiger_1Mb_bins_cut_", 5, "_", Sys.Date()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genoFile}{Genotype matrix output from this package}
  \item{phenoFile}{A datafram with 3 or more column as ID,sex,family id, phentpye1,...ID names have to match the ID names in genotype file}
  \item{output.name}{full path and file name for output file}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function write out a .csv file with name and path specified by output.name argument

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
