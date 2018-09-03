\name{wrap_get_density}
\alias{wrap_get_density}

\title{
A function to get the number of markers along consecutive genomic bins from input files for TIGER
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wrap_get_density(inputfile,binsize,cut=T,chr.match,cutoff=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputfile}{ A txt input file for TIGER pipeline, usually strat with chr and pos followed by Alleles from one line, allele call,and alleles from another line,and corresponding allele call}
  \item{binsize}{The size of the genomic bins}
  \item{cut}{ If the bins with too few markers should be discarded}
  \item{chr.match}{ a file matching Contig names to chromsome,}
  \item{cutoff}{ if the number of markers in a bin is smaller than cutoff, this bin would be discarded}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
return a vector with number of markers for this input file.
%%  ~Describe the value returned
%%  If it is a LIST, use
%% \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Yanjun Zan, Thibaut Payen, Mette Lillie Leif Andersson, Paul B. Siegel and Örjan Carlborg;Whole-genome QTL mapping in experimental pedigrees from outbred founders utilizing low coverage individual based sequencing 2018
}
\author{
Yanjun Zan, Thibaut Payen
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
