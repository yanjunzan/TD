\name{Extract_Double_co_all_chr}
\alias{Extract_Double_co_all_chr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the number of double genotype swith ,which is supposed to approximate double cross over, more accurated double cross over estimation should be done in Rqtl.
}
\description{
more accurated double cross over estimation should be done in Rqtl.
}
\usage{
 Extract_Double_co_all_chr(id_all,all_vcf,chromosome = chr,gap=3e6,filter=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id_all}{the Id of all the individuals after TIGER run}
  \item{filter}{ if the double recombination within a certain distance(gap) should be filtered }
  \item{all_vcf}{links to the TIGER output}
  \item{chromosome}{all the chromsomes in TIGER run}
  \item{gap}{ genotype swithing twich within this distance will be filtered out}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{co }{A list contain number of double genotype swich from chromsome 1 to N, and their sum for all the individuals}
  \item{dis }{the distances between those double switches for each individual}
%% ...
}
\references{
Yanjun Zan, Thibaut Payen, Mette Lillie Leif Andersson, Paul B. Siegel and Örjan Carlborg;Whole-genome QTL mapping in experimental pedigrees from outbred founders utilizing low coverage individual based sequencing 2018
}
\author{
Yanjun Zan, Thibaut Payen,
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{}

% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
