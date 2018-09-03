\name{Extract_co_all_chr}
\alias{Extract_co_all_chr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the number of  genotype swith ,which is supposed to approximate  cross over, more accurated double cross over estimation should be done in Rqtl.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Extract_co_all_chr(chromosome = chr,id_all = id_all,all_vcf = all_vcf,gap=3e6,filter=T)
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
  \item{co }{A list contain number of double genotype swich from chromsome 1 to N, and their sum for all the individuals}

%%  \item{comp1 }{Description of 'comp1'}
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
