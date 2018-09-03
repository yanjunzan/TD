\name{Extract_all}
\alias{Extract_all}
\title{
Extract the imputed genotypes  and  compile it into a recombination block genotype matrix.

}
\description{
In this matrix, each row represents one individual and each column a block in the genome in which no crossover event was identified in the population.
}
\usage{
Extract_all(chromosome = chr,id_all = id_all,all_vcf = all_vcf,gap=3e6,filter = T )
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
return a genotype matrix, with the first two rows representing chr and mid-pos the bins. Starting from the third row, each row represents one individual and each column is a block in the genome in which no crossover event was identified in the population.
%%  ~Describe the value returned
%%  If it is a LIST, use
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
\code{\link{Extract_co_all_chr}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
