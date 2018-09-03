\name{Mat2geno}
\alias{Mat2geno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format the genotype matrix to consectutive bins along the genomes
}
\description{
This is a function takes the output from Extract_all and reforat the data.frame to a new data.frame with columns within the same genomic bins are avearged.
}
\usage{
Mat2geno(out = out_3,id = id_all[index.keep],length.chr = length.chr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{binsize}{ size of the consecutive bins}
  \item{out_3}{ output matrix from EXtract_all}
  \item{id_all}{all or a subset of the inviduals to keep in the final matrix}
  \item{length.chr}{a vector contain the length of each chromsome in Mb}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
return a matrix with each row as a individual  and each column as one consecutive genomic bin
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
