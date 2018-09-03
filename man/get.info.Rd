\name{get.info}
\alias{get.info}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to generate information on how many windows are there on each chromsome/contig and theire strt and end pos in a linearised index
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.info(chroms,chroms.len,bin.size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chroms}{ a vector of chromsome/contig names}
  \item{chroms.len}{ a vector of chromsome/contig length}
  \item{chroms.len}{ a vector of chromsome/contig length}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Return a list with
 \item{num.bin }{A vector with number of bins on each contig}
 \item{index }{A data.frame with start and end of each contig in a linearised index}
 \item{loca }{A vector of mid position of each bin}
 \item{loca.chr }{A vector of names for each bin}
 \item{chr.loca}{ a contanated version of chr-pos, used for providing a informative name for each bin}
}
\references{
Yanjun Zan, Thibaut Payen, Leif Andersson, Paul B. Siegel and Örjan Carlborg;Whole-genome QTL mapping in experimental pedigrees from outbred founders utilizing low coverage individual based sequencing 2018
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
