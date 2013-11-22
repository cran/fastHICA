\name{similarity_hica}
\alias{similarity_hica}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate of the similarity matrix
}
\description{
This function provides an estimate of the similarity matrix of the original data, before performing HICA algorithm.
}
\usage{
similarity_hica(X, dim.subset = 512)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data matrix with nrow(X) observations and ncol(X) variables.
}
  \item{dim.subset}{
The dimension of the subset used for the evaluation of the similarity index (i.e., distance correlation). If this it is greater than nrow(X) all the observations are used, unless a random subset of \code{dim.subset} observations is used. The default value is set to 512.
}
}
\details{
This function is auxiliary for the \code{\link{basis_hica}} function. Indeed its output is the estimate of the similarity matrix at the first step of the algorithm.
}
\value{
\item{similarity_matrix}{similarity matrix of the original data.}
\item{subset}{subset used for the evaluation of distance correlation between variables.}
}
\references{
P. Secchi, S. Vantini, and P. Zanini (2014). Hierarchical Independent Component Analysis: a multi-resolution non-orthogonal data-driven basis. \emph{MOX-report 01/2014}, Politecnico di Milano.
}
\author{
Piercesare Secchi, Simone Vantini, and Paolo Zanini.
}
\note{
The distance correlation is evaluated through the function \code{dcor} of the package "energy". It becomes computationally unfeasible if the number of observations is too large. For this reason it is possibile to choose the dimension of the subsample to be used in the evaluation of the similarity matrix. By default the dimension is set to 512.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{basis_hica}}, \code{\link{energy_hica}}, \code{\link{extract_hica}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
