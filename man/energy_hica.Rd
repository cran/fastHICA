\name{energy_hica}
\alias{energy_hica}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Energy criterion
}
\description{
This function implements the energy criterion defined in Secchi, Vantini, and Zanini (2013).
}
\usage{
energy_hica(HICA.obj, maxcomp = 1, nlevel = 1, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HICA.obj}{
An object provided by the function \code{\link{basis_hica}}.
}
  \item{maxcomp}{
The maximum space dimension considered.
}
  \item{nlevel}{
The number of levels analyzed. Specifically the levels from p-\code{nlevel} to p-1 are analyzed, where p is the number of variables.
}
  \item{plot}{
A logical value. If TRUE the energy is plotted.
}
}
\details{
This function computes the energy according the criterion presented in Secchi, Vantini and Zanini (2013). It is useful to find the best representation. It receives in input the output of the \code{\link{basis_hica}} function.
}
\value{
\item{energy}{A matrix with \code{maxcomp} rows and p-1 columns, where p is the number of variables. In position (i,j) it contains the energy of the best i-dimensional space for the jth level of the tree. Only the last \code{nlevel} columns are filled.}
\item{components}{A matrix with \code{maxcomp} rows and p-1 columns, where p is the number of variables. In position (i,j), it contains the index of th ith basis element for jth level of the tree. Only the last \code{nlevel} columns are filled.}
\item{HICA.obj}{The same object, output of the function \code{\link{basis_hica}}, provided in input.}
}
\references{
P. Secchi, S. Vantini, and P. Zanini (2014). Hierarchical Independent Component Analysis: a multi-resolution non-orthogonal data-driven basis. \emph{MOX-report 01/2014}, Politecnico di Milano.
}
\author{
Piercesare Secchi, Simone Vantini, and Paolo Zanini
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{basis_hica}}, \code{\link{similarity_hica}}, \code{\link{extract_hica}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
