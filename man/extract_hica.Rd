\name{extract_hica}
\alias{extract_hica}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extraction of score and loading matrices.
}
\description{
This function extracts the score matrix and the loading matrix given the dimension of the subspace considered and the level of the tree chosen. Furthermore it provides the cumulant energies for the subspace extracted.
}
\usage{
extract_hica(energy.obj, comp, level)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{energy.obj}{
An object provided by the function \code{\link{energy_hica}}.
}
  \item{comp}{
Dimension of the subspace.
}
  \item{level}{
Level of the tree.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{X}{data matrix.}
\item{S}{score data matrix.}
\item{C}{loading matrix. Each column represents a basis element.}
\item{cum.energy}{cumulant energy for the subspace extracted.}
}
\references{
P. Secchi, S. Vantini, and P. Zanini (2014). Hierarchical Independent Component Analysis: a multi-resolution non-orthogonal data-driven basis. \emph{MOX-report 01/2014}, Politecnico di Milano.
}
\author{
Piercesare Secchi, Simone Vantini, and Paolo Zanini.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{basis_hica}}, \code{\link{similarity_hica}}, \code{\link{energy_hica}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
