\name{basis_hica}
\alias{basis_hica}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construction of the HICA basis
}
\description{
This function builds the HICA tree up to a prespecified height providing the corresponding non-orthogonal bases.
}
\usage{
basis_hica(X, maxlev = dim(X)[2] - 1, dim.subset = 512)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data matrix with nrow(X) observations and ncol(X) variables.
}
  \item{maxlev}{
The maximum level of the tree. This must be an integer between 1 and ncol(X)-1. The default value is set to ncol(X)-1.
}
  \item{dim.subset}{
The dimension of the subset used for the evaluation of the similarity index (i.e., distance correlation). If this it is greater than nrow(X) all the observations are used, unless a random subsample of \code{dim.subset} observations is used. The default value is set to 512.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 \item{X}{
data matrix.
}
 \item{basis}{
a list with \code{maxlev} elements. The ith element of the list contains the basis matrix provided at level i of the tree. Each column of the basis matrix represent a basis element.
}
 \item{aggregation}{
a matrix with \code{maxlev} rows and 3 columns. At each row the first two columns contain the variable indeces merged at the corresponding level of the tree. In the third column the distance correlation of the two merged variables is recorded.}
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
\code{\link{energy_hica}}, \code{\link{similarity_hica}}, \code{\link{extract_hica}}
}
\examples{
## Not run:

##########################################################
# Example - Independent sources and overlapping loadings #
##########################################################

c1=c(0,0,0,0,1,1)
c2=c(1,1,1,1,0,0)
c3=c(1,1,0,0,0,0)

s1=runif(400,0,20)
s2=runif(400,0,20)
s3=runif(400,0,20)

# Here we generate the simulated dataset

X=s1\%*\%t(c1)+s2\%*\%t(c2)+s3\%*\%t(c3)+rnorm(6*400,0,1)

# Here we perform HICA on the simulated dataset

basis=basis_hica(X,5)

# Here we plot the 3 main components of HICA basis 
# (according to the energy criterium) for 4th level

energy=energy_hica(basis,6,5,plot=TRUE)
ex4=extract_hica(energy,3,4)
loa4=ex4$C

par( mfrow = c(3,1))
barplot(loa4[,1], ylim = c(-1, 1),main="HICA transform - Level 4",
ylab="1st component",xlab="Coordinate",names.arg=1:6,col="red",mgp=c(2.5,1,0))
barplot(loa4[,2], ylim = c(-1, 1),ylab="2nd component",
xlab="Coordinate",names.arg=1:6,col="green",mgp=c(2.5,1,0))
barplot(loa4[,3], ylim = c(-1, 1),ylab="3rd component",
xlab="Coordinate",names.arg=1:6,col="blue",mgp=c(2.5,1,0))

## End (Not run)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
