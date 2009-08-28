\name{mdTest}
\alias{mdTest}
\alias{mdTest,Gmp,vector-method}
\alias{mdTest,Gmp.mul,vector-method}
\alias{mdTest,matrix-method}
\alias{mdTest,matrix,ANY-method}
\title{Mahalanobis Distance Test for GMPMs}
\description{
Performs multiparameter tests on GMPM objects (or matrices) using
Mahalanobis distances calculated on the permutation matrix.
}
\usage{
mdTest(x, y)
}
\arguments{
  \item{x}{ Either a \code{Gmp} object or a permutation matrix}
  \item{y}{ Either a list of vectors, or a single vector.  The vector
    (or vectors) should index the columns of the matrix on which to
    perform the test.  They can either be numeric indices or the names
    of the matrix columns.  If \code{y} is a list, then a separate test
    will be performed for each vector in the list. }
}
\value{
  A \code{\linkS4class{Mdtest}} object
}
\seealso{ S4 classes \code{\linkS4class{Mdtest}},
  \code{\linkS4class{Gmp}}, and function \code{\link{gmp}} }
\examples{
# create a within-subjects design
# we are interested in main effect of A (3-level variable)
# this requires a multivariate test
df1 <- data.frame(SubjID=rep(1:12,each=6),
                  A=rep(c("Aa","Ab","Ac"),each=2),
                  B=c("Ba","Bb"))
df1$Y <- rep(c(0,1,-1),each=2) + rnorm(72)

# parametric analysis
summary(aov(Y ~ A*B + Error(SubjID), df1))

# randomization analysis
df1.gmp <- gmp(Y ~ A*B | SubjID, gaussian, df1, c("A","B"))

mdTest(df1.gmp, c("A1","A2"))
mdTest(df1.gmp, c(2,3))

# demonstrate "list" usage
mdTest(df1.gmp, list(omnibus=c("A1","A2","B1"),mainEffA=c("A1","A2")))

# retrieve the permutation matrix
pmx <- getPermMx(df1.gmp)

# perform test using permutation matrix instead
mdTest(pmx, c("A1","A2"))

}
\keyword{multivariate}
\keyword{htest}
\keyword{robust}
