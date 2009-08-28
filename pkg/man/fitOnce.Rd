\name{fitOnce}
\alias{fitOnce}
\alias{fitOnce,Gmp-method}
\alias{fitOnce,Gmp.mul-method}
\title{Call Underlying Regression Function on Original Data}

\description{
  \code{fitOnce} calls the underlying regression function for a GMPM
  object on the stored model frame and returns the results of the call.
}
\usage{
fitOnce(object)
}
\arguments{
  \item{object}{A \code{\linkS4class{Gmp}} object}
}
\value{

  The results of a call to the underlying function.  This depends on the
  type of Gmp object (\code{\link{multinom}} for Gmp.mul objects,
  \code{\link{glm}} otherwise).

}
\seealso{ \code{\link{glm}}, \code{\link{multinom}}, \code{\link{gmp}} }
