\name{getFactorCodes}
\alias{getFactorCodes}
\alias{getFactorCodes,Gmp-method}
\title{Get Factor Codes from Gmp Object}
\usage{
getFactorCodes(object)
}
\description{
  Retrieves the factor codes for independent variables from a Gmp
  object.  The underlying coding scheme for GMPMs is deviation coding.
}
\arguments{
  \item{object}{A Gmp object, created by \link{gmpCreate} or \link{gmp}.}
}

\value{
  A list, with named elements corresponding to the independent
  variables, and a matrix containing the codes for that variable.
}

\author{Dale J. Barr <dale.barr@ucr.edu>}
\seealso{\code{\link{gmp}}}
