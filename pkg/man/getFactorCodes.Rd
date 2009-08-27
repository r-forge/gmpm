\name{getFactorCodes}
\alias{getFactorCodes}
\title{Get Factor Codes from Gmp Object}
\description{}

\usage{
getFactorCodes(object)
}
\description{
  Retrieves the factor codes from a Gmp object.
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
