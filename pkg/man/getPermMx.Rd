\name{getPermMx}
\alias{getPermMx}
\title{Get the Permutation Matrix from a GMPM object}
\description{
Retrieves the permutation matrix from a Gmp object.
}
\usage{
getPermMx(x)
}
\arguments{
  \item{x}{A \code{Gmp} object}
}

\value{
  A matrix object.  The first row of the permutation matrix includes the parameter
  estimates from the original fit.  Each remaining row has estimates
  from individual Monte Carlo runs.
}
