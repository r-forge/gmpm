\name{permSpace}
\alias{permSpace}
\title{Calculate the Number of Possible Permutations}
\description{
Calculates the number of possible permutation based on an analysis of
the design.
}
\usage{
permSpace(object)
}
\arguments{
  \item{object}{A \code{Gmp} object}
}
\value{
  A number, possibly +Inf the number of permutations is larger than the
  largest number that can be represented in R.
}
