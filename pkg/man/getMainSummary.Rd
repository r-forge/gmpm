\name{getMainSummary}
\alias{getMainSummary}
\title{Get Main Summary of Results for GMPM object}
\description{
Provides a main summary of p-values for a \code{Gmp} object.
}
\usage{
getMainSummary(x, byCovar = FALSE)
}
\arguments{
  \item{x}{A \code{Gmp} object}

  \item{byCovar}{A boolean (TRUE, FALSE) value indicating whether
    results are to be given for each individual covariate separately, or
    for all covariates combined.  Default is FALSE.}

}
\details{
This function will intelligently determine whether a single-parameter or
multi-parameter test is needed for each independent variable in the
experiment, and perform a corresponding test.
}
\value{
  A list of data.frame objects.  In the normal case, a list with a
  single element with the main results.  In the multinomial case, a list
  with as many elements as there are comparisons among DV levels.
}
\seealso{ \code{\link{getRegSummary}} }
