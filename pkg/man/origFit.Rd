\name{origFit}
\alias{origFit}
\title{Get Original Regression Fit}
\description{
Get the regression fit to the data as originally labeled.
}
\usage{
origFit(object)
}
\arguments{
  \item{object}{A \code{Gmp} object}
}
\value{
Results of a call to the underlying regression function, typically a
call to \code{\link{glm}} or to \code{\link{multinom}} in the case of
multinomial data.
}
