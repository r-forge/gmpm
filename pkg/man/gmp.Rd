\name{gmp}
\alias{gmp}
\alias{gmpCreate}
\alias{gmpFit}
\title{Creating and Fitting Generalized Multilevel Permutation Models (GMPMs)}
\description{
\code{gmpCreate} creates a GMPM object without performing permutation runs.
\code{gmpFit} performs permutation runs, and \code{gmp} is a convenience
function that performs a call to gmpCreate followed by a call to gmpFit.
}

\usage{

gmpCreate(formula, family, data=parent.frame(), ivars,
gmpControl=gmpCtrl(), ...)

gmpFit(object, gmpControl)

gmp(formula, family, data = parent.frame(), ivars, gmpControl = gmpCtrl(), ...)

}

\arguments{

  \item{formula}{A symbolic formula describing the regression model to
    be fitted, including a multilevel blocking factor following the '|'
    symbol.  Details of model specification can be found in the section
    \sQuote{Details}.}

  \item{family}{One of the following options: a valid glm class (e.g.,
    binomial, gaussian, poisson, etc.), \dQuote{multinomial}, and
    \dQuote{user} (for a user-defined model).  Options
    \dQuote{multinomial} and \dQuote{user} must be entered within
    quotation marks.}
    
  \item{data}{A data frame containing the data model will be fit to.  If
  left unspecified, gmp will search for the data within the parent
  environment.}

  \item{ivars}{Names of all independent variables (IVs) included as
  predictors in the model formula, entered as a vector of character
  strings.  All other predictors in the model will be considered
  \sQuote{covariates} and will not be subject to randomization.}

  \item{gmpControl}{A gmpControl object used to control the fitting
    function.  This is a list one or more of the following options:

    \itemize{
      \item{maxruns}{The number of Monte Carlo runs to be performed.}

      \item{report.interval}{The run interval at which to intermittently
	report, during the fitting process, the number of runs completed and
	the estimated time remaining.}
      }
    }

    \item{object}{A Gmp object}
    
  \item{\dots}{Other arguments to be passed on to the underlying function
    that performs the actual regression (see below).}
}

\details{

  Generalized Multilevel Permutation Models use regression to fit a
  model to the data under the original labeling of experimental
  conditions, and then compares the parameter estimates from this
  \dQuote{original fit} against null-hypothesis distributions obtained
  from fits to random relabelings of the data set.  These functions are
  especially useful for time-series data with difficult or unknown
  sources of dependency.  They can also be used for multivariate
  analyses.  Currently, p-values are provided for categorical predictors
  and interactions between categorical predictors and continuous
  covariates, but not for continuous covariates themselves.

  The \sQuote{formula} option to the model specifies the form of the
  model to be fit.  The left side of the formula identifies the
  dependent variable (DV), which is usually a vector of: real numbers
  (for continuous data), whole numbers (for counts), or 0s and 1s (for
  binomial data).  It can also be a factor with two-levels (binomial
  data) or more (multinomial data).  Lastly, one can also use the
  \code{cbind} syntax for binomial or multinomial DVs (see example
  below).

  The right hand side of the formula specifies the predictors in the
  regression model, using standard R format (see \alias{glm} for further
  information).  The predictors can include both design variables (IVs)
  and covariates.  Any IVs to be randomized must be identified in the
  \sQuote{ivars} option; all other variables are assumed to be
  covariates, and will not be randomized.  All IVs will be treated as
  categorical factors; continuous IVs are not allowed.  The IVs are
  internally coded using \sQuote{deviation} coding, such that the codes
  sum to zero.  The function \code{\link{getFactorCodes}} can be used to
  retrieve this internal coding.

  The \sQuote{family} option specifies the type of dependent variable.
  This option determines which of various regression subroutines are
  used to fit the model.  All of the options for generalized linear
  models (glms) are available (see \code{\link{family}} for details); if
  you choose one of these options then gmp will call the underlying
  regression function \code{\link{glm}} to perform the actual
  regression.  If the option \sQuote{multinomial} is chosen, then gmp
  will fit a multinomial model using \code{\link{multinom}} in package
  \pkg{nnet}.  In case tweaks to these underlying functions are
  necessary, it is possible to specify options to \code{gmp} or
  \code{gmpCreate} that will be passed along during the fitting process.
  For example, with multinomial models it is often necessary to increase
  the maximum number of iterations, which can be done by the argument
  \sQuote{maxit}.  User-defined models are not yet implemented, but left
  for future development.
}

\value{
  Each function returns a \alias{Gmp} object.
}

\references{ Barr, D. J. (in preparation).  Generalized Multilevel
  Permutation Models. }

\author{Dale J. Barr <dale.barr@ucr.edu>}

\seealso{\code{\link{kb07}}}

\examples{
# example using visual-world data with time as a continuous covariate
data(kb07)

# first let's fit the data using a binomial model
#
kb07.binom.gmp <- gmpCreate(cbind(T,NT) ~ t1*Speaker*Load | SubjID,
                       "binomial", kb07, ivars=c("Speaker","Load"))

# you will need to increase the number of runs (maxruns)
# to, say, 999 to get sensible results

kb07.binom.gmp <- gmpFit(kb07.binom.gmp, gmpCtrl(maxruns=19))

summary(kb07.binom.gmp)

# now let's do a more powerful multinomial analysis
# where we break out looks to blank regions and blinks
# into a separate "junk" category (X).  Region O
# ('other') will be our baseline category.  The main
# comparison, then is Target vs Other (TvO in the gmp output)

kb07.mnom.gmp <- gmpCreate(cbind(O,T,X) ~ t1*Speaker*Load | SubjID,
                     "multinomial", kb07, ivars=c("Speaker","Load"))

# again, you'll need to increase maxruns to get something sensible
kb07.mnom.gmp <- gmpFit(kb07.mnom.gmp, gmpCtrl(maxruns=19))

summary(kb07.mnom.gmp)
}
\keyword{htest}
\keyword{robust}
\keyword{design}
\keyword{ts}
\keyword{nonparametric}
\keyword{multivariate}
