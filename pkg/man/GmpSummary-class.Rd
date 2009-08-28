\name{GmpSummary-class}
\docType{class}
\alias{GmpSummary-class}
\alias{initialize,GmpSummary-method}
\alias{show,GmpSummary-method}
\title{Class GmpSummary}
\description{Summary information for a GMPM.}

\section{Objects from the Class}{

  Objects can be created by a call to \code{\link{summary}}.

}

\section{Slots}{
  \describe{

    \item{\code{gmpInfo}:}{Object of class \code{"list"}; basic
      information about the model.}
    
    \item{\code{gmpMainSum}:}{Object of class \code{"list"}; the main
      summary of p-values, including multiparameter tests where
      appropriate.}
    
    \item{\code{gmpRegSum}:}{Object of class \code{"list"}; optional
    summary of the regression, not included by default (see
    \sQuote{showReg} below).}
    
    \item{\code{showReg}:}{Object of class \code{"logical"} whether or
      not the regression summary is to be calculated and displayed.}
    
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "GmpSummary")}: ... }
    \item{show}{\code{signature(object = "GmpSummary")}: ... }
	 }
}
\seealso{
  \code{\link{gmp}}, \code{\link{Gmp-class}}
}
\keyword{classes}
