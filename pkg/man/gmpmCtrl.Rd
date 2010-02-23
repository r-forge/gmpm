\name{gmpmCtrl}
\alias{gmpmCtrl}
\title{Return a gmpmCtrl Object with Default Settings}
\description{
This is a convenience function that returns a gmpmControl object for
controlling Monte Carlo runs, with default settings.  (Will probably be
deprecated in future versions of package \pkg{gmpm}.)
}
\usage{
gmpmCtrl(maxruns = 999, report.interval = 10, outfile = NULL)
}
\arguments{
  \item{maxruns}{total number of Monte Carlo runs}
  \item{report.interval}{interval at which to report back to user}
  \item{outfile}{(for future development)}
}
\value{
  A list, with elements as indicated above.
}
