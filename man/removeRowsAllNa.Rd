% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre_treat.R
\name{removeRowsAllNa}
\alias{removeRowsAllNa}
\title{removeRowsAllNa}
\usage{
removeRowsAllNa(object)
}
\arguments{
\item{object}{a data.frame or matrix}
}
\description{
Delete rows that are all NA
}
\examples{
\dontrun{
 x <- matrix(c(rep(NA,3),rep(1,6)),nrow = 3)
 x <- removeRowsAllNa(x)
}
}
