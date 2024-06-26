% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre_treat.R
\name{remove_outliers}
\alias{remove_outliers}
\title{remove_outliers}
\usage{
remove_outliers(x, n, na.rm = TRUE, ...)
}
\arguments{
\item{x}{numeric vector}

\item{n}{the coefficient of IQR}

\item{na.rm}{logical; if true, any NA and NaN's are removed from x before the quantiles are computed}

\item{...}{further arguments passed to or from other methods}
}
\description{
remove missing values based on IQR principles
}
