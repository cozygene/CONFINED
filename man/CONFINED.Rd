% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CONFINED.R
\name{CONFINED}
\alias{CONFINED}
\title{Perform the CONFINED algorithm}
\usage{
CONFINED(X1, X2, t, k, thresh = 0.95, outfile = "", saveOP = TRUE)
}
\arguments{
\item{X1}{Input matrices. \eqn{m} by \eqn{n_1} matrix. Must have the same number of rows as \code{X2}}

\item{X2}{Input matrices. \eqn{m} by \eqn{n_2} matrix. Must have the same number of rows as \code{X1}}

\item{t}{Number of rows (e.g. methylation sites) to use}

\item{k}{Number of CONFINED components to produce (up to min\{\eqn{n_1}, \eqn{n_2 }\})}

\item{thresh}{Correlation threshold for selecting the number of canonical variables to use when generating the low-rank approximations in the feature selection step}

\item{outfile}{Prefix for saving the results}

\item{saveOP}{Boolean flag for saving the components/feature ranks to a txt file}
}
\value{
X1comps  -  \emph{CONFINED} components for \code{X1}. These capture shared variability between \code{X1} and \code{X2}

X2comps  -  \emph{CONFINED} components for \code{X2}. These capture shared variability between \code{X1} and \code{X2}
}
\description{
Generates components that capture sources of variability that are
shared between two datasets. See \code{https://github.com/cozygene/CONFINED} for more info.
}
