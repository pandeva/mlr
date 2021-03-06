% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createDummyFeatures.R
\name{createDummyFeatures}
\alias{createDummyFeatures}
\title{Generate dummy variables for factor features.}
\usage{
createDummyFeatures(obj, target = character(0L), method = "1-of-n",
  cols = NULL)
}
\arguments{
\item{obj}{[\code{data.frame} | \code{\link{Task}}]\cr
Input data.}

\item{target}{[\code{character(1)} | \code{character(2)} | \code{character(n.classes)}]\cr
Name(s) of the target variable(s).
Only used when \code{obj} is a data.frame, otherwise ignored.
If survival analysis is applicable, these are the names of the survival time and event columns,
so it has length 2.
For multilabel classification these are the names of logical columns that indicate whether
a class label is present and the number of target variables corresponds to the number of
classes.}

\item{method}{[\code{character(1)}]\cr
Available are:\cr
\dQuote{1-of-n}: For n factor levels there will be n dummy variables.\cr
\dQuote{reference}: There will be n-1 dummy variables leaving out the first factor level of each variable.\cr}

\item{cols}{[\code{character}]\cr
Columns to create dummy features for. Default is to use all columns.}
}
\value{
[\code{data.frame} | \code{\link{Task}}]. Same type as \code{obj}.
}
\description{
Replace all factor features with their dummy variables. Internally \code{\link{model.matrix}} is used.
Non factor features will be left untouched and passed to the result.
}
\seealso{
Other eda_and_preprocess: \code{\link{capLargeValues}},
  \code{\link{dropFeatures}},
  \code{\link{mergeSmallFactorLevels}},
  \code{\link{normalizeFeatures}},
  \code{\link{removeConstantFeatures}},
  \code{\link{summarizeColumns}}
}
