% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{pairwise.sum.decreasing}
\alias{pairwise.sum.decreasing}
\title{Return the component-wise sum of decreasing vectors in which NA is taken to mean that the vector has stopped decreasing
and we can use the last non-na element. Where both are NA, leave as NA.}
\usage{
pairwise.sum.decreasing(x, y)
}
\description{
Return the component-wise sum of decreasing vectors in which NA is taken to mean that the vector has stopped decreasing
and we can use the last non-na element. Where both are NA, leave as NA.
}
