% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{get_glm_coefs}
\alias{get_glm_coefs}
\title{get glm estimates}
\usage{
get_glm_coefs(data)
}
\arguments{
\item{data}{list with AE count (A), confounders (X), vaccine indication (V), number of each categories}
}
\value{
list of glm estimated beta, alpha coefficients with their p-values and adsusted q-values
}
\description{
get glm estimates
}
