% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toy_model.R
\name{toy_model}
\alias{toy_model}
\title{Meta-genomic Toy Model}
\usage{
toy_model(
  n,
  cor,
  M,
  qdist,
  param = NULL,
  method = "pearson",
  force.positive = FALSE
)
}
\arguments{
\item{n}{number of observations.}

\item{cor}{correlation structure of final data, it can be a matrix or a number.
If it's a matrix it must be symmetric, with diagonal elements all equal to 1 and
with all values in range [-1,1]. The size of cor matrix implies the dimensions
of the resulting simulated data. Instead, if cor is a number, it
indicates the correlation between only two variables.}

\item{M}{magnification factor, real positive number that modify the
heterogeneity of the samples. In practice the first variable is multiplied by
the factor M.}

\item{qdist}{quantile function of the selected distribution, it must the vector
of probability as first argument called p. The other parameters names must match with
the ones present in param list.}

\item{param}{param named array or matrix/data.frame with parameters of the selected target
distribution indicated in qdist. If params is a matrix then it must contain
the parameters of the final distribution along the columns and rows express
each dimension equal to the dimension of the cor matrix, instead if it is a
vector the parameters are repeated for each dimension. The names of the columns,
or of the vector elements, must be equal to the function parameters of qdist.
If some parameters are missing then the default values of the qdist function
are chosen. Default value is NULL that indicates to use always the default
parameters of qdist.}

\item{method}{type of correlation used. Possible choices are "pearson",
"kendall", "spearman" (default pearson).}

\item{force.positive}{logical, indicates when to force all generated NorTA
data to positive numbers, adding the minimum to all others (this passage does
not affect correlations).}
}
\value{
returns an object of class "toy_model" containing:
\itemize{
\item cor input correlation matrix.
\item normal generated Gaussian data.
\item cor_normal correlation matrix of normal.
\item NorTA generated data of choosen distribution.
\item cor_NorTA correlation matrix of NorTA.
\item L1 relative abundances of NorTA.
\item cor_L1 correlation matrix of L1.
\item CLR matrix of abundances of clr transformed data from NorTA.
\item cor_CLR correlation matrix of CLR.
}
}
\description{
numerical simulations for meta-genomic data with user-chosen
correlation structure and sample heterogeneity.
}
