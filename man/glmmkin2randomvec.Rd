% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StocSum.GE.R, R/StocSum.R
\name{glmmkin2randomvec}
\alias{glmmkin2randomvec}
\title{Create random vectors for a glmmkin object}
\usage{
glmmkin2randomvec(
  obj,
  Z = NULL,
  N.randomvec = 1000,
  group.idx = NULL,
  cluster.idx = NULL,
  robust = FALSE
)

glmmkin2randomvec(
  obj,
  Z = NULL,
  N.randomvec = 1000,
  group.idx = NULL,
  cluster.idx = NULL,
  robust = FALSE
)
}
\arguments{
\item{obj}{A glmmkin object.}

\item{Z}{A list of design matrices for the random effects. The length must match the number of variance components.}

\item{N.randomvec}{The number of random vectors to generate (default = 1000).}

\item{group.idx}{A length N index vector showing which observation belongs to which variance group, for heteroscedastic linear mixed models (default = NULL for homoscedastic linear mixed models).}

\item{cluster.idx}{A length N index vector showing which observation belongs to which cluster (default = NULL for no clusters).}

\item{robust}{A logical switch: whether robust variance should be used (default = FALSE).}
}
\value{
A list of class glmmkin.randomvec
\item{theta}{inherited from the glmmkin object. A vector or a list of variance component parameter estimates.}
\item{scaled.residuals}{inherited from the glmmkin object. A vector or a matrix for the scaled residuals, calculated as the original residuals divided by the dispersion parameter (in heteroscedastic linear mixed models, corresponding residual variance estimates by each group).}
{random.vectors}{A random matrix with dimensions equal to the sample size multiplied by \code{N.randomvec}.}
\item{X}{inherited from the glmmkin object. Model matrix for the fixed effects.}
\item{id_include}{inherited from the glmmkin object. A vector indivating the samples included in model fit.}

A list of class glmmkin.randomvec
\item{theta}{inherited from the glmmkin object. A vector or a list of variance component parameter estimates.}
\item{scaled.residuals}{inherited from the glmmkin object. A vector or a matrix for the scaled residuals, calculated as the original residuals divided by the dispersion parameter (in heteroscedastic linear mixed models, corresponding residual variance estimates by each group).}
\item{random.vectors}{A random matrix with dimensions equal to the sample size multiplied by \code{N.randomvec}.}
\item{X}{inherited from the glmmkin object. Model matrix for the fixed effects.}
\item{id_include}{inherited from the glmmkin object. A vector indivating the samples included in model fit.}
}
\description{
Generate random vectors from multivariate normal distribution with mean 0 and covariance matrix P.

Generate random vectors from multivariate normal distribution with mean 0 and covariance matrix P.
}
\examples{
\donttest{
library(StocSum)
library(GMMAT)
data(example)
attach(example)
seed <- 12345
set.seed(seed)
GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "StocSum")
GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
nullmod <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id", family = binomial(link = "logit"))
if(!is.null(nullmod$P)){
  obj <- glmmkin2randomvec(nullmod)
}else{
  kinship.chol <- chol(GRM)
  obj<-glmmkin2randomvec(nullmod, Z = list(t(kinship.chol)))
}
}
\donttest{
library(StocSum)
library(GMMAT)
library(data.table)
data(example)
attach(example)
seed <- 12345
set.seed(seed)
GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "StocSum")
GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
nullmod <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id", family = binomial(link = "logit"))
}
}
\author{
Han Chen, Nannan Wang
}
\keyword{random}
\keyword{vector}
