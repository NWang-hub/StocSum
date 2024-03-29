% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StocSum.GE.R
\name{GE.svt.pval}
\alias{GE.svt.pval}
\title{GLMM based single variant tests for gene-environment interactions}
\usage{
GE.svt.pval(
  meta.files.prefix,
  out.file,
  n.files = rep(1, length(meta.files.prefix)),
  n.pheno = 1,
  interaction,
  interaction.covariates = NULL,
  MAF.range = c(1e-07, 0.5),
  miss.cutoff = 1,
  auto.flip = FALSE,
  nperbatch = 100
)
}
\arguments{
\item{out.file}{The output file name.}

\item{n.files}{An integer vector with the same length as meta.files.prefix, indicating how many sets of intermediate files (.score.* and .var.*) are expected from each cohort, usually as the result of multi-threading in creating the intermediate files (default = rep(1, length(meta.files.prefix))).}

\item{n.pheno}{An integer indicating the number of phenotypes in multiple phenotype analysis (for single phenotype analysis, \code{n.pheno = 1}) (default = 1).}

\item{interaction}{A numeric or a character vector indicating the environmental factors. If a numberic vector, it represents which indices in the order of covariates are the environmental factors; if a character vector, it represents the variable names of the environmental factors.}

\item{interaction.covariates}{A numeric or a character vector indicating the interaction covariates. If a numeric vector, it represents which indices in the order of covariates are the interaction covariates; if a character vector, it represents the variable names of the interaction covariates.}

\item{MAF.range}{A numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).}

\item{miss.cutoff}{The maximum missing rate allowed for a variant to be included (default = 1, including all variants).}

\item{auto.flip}{A logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.}

\item{nperbatch}{An integer for how many SNPs should be tested in a batch (default = 100). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user's system.}

\item{meta.file.prefix}{A character vector for prefix of intermediate files (*.sample.* and *.resample.*).}
}
\description{
Use the summary statistic and stochastic statistics from \code{GE.stat} to perform single variant main effect tests, gene-environment interaction tests, and joint tests for association with genotypes in a GDS file (.gds). 7 user-defined tests are included: Main effect variance component test (MV), Main effect hybrid test of burden and variance component test using Fisher's method (MF), Interaction variance component test (IV), Interaction hybrid test of burden and variance component test using Fisher's method (IF), Joint variance component test (JV), Joint hybrid test of burden and variance component test using Fisher's method (JF), and Joint hybrid test of burden and variance component test using double Fisher's procedures (JD).
}
\examples{
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
if(!is.null(nullmod$P)){
  obj <- glmmkin2randomvec(nullmod)
}else{
  kinship.chol <- chol(GRM)
  obj<-glmmkin2randomvec(nullmod, Z = list(t(kinship.chol)))
}
out.prefix <- "test.GE.svt"
gdsfile <- system.file("extdata", "geno.gds", package = "StocSum")
interaction <- c("sex")
GE.stat(obj, interaction = interaction, geno.file = gdsfile, meta.file.prefix = out.prefix,MAF.range=c(0,0.5), miss.cutoff = 1)
out.file<-paste0(out.prefix,".out")
GE.svt.pval(meta.files.prefix = out.prefix, out.file, n.files = 1, n.pheno = 1, interaction=interaction, MAF.range=c(1e-7,0.5), miss.cutoff = 1, auto.flip=FALSE)
}
}
\seealso{
\code{glmmkin2randomvec}, \code{GE.stat}
}
\author{
Han Chen, Nannan Wang
}
\keyword{gene-environment}
\keyword{generalized}
\keyword{interaction}
\keyword{linear}
\keyword{mixed}
\keyword{model,}
