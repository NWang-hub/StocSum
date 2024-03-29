% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StocSum.R
\name{Cond.svt.pval}
\alias{Cond.svt.pval}
\title{Conditional analysis}
\usage{
Cond.svt.pval(
  meta.files.prefix,
  n.files = rep(1, length(meta.files.prefix)),
  tagChr,
  StartPos,
  EndPos,
  tagPos,
  MAF.range = c(1e-07, 0.5),
  miss.cutoff = 1,
  auto.flip = FALSE,
  nperbatch = 10000,
  tol = 1e-05
)
}
\arguments{
\item{n.files}{An integer vector with the same length as meta.files.prefix, indicating how many sets of intermediate files (.score.* and .var.*) are expected from each cohort, usually as the result of multi-threading in creating the intermediate files (default = rep(1, length(meta.files.prefix))).}

\item{tagChr}{The chromosome of associated variants condition on.}

\item{StartPos}{The start position of variants to be tested.}

\item{EndPos}{The end position of variants to be tested.}

\item{tagPos}{The position of associated variants condition on.}

\item{MAF.range}{A numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).}

\item{miss.cutoff}{The maximum missing rate allowed for a variant to be included (default = 1, including all variants).}

\item{auto.flip}{A logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.}

\item{nperbatch}{An integer for how many SNPs should be tested in a batch (default = 10000). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.}

\item{tol}{The threshold for the variance of a variant. If the variance of a variant less than the tolerance, its association test p-value will be NA (default = 1e-5).}

\item{meta.file.prefix}{A character vector for prefix of intermediate files (*.sample.* and *.resample.*).}
}
\value{
\code{Cond.svt.pval} returns a data frame with the following components:
\item{SNP}{SNP name.}
\item{chr}{chromosome name.}
\item{pos}{the genome location of SNP.}
\item{ref}{allele of reference.}
\item{alt}{alternative allele.}
\item{N}{total sample size.}
\item{missrate}{missing rate of variants.}
\item{altfreq}{alternative allele frequency.}
\item{SCORE}{the summary score of the alternaive allele.}
\item{VAL}{the variance of the summary score.}
\item{PVAL}{the conditional p-value of the suammry score.}
}
\description{
Conditional analysis.
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
out.prefix <- "test"
gdsfile <- system.file("extdata", "geno.gds", package = "StocSum")
G.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix,MAF.range=c(0,0.5), miss.cutoff = 1)
out <- Cond.svt.pval(out.prefix, n.files = 1, tagChr = 1, StartPos = 1, EndPos = 100, tagPos = 82, MAF.range=c(0,0.5), miss.cutoff = 1)
print(out)
}
}
\seealso{
\code{glmmkin2randomvec}, \code{G.stat}
}
\author{
Han Chen, Nannan Wang
}
\keyword{analysis}
\keyword{conditional}
