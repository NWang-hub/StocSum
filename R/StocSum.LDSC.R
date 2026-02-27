

#' Calculate LDscore
#' @description Use summary statistics and stochastic statistics from G.stat to estimate LDscore.
#' @param meta.file.prefix A character for prefix of intermediate files (*.sample.* and *.resample.*). 
#' @param n.files An integer indicating how many sets of intermediate files (.score.* and .var.*), usually as the result of multi-threading in creating the intermediate files (default = 1)).
#' @param saveGDS The format of intermediate files (*.sample.1 and *.resample.1). If "TRUE", then the intermiate files will be saved as GDS format to save storage, otherwise the files will be saved as binary format.
#' @param N.randomvec The number of random vectors to generate (default = 1000). Same as the values in last step of "glmmkin2randomvec"
#' @param N.randomvec The number of replicates to simulate the random vectors in \code{StocSum.LDSC.glmmkin2randomvec} (default = 1000).
#' @param MAF.range A numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
#' @param miss.cutoff The maximum missing rate allowed for a variant to be included (default = 1, including all variants).
#' @param wind.b A positive integer define the window size in bases (default = 1000000).
#' @param use.minor.allele A logical switch for whether to use the minor allele (instead of the alt allele) as the coding allele (default = FALSE). It does not change SKAT results, but Burden (as well as SKAT-O and hybrid test to combine the burden test and SKAT) will be affected. Along with the MAF filter, this option is useful for combining rare mutations, assuming rare allele effects are in the same direction. Use with caution, as major/minor alleles may flip in different cohorts. In that case, minor allele will be determined based on the allele frequency in the combined samples.
#' @param auto.flip A logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
#' @param nperbatch An integer for how many SNPs to be included in a batch (default = 10000). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.
#' @param ncores A positive integer indicating the number of cores to be used in parallel computing (default = 1).
#' @return \code{LDSC.win} returns a data frame with the following components:
#' \item{chr}{chromosome name.}
#' \item{pos}{the genome location of SNP.}
#' \item{N}{total sample size.}
#' \item{altfreq}{alternative allele frequency.}
#' \item{LDscore}{LD score for each SNP.}
#' @reference 
#' @author Han Chen, Nannan Wang
#' @seealso \code{LDSC.glmmkin2randomvec}, \code{G.stat}
#' @examples
#' \donttest{
#' library(StocSum)
#' library(GMMAT)
#' library(data.table)
#' data(example)
#' attach(example)
#' seed <- 12345
#' set.seed(seed)
#' GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "StocSum")
#' GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
#' nullmod <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id", family = binomial(link = "logit"))
#' obj <- LDSC.glmmkin2randomvec(nullmod)
#' out.prefix <- "test"
#' gdsfile <- system.file("extdata", "geno.gds", package = "StocSum")
#' G.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix, MAF.range=c(0,0.5), miss.cutoff = 1)
#' out<-LDSC.win(out.prefix, use.minor.allele = FALSE, auto.flip = FALSE, wind.b = 1000000, nperbatch = 10000)
#' }
#' @keywords LD Score
#' @export

MAF.weights.beta.fun <- function(freq, beta1, beta2) {
  freq[freq > 0.5] <- 1 - freq[freq > 0.5]
  ifelse(freq <= 0, 0, dbeta(freq, beta1, beta2))
}

LDSC.win <- function(meta.files.prefix, n.files = 1, saveGDS = TRUE, MAF.range = c(1e-2, 0.5), miss.cutoff = 1, wind.b=1000000, extract=NULL, exclude=NULL, use.minor.allele = FALSE, auto.flip = FALSE, nperbatch = 10000, ncores = 1)
{
  MAF.weights.beta <- c(0.5, 0.5)
  if (!is.null(extract)) {
    extract<-read.table(extract,header=F,sep="\t")
    extract<-extract$V1
  }
  if (!is.null(exclude)) {
    exclude<-read.table(exclude,header=F,sep="\t")
    exclude<-exclude$V1
  }
  if(.Platform$endian!="little") stop("Error: platform must be little endian.")
  group.info <- NULL
  for(j in 1:n.files) {
    tmp <- fread(paste0(meta.files.prefix, "_part", j,".freq"), header=TRUE, data.table = FALSE)
    if (class(tmp) == "try-error") {
      stop(paste0("Error: cannot read ", meta.files.prefix, "_part", j, ".freq!"))
    }
    tmp$file <- j
    tmp$variant.idx <- 1:nrow(tmp)
    if(nrow(tmp) > 0) group.info <- rbind(group.info, tmp)
    rm(tmp)
  }
  group.info <- group.info[order(group.info$chr, group.info$pos), ]
  if (!("missrate" %in% colnames(group.info))) group.info$missrate<-0.01
  n.groups.all<-nrow(group.info)
  group.info$idx <- 1:nrow(group.info)
  group.info$pos.end<-group.info$pos+wind.b
  group.info$pos.end[group.info$pos.end>group.info$pos[n.groups.all]]<-group.info$pos[n.groups.all]
  group.info$pos.start<-group.info$pos-wind.b
  group.info$pos.start[group.info$pos.start<group.info$pos[1]]<-group.info$pos[1]
  snp_group1<-paste(group.info$chr,group.info$pos,group.info$ref,group.info$alt,sep=":")
  snp_group2<-paste(group.info$chr,group.info$pos,group.info$alt,group.info$ref,sep=":")
  snp2<-c()
  if (!is.null(extract)) {
    snp2_1<-intersect(extract,snp_group1)
    snp2_2<-intersect(extract,snp_group2)
    if(length(snp2_1)>=length(snp2_2)) {
      snp_group<-snp_group1
      snp2<-snp2_1
    } else {
      snp_group<-snp_group2
      snp2<-snp2_2
    } 
    rm(snp2_1, snp2_2)
  } else if (!is.null(exclude)) {
    snp2_1<-intersect(exclude,snp_group1)
    snp2_2<-intersect(exclude,snp_group2)
    if(length(snp2_1)>=length(snp2_2)) {
      snp_group<-snp_group1
      snp2<-snp2_1 
    } else {
      snp_group<-snp_group2 
      snp2<-snp2_2
    }
    rm(snp2_1, snp2_2)
  }
  if (!is.null(extract) & length(snp2)==0) stop(paste0("Error: no overlap variants between extract and genotype"))
  if (!is.null(exclude) & length(snp2)==0) stop(paste0("Error: no overlap variants between exclude and genotype"))
  if(!is.null(extract)) extract.idx<-group.info$idx[match(extract,snp_group)]
  if(!is.null(exclude)) exclude.idx<-group.info$idx[match(exclude,snp_group)]
  rm(snp_group, snp_group1, snp_group2)
  scores <- cons <- vector("list", 1)
  N.resampling <- rep(0, 1)
  if (!saveGDS){
    cons <- file(paste0(meta.files.prefix, "_part1.stoc"), "rb")
    N.resampling <- readBin(cons, what = "integer", n = 1, size = 4)
  } else {
    cons <- paste0(meta.files.prefix, "_part1")
    fnamej<-paste0(cons,".stoc")             
    N.resampling <-readgds(fnamej,dtype="Nrandomvec")
  }
  N.randomvec<-N.resampling
  current.lines <- current.cons <- rep(1, 1)
  all.out <- NULL
  ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
  if(ncores > 1) {
    doMC::registerDoMC(cores = ncores)
    n.groups.percore <- (n.groups.all-1) %/% ncores + 1
    n.groups.percore_1 <- n.groups.percore * ncores - n.groups.all
    b <- NULL
    all.out <- foreach(b=1:ncores, .combine=rbind, .multicombine = TRUE, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
      idx0 <- if(b <= n.groups.percore_1) group.info$idx[((b-1)*(n.groups.percore-1)+1):(b*(n.groups.percore-1))] else group.info$idx[(n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1-1)*n.groups.percore+1):(n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1)*n.groups.percore)]
      n.groups <- length(idx0)
      nbatch.flush <- (n.groups-1) %/% nperbatch + 1
      for(i in 1:nbatch.flush) {
        ii<-idx0[((i-1)*nperbatch+1):min((i*nperbatch),n.groups)]
        itmp.idx<-group.info$idx[ii]
        dis<-group.info$pos-group.info[group.info$idx==itmp.idx[1],]$pos.start
        tmp.idx1<-group.info$idx[dis==min(dis[dis>=0])]
        dis<-group.info$pos-group.info[group.info$idx==itmp.idx[length(itmp.idx)],]$pos.end
        tmp.idx2<-group.info$idx[dis==min(dis[dis>=0])]
        tmp.idx<-group.info$idx[tmp.idx1:tmp.idx2]
        rm(dis)
        U.list <- V.list <- vector("list", 1)
        variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
        tmp.scores <- group.info[group.info$idx %in% tmp.idx, , drop = FALSE]
        if(any(tmp.include <- !is.na(tmp.scores$altfreq))) {
          U.list <- tmp.scores[tmp.include, , drop = FALSE]
          tmp.V <- matrix(NA, sum(tmp.include), N.resampling)
          if (!saveGDS){
            for(ij in 1:sum(tmp.include)){
              if(U.list$file[ij]!=current.cons) {
                close(cons)
                current.cons <- U.list$file[ij]
                cons <- file(paste0(meta.files.prefix, "_part", current.cons,".stoc"), "rb")
                tmp.N.resampling <- readBin(cons, what = "integer", n = 1, size = 4)
                if(tmp.N.resampling != N.resampling) stop(paste0("Error: N.resampling in ", meta.files.prefix, "_stoc_part", current.cons, " does not match that in ",meta.files.prefix, "_stoc_part.1"))
              }
              current.cons <- U.list$file[ij]
              cons <- file(paste0(meta.files.prefix, "_part", current.cons,".stoc"), "rb")
              ind <- 4+4*N.resampling*(U.list$variant.idx[ij]-1)
              seek(cons, where = ind, origin = "start", rw = "read")
              tmp.V[ij,] <- readBin(cons, what = "numeric", n = N.resampling, size = 4)
            }    
          } else if (saveGDS) {
            jfile_num<-unique(U.list$file)
            for (kk in jfile_num){
              kk_Uloc<-which(U.list$file==kk)
              current.cons <- kk
              cons <- paste0(meta.files.prefix, "_part", current.cons)
              fnamej<-paste0(cons,".stoc")
              tmp.N.resampling <- readgds(fnamej,dtype="Nrandomvec")
              if(tmp.N.resampling != N.resampling) stop(paste0("Error: N.resampling in ", meta.files.prefix, "_stoc_part", current.cons, " does not match that in ",meta.files.prefix, "_stoc_part1"))
              loc_Ugds<-U.list$variant.idx[kk_Uloc]
              tmp.V_kk <- readgds(fnamej,dtype="U",start=c(1,min(loc_Ugds)),count=c(-1,max(loc_Ugds)-min(loc_Ugds)+1))
              tmp.V[kk_Uloc,]<-t(tmp.V_kk[,loc_Ugds-min(loc_Ugds)+1])
            }
          }
          V.list <- tmp.V / sqrt(N.resampling)
          rm(tmp.V)
          variant.indices <- c(variant.indices, U.list$idx)
          tmp.N <- c(tmp.N, U.list$N)
          tmp.Nmiss <- c(tmp.Nmiss, U.list$N * U.list$missrate/(1-U.list$missrate))
          tmp.AC <- c(tmp.AC, 2*U.list$N*U.list$altfreq)
        }
        tmp.variant.indices <- variant.indices
        if (!is.null(extract)) variant.indices<-intersect(variant.indices,extract.idx)
        if (!is.null(exclude)) variant.indices<-setdiff(variant.indices,exclude.idx)
        if(length(variant.indices) == 0) next
        variant.indices <- sort(unique(variant.indices))
        N <- sapply(variant.indices, function(x) sum(tmp.N[tmp.variant.indices==x]))
        Nmiss <- sapply(variant.indices, function(x) sum(tmp.Nmiss[tmp.variant.indices==x]))
        AF <- sapply(variant.indices, function(x) sum(tmp.AC[tmp.variant.indices==x]))/2/N
        include <- (Nmiss/(N+Nmiss) <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
        rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices)
        if(sum(include) == 0) next
        variant.indices <- variant.indices[include]
        N<-N[include]
        AF<-AF[include]
        itmp.idx<-intersect(itmp.idx,variant.indices)
        out <- group.info[match(itmp.idx,group.info$idx),]
        n.p <- length(variant.indices)
        V <- matrix(0, n.p, sum(N.resampling))
        if(!is.null(U.list) & !is.null(V.list)) {
          IDX <- match(U.list$idx,variant.indices)
          if(sum(!is.na(IDX)) == 0) next
          IDX2 <- which(!is.na(IDX))
          IDX <- IDX[IDX2]
          V[IDX, 1:N.resampling] <- V.list[IDX2,]
        }
        n.batchi<-length(itmp.idx)
        tmp.weight <- MAF.weights.beta.fun(AF, MAF.weights.beta[1], MAF.weights.beta[2])
        if(use.minor.allele) tmp.weight[AF > 0.5] <- -tmp.weight[AF > 0.5]
        # V <- V*tmp.weight*(gamma(MAF.weights.beta[1])*gamma(MAF.weights.beta[2])/gamma(MAF.weights.beta[1]+MAF.weights.beta[2])/sqrt(2))
         rowvar_V <- apply(V, 1, var)
         rowmean_V <- apply(V, 1, mean)
         V<-(V-rowmean_V)/sqrt(rowvar_V)
        
        jloc <- match(itmp.idx,variant.indices)
        V_self <- V[jloc,]
        jN <-N[jloc]
        pos_self <- group.info$pos[itmp.idx]
        pos_full <- group.info$pos[variant.indices]
        mask <- outer(pos_self,pos_full,function(a, b) abs(a - b))
        mask <- ifelse(mask < wind.b, 1, 0)
        mat <- tcrossprod(V_self,V) 
        mat <- mat * mask/(N.randomvec-1)
        nvariants <- rowSums(mask)
        jLDscore <- (1+1/(N.randomvec-2)+1/(jN-2))*rowSums(mat^2)-nvariants/(N.randomvec-2)-nvariants/(jN-2)
        out$LDscore<-jLDscore
        all.out <- rbind(all.out, out)
      }
      all.out<-all.out[,c("chr","pos","ref","alt","N","missrate","altfreq","LDscore")]
      return(all.out)
    }
  } else { # use a single core
    n.groups <- n.groups.all  #length(variant.id1)
    nbatch.flush <- (n.groups-1) %/% nperbatch + 1
    for(i in 1:nbatch.flush) {
      itmp.idx<-group.info$idx[((i-1)*nperbatch+1):min((i*nperbatch),n.groups)]
      dis<-group.info$pos-group.info[group.info$idx==itmp.idx[1],]$pos.start
      tmp.idx1<-group.info$idx[dis==min(dis[dis>=0])]
      dis<-group.info$pos-group.info[group.info$idx==itmp.idx[length(itmp.idx)],]$pos.end
      tmp.idx2<-group.info$idx[dis==min(dis[dis>=0])]
      tmp.idx<-group.info$idx[tmp.idx1:tmp.idx2]
      rm(dis)
      U.list <- V.list <- vector("list", 1)
      variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
      tmp.scores <- group.info[group.info$idx %in% tmp.idx, , drop = FALSE]
      if(any(tmp.include <- !is.na(tmp.scores$altfreq))) {
        U.list <- tmp.scores[tmp.include, , drop = FALSE]
        tmp.V <- matrix(NA, sum(tmp.include), N.resampling)
        if (!saveGDS){
          for(ij in 1:sum(tmp.include)){
            if(U.list$file[ij]!=current.cons) {
              close(cons)
              current.cons <- U.list$file[ij]
              cons <- file(paste0(meta.files.prefix, "_part", current.cons,".stoc"), "rb")
              tmp.N.resampling <- readBin(cons, what = "integer", n = 1, size = 4)
              if(tmp.N.resampling != N.resampling) stop(paste0("Error: N.resampling in ", meta.files.prefix, "_part", current.cons, " does not match that in ",meta.files.prefix, "_stoc_part1"))
              current.lines <- 1
            }
            if(U.list$variant.idx[ij]!=current.lines) seek(cons, where = 4*N.resampling*(U.list$variant.idx[ij]-current.lines), origin = "current", rw = "read")
            tmp.V[ij,] <- readBin(cons, what = "numeric", n = N.resampling, size = 4)
            current.lines <- U.list$variant.idx[ij]+1
          }    
        } else if (saveGDS) {
          jfile_num<-unique(U.list$file)
          for (kk in jfile_num){
            kk_Uloc<-which(U.list$file==kk)
            current.cons <- kk
            cons <- paste0(meta.files.prefix, "_part", current.cons)
            fnamej<-paste0(cons,".stoc")
            tmp.N.resampling <- readgds(fnamej,dtype="Nrandomvec")
            if(tmp.N.resampling != N.resampling) stop(paste0("Error: N.resampling in ", meta.files.prefix, "_stoc_part", current.cons, " does not match that in ",meta.files.prefix, "_stoc_part1"))
            loc_Ugds<-U.list$variant.idx[kk_Uloc]
            tmp.V_kk <- readgds(fnamej,dtype="U",start=c(1,min(loc_Ugds)),count=c(-1,max(loc_Ugds)-min(loc_Ugds)+1))
            tmp.V[kk_Uloc,]<-t(tmp.V_kk[,loc_Ugds-min(loc_Ugds)+1])
          }
        }
        V.list <- tmp.V / sqrt(N.resampling)
        rm(tmp.V)
        variant.indices <- c(variant.indices, U.list$idx)
        tmp.N <- c(tmp.N, U.list$N)
        tmp.Nmiss <- c(tmp.Nmiss, U.list$N * U.list$missrate/(1-U.list$missrate))
        tmp.AC <- c(tmp.AC, 2*U.list$N*U.list$altfreq)
      }
      tmp.variant.indices <- variant.indices
      if (!is.null(extract)) variant.indices<-intersect(variant.indices,extract.idx)
      if (!is.null(exclude)) variant.indices<-setdiff(variant.indices,exclude.idx)
      if(length(variant.indices) == 0) next
      variant.indices <- sort(unique(variant.indices))
      N <- sapply(variant.indices, function(x) sum(tmp.N[tmp.variant.indices==x]))
      Nmiss <- sapply(variant.indices, function(x) sum(tmp.Nmiss[tmp.variant.indices==x]))
      AF <- sapply(variant.indices, function(x) sum(tmp.AC[tmp.variant.indices==x]))/2/N
      include <- (Nmiss/(N+Nmiss) <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
      rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices)
      if(sum(include) == 0) next
      variant.indices <- variant.indices[include]
      N<-N[include]
      AF<-AF[include]
      itmp.idx<-intersect(itmp.idx,variant.indices)
      out <- group.info[match(itmp.idx,group.info$idx),]
      n.p <- length(variant.indices)
      V <- matrix(0, n.p, sum(N.resampling))
      if(!is.null(U.list) & !is.null(V.list)) {
        IDX <- match(U.list$idx,variant.indices)
        if(sum(!is.na(IDX)) == 0) next
        IDX2 <- which(!is.na(IDX))
        IDX <- IDX[IDX2]
        V[IDX, 1:N.resampling] <- V.list[IDX2,]
      }
      n.batchi<-length(itmp.idx)
      tmp.weight <- MAF.weights.beta.fun(AF, MAF.weights.beta[1], MAF.weights.beta[2])
      if(use.minor.allele) tmp.weight[AF > 0.5] <- -tmp.weight[AF > 0.5]
      # V <- V*tmp.weight*(gamma(MAF.weights.beta[1])*gamma(MAF.weights.beta[2])/gamma(MAF.weights.beta[1]+MAF.weights.beta[2])/sqrt(2))
      rowvar_V <- apply(V, 1, var)
      rowmean_V <- apply(V, 1, mean)
      V<-2*(V-rowmean_V)/sqrt(rowvar_V)
      
      jLDscore<-rep(NA,n.batchi)
      jloc <- match(itmp.idx,variant.indices)
      V_self <- V[jloc,]
      jN <-N[jloc]
      pos_self <- group.info$pos[itmp.idx]
      pos_full <- group.info$pos[variant.indices]
      mask <- outer(pos_self,pos_full,function(a, b) abs(a - b))
      mask <- ifelse(mask < wind.b, 1, 0)
      mat <- tcrossprod(V_self,V) 
      mat <- mat * mask/(N.randomvec-1)
      nvariants <- rowSums(mask)
      jLDscore <- (1+1/(N.randomvec-2)+1/(jN-2))*rowSums(mat^2)-nvariants/(N.randomvec-2)-nvariants/(jN-2)
      out$LDscore<-jLDscore
      all.out <- rbind(all.out, out)
    }
  }
  if (!saveGDS){
    close(cons)
  }
  all.out<-all.out[,c("chr","pos","ref","alt","N","missrate","altfreq","LDscore")]
  return(all.out)
}

