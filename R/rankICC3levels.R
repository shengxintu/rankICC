#' Rank-based ICC with three hierarchies
#'
#' \code{rankICC3levels} computes the rank-based intraclass correlation coefficient (ICC) with three hierarchies. Starting from the innermost level, the three levels are named level 1, level 2, and level 3. The rank-based ICC at level 2 evaluates the correlation between a random pair from the same level-2 unit. The rank-based ICC at level 3 evaluates the correlation between a random pair from the same level-3 unit but different level-2 units.
#'
#' @param x a numeric vector of observations.
#' @param level2 a index vector of the level-2 units.
#' @param level3 a index vector of the level-3 units.
#' @param weights a character string indicating which weighting method is used. Or an optional vector of user-defined weights to be used. Should be one of the strings \code{"level1"}, \code{"level2"}, \code{"level3"}, or a numeric vector. See Details.
#' @param conf.int numeric specifying confidence interval coverage.
#' @param fisher logical indicating whether to apply fisher transformation to compute confidence intervals.
#' @details \code{"level1"} assigns equal weights to level-1 units; \eqn{p_{ijk}=1/(\sum_{i=1}^n\sum_{j=1}^{n_i}m_{ij})}, where \eqn{n} is the total number of level-3 units, \eqn{n_i} is the number of level-2 units in the \eqn{i}th level-3 unit, and \eqn{m_{ij}} is the number of level-1 units in the \eqn{j}th level-2 unit and the \eqn{i}th level-3 unit. \code{"level2"} assigns equal weights to level-2 units; \eqn{p_{ijk}=1/(m_{ij}\sum_{i=1}^n n_{i})}. \code{"level3"} assigns equal weights to level-3 units; \eqn{p_{ijk}=1/(nn_{i}m_{ij})}.
#' @return a matrix with two rows. The first row is for rank-based ICC at level 2 and the second row is for rank-based ICC at level 3. Each row with following components.
#' \tabular{ll}{
#'   \code{rankICC} \tab the rank-based ICC. \cr
#'   \tab \cr
#'   \code{SE} \tab the standard error. \cr
#'   \tab \cr
#'   \code{Lower, Upper} \tab the lower and upper bound of the confidence interval.\cr
#' }
#' @examples
#' \dontrun{
#' rankICC3levels(x, level2, level3, weights = "level2")
#' }
#' @export

rankICC3levels <- function(x, level2, level3, weights = c("level1", "level2", "level3"),
                        conf.int = 0.95, fisher = FALSE){
    level3 <- factor(level3, levels = unique(level3))
    level2 <- as.character(level2)
    idx <- order(level3, level2)
    x <- x[idx]
    level3 <- level3[idx]
    level2 <- level2[idx]
    level32 <- paste(level3, level2, sep="-")
    level32 <- factor(level32, levels = unique(level32))
    #number of obs in each level3
    ki <- tabulate(level3)
    #number of obs in each level2
    kij <- unlist(tapply(level2, level3, table))
    #number of level2 in each level3
    mi <- unlist(lapply(tapply(level2, level3, unique), length))
    #number of obs
    N <- length(level3)
    if(is.numeric(weights)){
      if(length(weights) != length(x)) stop("user-defined weights do not have the same length as the observations")
      else{
          pijk <- weights[idx]
          pij <- rep(tapply(pijk, level32, sum), kij)
          pi <- rep(tapply(pijk, level3, sum), ki)
      }
    }
    else if(weights[1]=="level3"){
        n3 <- length(unique(level3))
        pi <- rep(1 / n3, N)
        pij <- rep(rep(1 / mi, mi), kij) * pi
        pijk <-  pij / (unname(rep(kij, kij)))
      }
    else if(weights[1]=="level2"){
        pi <- rep(mi / sum(mi), ki)
        pij <- rep(1 / sum(mi), N)
        pijk <-  pij / (unname(rep(kij, kij)))
      }
    else if(weights[1]=="level1"){
        pi <- rep(ki / N, ki)
        pijk <-  rep(1 / N, N)
        pij <- rep(tapply(pijk, level32, sum), kij)
      }
    else stop("a wrong weighting method name entered!")
    #CDFs of observations
    ef <- emp_CDF(x, pijk)
    #averaged CDF
    avg <- sum(ef * pijk)
    #total variance
    tv <- sum((ef - avg)^2 * pijk)
    cl <- unique(level3)
    n <- length(cl)
    output <- matrix(NA, ncol = 4, nrow = 2)
    colnames(output) <- c("rankICC", "SE", "Lower", "Upper")
    rownames(output) <- c("gamma2", "gamma3")
    #################gamma2
    #covariance
    l_l <- tapply(ef - avg, level32, I)
    l_pij <- as.list(tapply(pij, level32, unique))
    cv <- mapply(cov_CDF, l_l, l_pij, USE.NAMES = FALSE)
    l_idx <- unlist(lapply(tapply(level3, level32, I), unique))
    cv <- tapply(cv, l_idx, sum)
    #estimate
    output["gamma2", "rankICC"] <- sum(cv) / tv
    #the standard error and a confidence interval
    if(conf.int){
      an <- sum(cv) / n
      bn <- tv / n
      d1 <- cv / bn
      d2 <- c(unname(tapply((ef - avg) ^ 2 * pijk, level3, sum)) * -an / bn ^ 2)
      dat <- data.frame(x = x, level3 = level3, level2 = level2, level32 = level32, pij = pij, pijk = pijk, ef = ef)
      d3 <- vapply(cl, d3f.gamma2, numeric(1), dat, avg, n, an, bn)
      se2 <- sd((d1 + d2 + d3) / sqrt(n))
      output["gamma2", c("Lower", "Upper")] <- getCI(output["gamma2", "rankICC"], se2, conf.int, fisher)
      output["gamma2", "SE"] <- se2
    }
    ##################gamma3
    #covariance
    l_l <- tapply(ef - avg, level3, I)
    l_pi <- as.list(tapply(pi, level3, unique))
    l_ij <- tapply(level2, level3, I)
    cv <- mapply(cov_CDF_multi, l_l, l_pi, l_ij, USE.NAMES = FALSE)
    #estimate
    output["gamma3", "rankICC"] <- sum(cv) / tv
    #the standard error and a confidence interval
    if(conf.int){
      an <- sum(cv) / n
      bn <- tv / n
      d1 <- cv / bn
      d2 <- c(unname(tapply((ef - avg) ^ 2 * pijk, level3, sum)) * -an / bn ^ 2)
      dat <- data.frame(x = x, level3 = level3, level2 = level2, pi = pi, pijk = pijk, ef = ef)
      d3 <- vapply(cl, d3f.gamma3, numeric(1), dat, cl, avg, n, an, bn)
      se3 <- sd((d1 + d2 + d3) / sqrt(n))
      output["gamma3", c("Lower", "Upper")] <- getCI(output["gamma3", "rankICC"], se3, conf.int, fisher)
      output["gamma3", "SE"] <- se3
    }
    return(output)
}

#Functions used for obtaining estimates for three-level data
cov_CDF_multi <- function(ef, pci, idx){
  s <- 0
  for(i in seq_along(idx)){
    s <- s + sum(ef[i] * ef[idx != idx[i]])
  }
  kij <- table(idx)
  wi <- (sum(kij))^2 - sum(kij^2)
  s <- s * pci / wi
  return(s)
}

d3f.gamma2 <- function(ip, dat, avg, n, an, bn) {
  l_rowix <- tapply(seq(nrow(dat)), dat[,'level3'], I)
  l_rowijx <- tapply(seq(nrow(dat)), dat[,'level32'], I)
  allL1 <- dat[,'ef'] - avg
  allNIJ <- lengths(l_rowijx)
  allPIJ <- tapply(dat[,'pij'], dat[,'level32'], unique)
  xi <- dat[l_rowix[[ip]],'x']
  pijk <- dat[l_rowix[[ip]],'pijk']
  allLTE <- outer_lte(xi, dat[,'x'])
  allCMP <- colSums(((allLTE[[1]] | allLTE[[2]]) + allLTE[[2]]) / 2 * pijk)
  t1 <- numeric(length(unique(dat[,'level32'])))
  dd1 <- 0
  for(ij in seq_along(t1)) {
    ijx <- l_rowijx[[ij]]
    l1 <- allL1[ijx]
    l2 <- allCMP[ijx]
    nij <- unname(allNIJ[ij])
    Pij <- allPIJ[[ij]]

    s1 <- sum(unlist(lapply(l1, `+`, l1)))
    s2 <- sum(l1 * 2)
    t1[ij] <- (s1 - s2) / (nij * (nij - 1)) * Pij

    s1 <- sum(unlist(lapply(l1, `*`, l2)))
    s2 <- sum(l1 * l2)
    dd1 <- dd1 + (s1 - s2) * 2 * Pij / (nij * (nij - 1))
  }

  fi <- dat[l_rowix[[ip]],'ef']
  t2 <- allLTE[[1]] / 2 + allLTE[[2]]
  calc1 <- colSums(t2 * pijk) * dat[,'pijk']
  cni <- sum(pijk * fi) + sum(calc1) / n
  t3 <- sum(calc1 * allL1) * 2
  t4 <- cni * sum(dat[,'pijk'] * allL1) * 2
  dd1 / bn - sum(t1) * cni / bn - an / (bn^2) * (t3 - t4)
}

d3f.gamma3 <- function(ip, dat, cl, avg, n, an, bn) {
  l_rowix <- tapply(seq(nrow(dat)), dat[,'level3'], I)
  allL1 <- dat[,'ef'] - avg
  allNI <- lengths(l_rowix)
  allPI <- tapply(dat[,'pi'], dat[,'level3'], unique)
  xi <- dat[l_rowix[[ip]],'x']
  pijk <- dat[l_rowix[[ip]],'pijk']
  allLTE <- outer_lte(xi, dat[,'x'])
  allCMP <- colSums(((allLTE[[1]] | allLTE[[2]]) + allLTE[[2]]) / 2 * pijk)
  dd1 <- t1 <- 0
  for(i in seq_along(cl)){
    ix <- l_rowix[[i]]
    l1 <- allL1[ix]
    l2 <- allCMP[ix]
    Pi <- allPI[[i]]
    cll2 <- dat[ix, 'level2']
    mij <- table(cll2)
    wi <- (sum(mij))^2 - sum(mij^2)
    s1 <- s2 <- 0
    for(j in unique(cll2)){
      l1j <- l1[cll2 == j]
      l2j <- l2[cll2 != j]
      s1 <- s1 + sum(unlist(lapply(l1j, `*`, l2j)))

      l1j <- l1[cll2 != j]
      l2j <- l2[cll2 == j]
      s1 <- s1 + sum(unlist(lapply(l1j, `*`, l2j)))

      s2 <- s2 + sum(unlist(lapply(l1[cll2 == j], `+`, l1[cll2 != j])))
    }
    dd1 <- dd1 + s1 * Pi / wi
    t1 <- t1 + s2 * Pi / wi
  }
  fi <- dat[l_rowix[[ip]],'ef']
  t2 <- allLTE[[1]] / 2 + allLTE[[2]]
  calc1 <- colSums(t2 * pijk) * dat[,'pijk']
  cni <- sum(pijk * fi) + sum(calc1) / n
  t3 <- sum(calc1 * allL1) * 2
  t4 <- cni * sum(dat[,'pijk'] * allL1) * 2
  dd1 / bn - t1 * cni / bn - an / (bn^2) * (t3 - t4)
}
