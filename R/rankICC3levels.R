#' Rank ICC with three hierarchies
#'
#' \code{rankICC3levels} computes the rank intraclass correlation coefficient (ICC) with three hierarchies. Starting from the innermost level, the three levels are named level 1, level 2, and level 3. The rank ICC at level 2 evaluates the rank correlation between a random pair from the same level-2 unit. The rank ICC at level 3 evaluates the rank correlation between a random pair from the same level-3 unit but different level-2 units.
#'
#' @param x a numeric or factor vector.
#' @param level2 a vector indicating level-2 membership.
#' @param level3 a vector indicating level-3 membership.
#' @param weights a character string indicating which weighting method is used. Or an optional vector of user-defined weights to be used. Should be one of the strings \code{"level1"}, \code{"level2"}, \code{"level3"}, or a numeric vector. Default is \code{"level1"}. See Details.
#' @param conf.int numeric specifying confidence interval level.
#' @param fisher logical, indicating whether to apply Fisher transformation to compute confidence intervals.
#' @param na.rm logical. Should missing values be removed?
#' @details \code{"level1"} assigns equal weights to level-1 units; \eqn{p_{ijk}=1/(\sum_{i=1}^n\sum_{j=1}^{n_i}m_{ij})}, where \eqn{n} is the total number of level-3 units, \eqn{n_i} is the number of level-2 units in the \eqn{i}th level-3 unit, and \eqn{m_{ij}} is the number of level-1 units in the \eqn{j}th level-2 unit and the \eqn{i}th level-3 unit. \code{"level2"} assigns equal weights to level-2 units; \eqn{p_{ijk}=1/(m_{ij}\sum_{i=1}^n n_{i})}. \code{"level3"} assigns equal weights to level-3 units; \eqn{p_{ijk}=1/(nn_{i}m_{ij})}.
#' @return a matrix with two rows. The first row is for rank ICC at level 2 and the second row is for rank ICC at level 3. Each row has the following components.
#' \tabular{ll}{
#'   \code{rankICC} \tab the rank ICC. \cr
#'   \tab \cr
#'   \code{SE} \tab the standard error. \cr
#'   \tab \cr
#'   \code{Lower, Upper} \tab the lower and upper bound of the confidence interval.\cr
#' }
#' @references Tu, S., Li, C., Zeng, D., and Shepherd, B. E. (2023). Rank intraclass correlation for clustered data. Statistics in Medicine 42, 4333-4348.
#' @examples
#' k <- 50; m <- 5
#' sigma.u <- 1; sigma.e <- 2
#' u <- rnorm(k, 5, sigma.u)
#' x1 <- matrix(NA, k, m)
#' for (i in 1:k){
#' x1[i,] <- u[i] + rnorm(5, 0, sigma.e)
#' }
#' x <- as.vector(t(x1))
#' level2 <- rep(1:k, each=5)
#' level3 <- round(level2 / 10)
#' rankICC3levels(x, level2, level3, weights = "level2")
#' @export
#' @importFrom stats complete.cases qnorm sd

rankICC3levels <- function(x, level2, level3, weights = c("level1", "level2", "level3"),
                        conf.int = 0.95, fisher = FALSE, na.rm = FALSE){
    if(!is.numeric(x) & !is.factor(x)) stop("x must be a numeric or factor vector")
    else if(length(level2) != length(x) | length(level3) != length(x)) stop("lengths of x and level membership differ")
    else if(is.numeric(weights) & length(weights) != length(x)) stop("lengths of x and user-defined weights differ")
    if(!is.numeric(weights) & !(weights[1] %in% c("level1", "level2", "level3"))) stop("a wrong method name entered!")
    idx_order <- seq_along(x)
    if(na.rm){
      idx_order <- complete.cases(x, level2, level3)
      x <- x[idx_order]
      level2 <- level2[idx_order]
      level3 <- level3[idx_order]
      if(is.numeric(weights)) weights <- weights[idx_order]
    }
    if(is.factor(x)) x <- as.numeric(x)
    g2 <- rankICC3levels_gamma2(x, level2, level3, weights = weights, conf.int = conf.int, fisher = fisher)
    g3 <- rankICC3levels_gamma3(x, level2, level3, weights = weights, conf.int = conf.int, fisher = fisher)
    output <- rbind(g2, g3)
    return(output)
}

#Functions used for obtaining estimates for three-level data
rankICC3levels_gamma2 <- function(x, level2, level3, weights = c("level1", "level2", "level3"),
                          conf.int = 0.95, fisher = FALSE){
  #remove clusters with one observation
  level32 <- paste(level3, level2, sep="-")
  kij <- table(level32)
  if(sum(kij==1)){
    idx <- names(kij[kij > 1])
    level2 <- level2[level32 %in% idx]
    level3 <- level3[level32 %in% idx]
    x <- x[level32 %in% idx]
    warning("level-2 units with only one observation were removed")
  }
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
    pijk <- weights[idx]
    pij <- rep(tapply(pijk, level32, sum), kij)
    pi <- rep(tapply(pijk, level3, sum), ki)
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
  #CDFs of observations
  ef <- emp_CDF(x, pijk)
  #averaged CDF
  avg <- sum(ef * pijk)
  #total variance
  tv <- sum((ef - avg)^2 * pijk)
  cl <- unique(level3)
  n <- length(cl)
  output <- matrix(NA, ncol = 4, nrow = 1)
  colnames(output) <- c("rankICC", "SE", "Lower", "Upper")
  rownames(output) <- c("gamma2")
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
    d3 <- vapply(cl, d3f_gamma2, numeric(1), dat, avg, n, an, bn)
    se2 <- sd((d1 + d2 + d3) / sqrt(n))
    output["gamma2", c("Lower", "Upper")] <- getCI(output["gamma2", "rankICC"], se2, conf.int, fisher)
    output["gamma2", "SE"] <- se2
  }
  return(output)
}


rankICC3levels_gamma3 <- function(x, level2, level3, weights = c("level1", "level2", "level3"),
                                  conf.int = 0.95, fisher = FALSE){
  #remove clusters with one observation
  level3 <- as.character(level3)
  mi <- lapply(tapply(level2, level3, unique), length)
  if(sum(mi==1)){
    idx <- names(mi[mi > 1])
    x <- x[level3 %in% idx]
    level2 <- level2[level3 %in% idx]
    level3 <- level3[level3 %in% idx]
    warning("level-3 units with only one level-2 unit were removed")
  }
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
    pijk <- weights[idx]
    pij <- rep(tapply(pijk, level32, sum), kij)
    pi <- rep(tapply(pijk, level3, sum), ki)
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
  #CDFs of observations
  ef <- emp_CDF(x, pijk)
  #averaged CDF
  avg <- sum(ef * pijk)
  #total variance
  tv <- sum((ef - avg)^2 * pijk)
  cl <- unique(level3)
  n <- length(cl)
  output <- matrix(NA, ncol = 4, nrow = 1)
  colnames(output) <- c("rankICC", "SE", "Lower", "Upper")
  rownames(output) <- c("gamma3")
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
    d3 <- vapply(cl, d3f_gamma3, numeric(1), dat, cl, avg, n, an, bn)
    se3 <- sd((d1 + d2 + d3) / sqrt(n))
    output["gamma3", c("Lower", "Upper")] <- getCI(output["gamma3", "rankICC"], se3, conf.int, fisher)
    output["gamma3", "SE"] <- se3
  }
  return(output)
}

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

d3f_gamma2 <- function(ip, dat, avg, n, an, bn) {
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

d3f_gamma3 <- function(ip, dat, cl, avg, n, an, bn) {
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
