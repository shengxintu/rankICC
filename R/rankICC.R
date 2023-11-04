#' Rank ICC with two hierarchies
#'
#' \code{rankICC} computes the rank intraclass correlation coefficient (ICC) of a two-level hierarchical distribution. It can be used with any orderable variable, including continuous and discrete variables. Different weighting methods are provided, including methods assigning equal weights to observations or to clusters.
#' @param x a numeric or factor vector.
#' @param cluster a vector of cluster index corresponding to \code{x}.
#' @param weights a character string indicating which weighting method is used. Or an optional vector of user-defined weights to be used. Should be one of the strings \code{"obs"}, \code{"clusters"}, \code{"ess"}, \code{"combination"}, or a numeric vector. Default is \code{"obs"}. See Details.
#' @param conf.int numeric specifying confidence interval level.
#' @param fisher logical, indicating whether to apply Fisher transformation to compute confidence intervals.
#' @param na.rm logical. Should missing values be removed?
#' @param ... additional arguments to be passed to the iteration function if \code{weights} is \code{"ess"} or \code{"combination"}. Specifying the tolerance via \code{"tol"} and the maximum iteration times via \code{"maxIter"}.
#' @details \code{"obs"} assigns equal weights to observations; \eqn{p_{ij} = 1/N}, where \var{N} is the total number of observations. \code{"clusters"} assigns equal weights to clusters; \eqn{p_{ij} = 1/(nk_i)}, where \var{n} is the total number of clusters and k_i is the cluster size. \code{"ess"} and \code{"combination"} implement iterations until convergence; \eqn{p_{ij}(\gamma_I)=1/(1+k_i\gamma_I)/\{\sum_{j=1}^n k_j/(1+k_j\gamma_I)\}} for \code{"ess"}, \eqn{p_{ij}(\gamma_I)=(1-\gamma_I)/N+\gamma_I/(nk_i)} for \code{"combination"}.
#' @return a vector with following components.
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
#' cluster <- rep(1:k, each=5)
#' rankICC(x, cluster, weights = "clusters")
#' rankICC(x, cluster, weights = "ess", tol = 1e-4, maxIter = 10)
#' @export
#' @importFrom stats complete.cases qnorm sd

rankICC <- function(x, cluster, weights = c("obs", "clusters", "ess", "combination"),
                           conf.int = 0.95, fisher = FALSE, na.rm = FALSE, ...){
  if(!is.numeric(x) & !is.factor(x)) stop("x must be a numeric or factor vector")
  else if(is.numeric(weights) & length(weights) != length(x)) stop("lengths of x and user-defined weights differ")
  else if(length(cluster) != length(x)) stop("lengths of x and cluster differ")
  idx <- seq_along(x)
  if(na.rm){
    idx <- complete.cases(x, cluster)
    x <- x[idx]
    cluster <- cluster[idx]
  }
  if(is.factor(x)) x <- as.numeric(x)
  if(is.numeric(weights)){
    weights <- weights[idx]
    output <- rankICCest(x, cluster, user_defined_weights = weights,
                                conf.int = conf.int, fisher = fisher)
  }
  else if(weights[1] %in% c("ess", "combination")) output <- rankICCest_iter(x, cluster, weights = weights[1],
                                                                               conf.int = conf.int, fisher = fisher, ...)
  else if(weights[1] %in% c("obs", "clusters")){
      ri <- ifelse(weights[1] == "obs", 0, 1)
      output <- rankICCest(x, cluster, ri = ri, conf.int = conf.int, fisher = fisher)
  }
  else stop("a wrong weighting method name entered!")
  return(output)
}

rankICCest <- function(x, cluster, ri = 0, opt_method = "ess",
                       user_defined_weights = NULL, conf.int = 0.95, fisher = FALSE){
  #Check whether there are clusters with only one observations
  cluster <- as.character(cluster)
  ki <- table(cluster)
  if(sum(ki==1)){
    idx <- names(ki[ki > 1])
    x <- x[cluster %in% idx]
    cluster <- cluster[cluster %in% idx]
    warning("clusters with only one observation were removed")
  }
  cluster <- factor(cluster, levels=unique(cluster))
  x <- x[order(cluster)]
  #set up the weights
  if(is.null(user_defined_weights)){
    cluster <- sort(cluster)
    ki <- tabulate(cluster)
    if(opt_method == "ess"){
      #effective sample size
      neffi <- ki / (1 + (ki - 1) * ri)
      #cluster weights
      Pi <- neffi / sum(neffi)
      #individual weights
      pij <- rep(Pi / ki, ki)
      pci <- rep(Pi, ki)
    }
    if(opt_method=="combination"){
      N <- length(x)
      n <- length(ki)
      size <- rep(ki, ki)
      pij <- (1-ri)/N + ri/size/n
      Pi <- (1-ri)/N * ki + ri/n
      pci <- pij * size
    }
  }
  else{
    user_defined_weights <- user_defined_weights/sum(user_defined_weights)
    pij <- user_defined_weights[order(cluster)]
    cluster <- sort(cluster)
    ki <- tabulate(cluster)
    Pi <- tapply(pij, cluster, sum)
    pci <- rep(Pi, ki)
  }
  #CDFs of observations
  ef <- emp_CDF(x, pij)
  dat <- data.frame(x = x, cluster = cluster, pci = pci, pij = pij, ef = ef)
  #averaged CDF
  avg <- sum(ef * pij)
  #total variance
  tv <- sum((ef - avg)^2 * pij)
  cl <- unique(cluster)
  #covariance
  l_l <- tapply(ef - avg, cluster, I)
  l_p <- as.list(Pi)
  cv <- mapply(cov_CDF, l_l, l_p, USE.NAMES = FALSE)
  #estimate
  est <- sum(cv) / tv
  output <- c(est, rep(NA, 3))
  names(output) <- c("rankICC", "SE", "Lower", "Upper")
  ####calculate the standard error and a confidence interval
  if(conf.int){
    dat <- data.frame(x = x, cluster = cluster, pci = pci, pij = pij, ef = ef)
    n <- length(cl)
    an <- sum(cv) / n
    bn <- tv / n
    d1 <- cv / bn
    d2 <- c(unname(tapply((ef - avg) ^ 2 * pij, cluster, sum)) * -an / bn ^ 2)
    d3 <- vapply(cl, d3f, numeric(1), dat, cl, avg, n, an, bn)
    se <- sd((d1 + d2 + d3) / sqrt(n))
    output[c("Lower", "Upper")] <- getCI(est, se, conf.int, fisher)
    output["SE"] <- se
  }
  return(output)
}

rankICCest_iter <- function(x, cluster, weights = c("ess", "combination"),
                            conf.int = 0.95, fisher = FALSE, tol = 1e-5, maxIter = 100){
  i <- 0; d <- 10; ri0 <- 0
  while(i < maxIter & d > tol){
      rnew <- rankICCest(x, cluster, ri = ri0, opt_method = weights[1], conf.int = 0)
      d <- abs(rnew["rankICC"] - ri0)
      ri0 <- rnew["rankICC"]
      i <- i + 1
  }
  if(d > tol) warning("algorithm did not converge!")
  output <- rankICCest(x, cluster, ri0, opt_method = weights[1], conf.int = conf.int, fisher = fisher)
  return(output)
}

emp_CDF <- function(x, pij, tol = 1e-7) {
  vapply(x, function(i) {
    c1 <- x - i
    isEq <- abs(c1) < tol
    isLT <- (!isEq & c1 < 0)
    sum(isLT * pij + (isEq * pij) / 2)
  }, numeric(1))
}

cov_CDF <- function(ef, pci){
  mat <- outer(ef, ef, `*`)
  s <- sum(mat) - sum(diag(mat))
  r <- s / (nrow(mat) * (nrow(mat)-1)) * pci
  return(r)
}

outer_lte <- function(a, b, tol = 1e-7) {
  m1 <- outer(a, b, `-`)
  isEq <- abs(m1) < tol
  isLT <- (!isEq & m1 < 0)
  list(isEq, isLT)
}

d3f <- function(ip, dat, cl, avg, n, an, bn) {
  l_rowix <- tapply(seq(nrow(dat)), dat[,'cluster'], I)
  allL1 <- dat[,'ef'] - avg
  allNI <- lengths(l_rowix)
  allPI <- tapply(dat[,'pci'], dat[,'cluster'], unique)
  xi <- dat[l_rowix[[ip]],'x']
  pij <- dat[l_rowix[[ip]],'pij']
  allLTE <- outer_lte(xi, dat[,'x'])
  allCMP <- colSums(((allLTE[[1]] | allLTE[[2]]) + allLTE[[2]]) / 2 * pij)#F(xi'j|xij'')
  t1 <- numeric(length(cl))
  d34 <- 0
  for(i in seq_along(t1)) {
    ix <- l_rowix[[i]]
    l1 <- allL1[ix]
    l2 <- allCMP[ix]
    ni <- unname(allNI[i])
    pi <- allPI[[i]]
    s1 <- sum(unlist(lapply(l1, `+`, l1)))
    s2 <- sum(l1 * 2)
    t1[i] <- (s1 - s2) / (ni * (ni - 1)) * pi
    s1 <- sum(unlist(lapply(l1, `*`, l2)))
    s2 <- sum(l1 * l2)
    d34 <- d34 + (s1 - s2) * 2 * pi / (ni * (ni - 1))
  }
  fi <- dat[l_rowix[[ip]],'ef']
  t2 <- allLTE[[1]] / 2 + allLTE[[2]]
  calc1 <- colSums(t2 * pij) * dat[,'pij']
  cni <- sum(pij * fi) + sum(calc1) / n
  t3 <- sum(calc1 * allL1) * 2
  t4 <- cni * sum(dat[,'pij'] * allL1) * 2
  d34 / bn - sum(t1) * cni / bn - an / (bn^2) * (t3 - t4)
}

getCI <- function(x, s, conf.int = 0.95, fisher = FALSE){
  if(fisher){
    #Fisher transformation
    tr <- log((1 + x)/(1 - x)) / 2
    ts <- s / ((1 + x)*(1 - x))
    l <- tr - qnorm(1 / 2 + conf.int / 2) * ts
    u <- tr + qnorm(1 / 2 + conf.int / 2) * ts
    l <- (exp(2 * l) - 1) / (exp(2 * l) + 1)#tanh()
    u <- (exp(2 * u) - 1) / (exp(2 * u) + 1)
  }
  else{
    l <- x - qnorm(1 / 2 + conf.int / 2) * s
    u <- x + qnorm(1 / 2 + conf.int / 2) * s
  }
  return(c(l, u))
}

