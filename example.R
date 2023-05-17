########The function to generate two-level clustered data
#n is the number of clusters
#cluster.limit is the range of cluster sizes. The cluster size follows a uniform distribution within the range.
#mu is the mean 
#vu is the between-cluster variance 
#vr is the within-cluster variance
generate.data <- function(seed=1234, n=30, cluster.limit=c(10,10), mu=1, vu=1,vr=1){
  set.seed(seed)
  #cluster size
  if(cluster.limit[1] == cluster.limit[2]){
      size.cluster <- rep(cluster.limit[1], n)
  }else size.cluster <- replicate(n, sample(cluster.limit[1]:cluster.limit[2], 1))
  #generate data for each cluster
  x <- list()
  cluster <- list()
  id <- list()
  #generate cluster means
  u <- rnorm(n, mean = mu, sd = sqrt(vu))
  for(i in 1:n){
    n.cluster <- size.cluster[i]
    r <- rnorm(n.cluster, mean = 0, sd = sqrt(vr))
    x[[i]] <- r+u[i] 
    cluster[[i]] <- rep(i, n.cluster)
    id[[i]] <- 1:n.cluster
  }
  dat <- data.frame("x"= unlist(x),
                    "cluster"=as.factor(unlist(cluster)),
                    "id"=as.factor(unlist(id)))
  return(dat)
}

#Examples of using the function in the rankICC package to obtain the rank ICC for two-level data 
dat <- generate.data(1234, n=30, cluster.limit=c(2,10), mu=1, vu=1, vr=2)
#Assigning equal weights to clusters
rankICC::rankICC(dat$x, dat$cluster, weights = "clusters")
#Assigning equal weights to observations 
rankICC::rankICC(dat$x, dat$cluster, weights = "obs")
#Iterative weighting based on the effective sample size
rankICC::rankICC(dat$x, dat$cluster, weights = "ess")
#Iterative weighting based on the combination of equal weights for clusters and equal weights for observations
rankICC::rankICC(dat$x, dat$cluster, weights = "combination")

########The function to generate three-level clustered data
#n1 is the range of the number of level-1 units. The number of level-1 units follows a uniform distribution within the range.
#n2 is the range of the number of level-2 units. The number of level-2 units follows a uniform distribution within the range.
#n3 is the number of level-3 units
#mu is the mean 
#v1 is the variance within a level-2 unit
#v2 is the variance between level-2 units within a level-3 unit
#v3 is the variance between level-3 units
generate.data.multi <- function(seed=1234, n1=c(10,10), n2=c(10,10), n3=30, mu=1, v1=1, v2=1, v3=1){
  set.seed(seed)
  #cluster size
  if(n2[1] == n2[2]){
    size <- rep(n2[1], n3)
  }else size <- replicate(n3, sample(n2[1]:n2[2], 1))
  u <- rnorm(n3, mean = mu, sd = sqrt(v3))
  dat <- list()
  for(i in 1:n3){
    dsub <- generate.data(n=size[i], cluster.limit=n1, mu=u[i], vu=v2, vr=v1)
    colnames(dsub)[2] <- "level2"
    dsub$level3 <- i
    dat[[i]] <- dsub
  }
  dat <- do.call(rbind, dat)
  dat$level3 <- as.factor(dat$level3)
  return(dat)
}
#Examples of using the function to obtain the rank ICC for three-level data 
dat <- generate.data.multi(seed=124, n1=c(2,10), n2=c(5,15), n3=30, mu=1, v1=1, v2=2, v3=2)
#Equal weights for level-1 units 
rankICC::rankICC3levels(dat$x, dat$level2, dat$level3, weights = "level1")
#Equal weights for level-2 units 
rankICC::rankICC3levels(dat$x, dat$level2, dat$level3, weights = "level2")
#Equal weights for level-3 units 
rankICC::rankICC3levels(dat$x, dat$level2, dat$level3, weights = "level3")

