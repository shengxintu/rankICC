#############################Codes for applications
##############Patient health questionnaire-9 
load("data/baseline_hops_dep_deid.rda")
df <- baseline_hops_dep_deid[complete.cases(baseline_hops_dep_deid),]
p <- summarise(group_by(df,Female,Male),length(Male))
colnames(p)[3] <- "Count"
#Figure 10
ggplot(p,aes(Female,Male))+ geom_point(aes(size=Count), shape=20)+ theme_classic()+
  theme(legend.margin = margin(-5,10, -5, 10))+ theme(text = element_text(size=15))+
  scale_size(breaks = c(1, 5, 20, 150))+xlab("PHQ-9 score (female)") + ylab("PHQ-9 score (male)") +
  scale_y_continuous(limits = c(0, 27), breaks = seq(0, 25, 5)) +
  scale_x_continuous(limits = c(0, 27), breaks = seq(0, 25, 5))
#Rank ICC at the couple and clinical site levels 
dt <- data.frame(x=c(df$Female,df$Male), 
                 level2=factor(rep(1:nrow(df), 2)), 
                 level3 = factor(rep(df$cluster, 2)), 
                 sex=rep(c("Female","Male"), each=nrow(df)))
dt <- dt %>% arrange(level3, level2)
rankICC::rankICC3levels(dt$x, dt$level2, dt$level3, weights = "level2")
#Rank ICC of females at the clinical level
rankICC::rankICC(dt[dt$sex=="Female", ]$x, dt[dt$sex=="Female", ]$level3)
#Rank ICC of males at the clinical level
rankICC::rankICC(dt[dt$sex=="Male", ]$x, dt[dt$sex=="Male", ]$level3)
#ICC estimation
dt <- dt %>% arrange(level3, level2)
m <- lmer(x ~ (1|level2)+(1|level3), data = dt, REML=T)
v <- as.data.frame(VarCorr(m))[,"vcov"]
#ICC at the couple level
(v[1]+v[2])/sum(v)
#ICC at the clinical site level 
v[2]/sum(v)
#ICC of females at the clinical level
m.f <- lmer(x ~ (1|level3), data = dt[dt$sex == "Female",], REML=T)
v <- as.data.frame(VarCorr(m.f))[,"vcov"]
v[1]/sum(v)
#ICC of males at the clinical level
m.m <- lmer(x ~ (1|level3), data = dt[dt$sex == "Male",], REML=T)
v <- as.data.frame(VarCorr(m.m))[,"vcov"]
v[1]/sum(v)


##############Albumin:Creatinine Ratio 
Final_Nigeria_R3_Data <- readxl::read_excel("data/Final_Nigeria_R3_Data.xlsx")
#Figure 8
ggplot(Final_Nigeria_R3_Data %>% filter(Dolutegravir == 1),
       aes(ACR_1, ACR_2))+ geom_point() + theme_classic()+ theme(text = element_text(size=15))+
  xlab("First uACR (mg/g)") + ylab("Second uACR (mg/g)")+
  scale_y_continuous(limits = c(0, 3000), breaks = seq(0,3000, 1000)) +
  scale_x_continuous(limits = c(0, 3000), breaks = seq(0,3000, 1000))
dat <- Final_Nigeria_R3_Data %>% 
  gather("ACR","value",c("ACR_1","ACR_2"))%>% 
  mutate(cluster = factor(rep(1:nrow(Final_Nigeria_R3_Data), 2))) %>% 
  filter(Dolutegravir == 1) %>% 
  arrange(cluster)
#Rank ICC of original data 
rankICC::rankICC(dat$value, dat$cluster)
#ICC of original data 
m <- lmer(x ~ (1|cluster), data = dat, REML = T)
v <- as.data.frame(VarCorr(m))[,"vcov"]
v[1]/sum(v)
#Rank ICC after extreme values removed 
idx <- dat %>% filter(value > 2000 & ACR == "ACR_2") %>% pull(cluster)
x <- dat %>% filter(cluster != idx) %>% pull(value)
cluster <- dat %>% filter(cluster != idx) %>% pull(cluster)
rankICC::rankICC(x, cluster)
#ICC after extreme values removed 
m <- lmer(value ~ (1|cluster), data = dat %>% filter(cluster != idx), REML = T)
v <- as.data.frame(VarCorr(m))[,"vcov"]
v[1]/sum(v)
#Rank ICC after log transformation
rankICC::rankICC(log(dat$value), dat$cluster)
#ICC after log-transformation 
m <- lmer(value ~ (1|cluster), data = dat %>% mutate(value = log(value)), REML = T)
v <- as.data.frame(VarCorr(m))[,"vcov"]
v[1]/sum(v)
#Rank ICC after square root transformation
rankICC::rankICC(sqrt(dat$value), dat$cluster)
#ICC after square root transformation 
m <- lmer(value ~ (1|cluster), data = dat %>% mutate(value = sqrt(value)), REML = T)
v <- as.data.frame(VarCorr(m))[,"vcov"]
v[1]/sum(v)

##############Status Epilepticus
SeizureTrack <- readxl::read_excel("data/SeizureTrack.xlsx",sheet = "Seiz prior Raw")
SeizureTrack <- SeizureTrack[!is.na(SeizureTrack$seize_prior_bm1),]
SeizureTrack$cluster <- stringr::str_extract(SeizureTrack$record_id,"[0-9]+")
#Figure 8
hist(SeizureTrack$seize_prior_bm1, breaks = 50, xlab="Numbers of seizures", main="", cex.lab=1.3)
dat <- data.frame(x= SeizureTrack$seize_prior_bm1, cluster = factor(SeizureTrack$cluster)) %>% arrange(cluster)
#Rank ICC based on assigning equal weights to clusters
rankICC::rankICC(dat$x, dat$cluster, weights = "clusters")
#Rank ICC based on assigning equal weights to observations
rankICC::rankICC(dat$x, dat$cluster, weights = "obs")
#Rank ICC based on iterative weighting of the effective sample sizes
rankICC::rankICC(dat$x, dat$cluster, weights = "ess")
#Rank ICC based on iterative weighting of the linear combination of equal clusters and equal observations
rankICC::rankICC(dat$x, dat$cluster, weights = "combination")
#ICC from generalized random effects model 
#linear
dat <- dat %>% group_by(cluster) %>% mutate(ki = n())
m <- lmer(x ~ (1|cluster), data = dat, REML=T)
v <- as.data.frame(VarCorr(m))[,"vcov"]
v[1]/sum(v)
#quasipoisson link
mq <- MASS::glmmPQL(x ~ 1, random = list(~1|cluster), family = "quasipoisson", data = dat)
omegaN <- as.numeric(VarCorr(mq)[2, 1])
lam2 <- exp(fixef(mq) + 0.5 * as.numeric(VarCorr(mq)[1, 1]))
VarOtN <- trigamma(lam2/omegaN)
#trigamma (usually preferred)
as.numeric(VarCorr(mq)[1, 1])/(as.numeric(VarCorr(mq)[1, 1]) + VarOtN)
#lognormal approximation
as.numeric(VarCorr(mq)[1, 1])/(as.numeric(VarCorr(mq)[1, 1]) + log(1+omegaN/lam2))
#delta method
as.numeric(VarCorr(mq)[1, 1])/(as.numeric(VarCorr(mq)[1, 1]) + omegaN/lam2)
#negative binomial link
mnb <- glmer.nb(x ~ 1 + (1 | cluster), data = dat)
thetaN <- getME(mnb, "glmer.nb.theta")
lambda <- exp(fixef(mnb) + 0.5 * as.numeric(VarCorr(mnb)))
VarOtN <- trigamma((1/lambda + 1/thetaN)^(-1)) # trigamma function
as.numeric(VarCorr(mnb))/(as.numeric(VarCorr(mnb)) + VarOtN)



#############################Simulation codes
#######Functions to generate clustered data
generate.data.s1 <- function(seed=1234, n=30, cluster.limit=c(10,10), mu=1, vu=1,vr=1, hlf=F){
  set.seed(seed)
  #cluster size
  if(hlf) size.cluster <- rep(cluster.limit, each = n/2)
  else{
    if(cluster.limit[1] == cluster.limit[2]){
      size.cluster <- rep(cluster.limit[1], n)
    }else size.cluster <- replicate(n, sample(cluster.limit[1]:cluster.limit[2], 1))
  }
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

generate.data.s3 <- function(seed=1234, n=30, cluster.limit=c(10,10), mu=1, vu=1,vr=1, hlf=F){
  set.seed(seed)
  #cluster size
  if(hlf) size.cluster <- rep(cluster.limit, each = n/2)
  else{
    if(cluster.limit[1] == cluster.limit[2]){
      size.cluster <- rep(cluster.limit[1], n)
    }else size.cluster <- replicate(n, sample(cluster.limit[1]:cluster.limit[2], 1))
  }
  #generate data for each cluster
  x <- list()
  cluster <- list()
  id <- list()
  #generate cluster means
  u <- exp(rnorm(n, mean = mu, sd = sqrt(vu)))
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


generate.data.s4 <- function(seed=1234, n=30, cluster.limit=c(2,30), mu=1, vu=1, vr1=0.5, vr2=1.5){
  set.seed(seed)
  #cluster size
  size.cluster <- rep(cluster.limit, each = n/2)
  #generate data for each cluster
  x <- list()
  cluster <- list()
  id <- list()
  #generate cluster means
  u <- rnorm(n, mean = mu, sd = sqrt(vu))
  n1 <- n/2
  for(i in 1:n1){
    n.cluster <- size.cluster[i]
    r <- rnorm(n.cluster, mean = 0, sd = sqrt(vr1))
    x[[i]] <- r+u[i] 
    cluster[[i]] <- rep(i, n.cluster)
    id[[i]] <- 1:n.cluster
  }
  for(i in (n1+1):n){
    n.cluster <- size.cluster[i]
    r <- rnorm(n.cluster, mean = 0, sd = sqrt(vr2))
    x[[i]] <- r+u[i] 
    cluster[[i]] <- rep(i, n.cluster)
    id[[i]] <- 1:n.cluster
  }  
  dat <- data.frame("x"= unlist(x),
                    "cluster"=as.factor(unlist(cluster)),
                    "id"=as.factor(unlist(id)))
  return(dat)
}


generate.data.s5 <- function(seed=1234, n=30, cluster.limit=c(10,10), mu=1, vu=1,vr=1, hlf=F, l=5){
  set.seed(seed)
  #cluster size
  if(hlf) size.cluster <- rep(cluster.limit, each = n/2)
  else{
    if(cluster.limit[1] == cluster.limit[2]){
      size.cluster <- rep(cluster.limit[1], n)
    }else size.cluster <- replicate(n, sample(cluster.limit[1]:cluster.limit[2], 1))
  }
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
  ###discretize
  x <- unlist(x)
  bi <- c(-Inf, quantile(x, seq(1/l, 1-1/l, 1/l)), Inf)
  xi <- as.numeric(cut(x, breaks = bi, labels = 1:l))  
  dat <- data.frame("x"= xi,
                    "cluster"=as.factor(unlist(cluster)),
                    "id"=as.factor(unlist(id)))
  return(dat)
}


########A function to obtain the midrank 
CDF.mean <- function(x, pij, tol = 1e-7) {
  vapply(x, function(i) {
    c1 <- x - i
    isEq <- abs(c1) < tol
    isLT <- (!isEq & c1 < 0)
    sum(isLT * pij + (isEq * pij) / 2)
  }, numeric(1))
}

########A function to obtain covariance 
CDF.cov <- function(ef, pci){
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

#########Calculate components for standard error 
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
    Pi <- allPI[[i]]
    s1 <- sum(unlist(lapply(l1, `+`, l1)))
    s2 <- sum(l1 * 2)
    t1[i] <- (s1 - s2) / (ni * (ni - 1)) * Pi
    
    s1 <- sum(unlist(lapply(l1, `*`, l2)))
    s2 <- sum(l1 * l2)
    d34 <- d34 + (s1 - s2) * 2 * Pi / (ni * (ni - 1))
  }
  
  fi <- dat[l_rowix[[ip]],'ef']
  t2 <- allLTE[[1]] / 2 + allLTE[[2]]
  calc1 <- colSums(t2 * pij) * dat[,'pij']
  cni <- sum(pij * fi) + sum(calc1) / n
  t3 <- sum(calc1 * allL1) * 2 
  t4 <- cni * sum(dat[,'pij'] * allL1) * 2 
  d34 / bn - sum(t1) * cni / bn - an / (bn^2) * (t3 - t4)
}

#########A function for estimates and standard error 
est.Wcdf <- function(x, cluster, ri, std=T, opt_method='eff') {
  
  cluster <- factor(cluster, levels=unique(cluster))
  x <- x[order(cluster)]
  cluster <- sort(cluster)
  ####obtain estimates
  #cluster size
  ki <- tabulate(cluster) 
  if(opt_method == "eff"){
    #effective sample size
    neffi <- ki / (1 + (ki - 1) * ri) 
    #cluster weights
    Pi <- neffi / sum(neffi)
    #individual weights
    pij <- rep(Pi / ki, ki)  
    pci <- rep(Pi, ki)
  }
  if(opt_method=="combo"){
    N <- length(x)
    n <- length(ki)
    size <- rep(ki, ki)
    pij <- (1-ri)/N + ri/size/n
    Pi <- (1-ri)/N * ki + ri/n
    pci <- pij * size
  }
  #CDFs of observations
  ef <- CDF.mean(x, pij) 
  dat <- data.frame(x = x, cluster = cluster, pci = pci, pij = pij, ef = ef)
  #averaged CDF
  avg <- sum(ef * pij) 
  #total variance
  tv <- sum((ef - avg)^2 * pij)  
  cl <- unique(cluster)
  
  #covariance
  l_l <- tapply(ef - avg, cluster, I) 
  l_p <- as.list(Pi)
  cv <- mapply(CDF.cov, l_l, l_p, USE.NAMES = FALSE) 
  #estimate
  est <- sum(cv) / tv 
  ####calculate standard error 
  if(std){
    n <- length(cl)
    an <- sum(cv) / n
    bn <- tv / n
    d1 <- cv / bn
    d2 <- c(unname(tapply((ef - avg) ^ 2 * pij, cluster, sum)) * -an / bn ^ 2)
    d3 <- vapply(cl, d3f, numeric(1), dat, cl, avg, n, an, bn)
    se <- sd((d1 + d2 + d3) / sqrt(n))
    output <- list(est = est, se = se)
  }
  else output <- est
  
  return(output)
}


est.iter <- function(x, cluster, ri=0, tol=1e-5, maxIter=100, std=T){
  i <- 0; d <- 10
  ri1 <- ri2 <- ri
  while(i < maxIter & d > tol){
    rnew <- est.Wcdf(x, cluster, ri1, opt_method = "eff", std=F)
    d <- abs(rnew - ri1)
    ri1 <- rnew
    i <- i + 1
  }
  n1 <- i
  i <- 0; d <- 10
  while(i < maxIter & d > tol){
    rnew <- est.Wcdf(x, cluster, ri2, opt_method = "combo", std=F)
    d <- abs(rnew - ri2)
    ri2 <- rnew
    i <- i + 1
  }
  n2 <- i  
  if(std){
    est1 <- est.Wcdf(x, cluster, ri1, opt_method = "eff")
    est2 <- est.Wcdf(x, cluster, ri2, opt_method = "combo")
    output <- list("est"=c(ri1, ri2),
                   "se"=c(est1$se, est2$se),
                   "Niter"=c(n1, n2))
  }
  else{output <- list("est"=c(ri1, ri2),
                      "Niter"=c(n1, n2))}
  return(output)
}


est.Wcdf.boot <- function(seed, dat, N=100, ri=1, subboot=T){
  set.seed(seed)
  cl <- unique(dat$cluster)
  n <- length(cl)
  count <- table(dat$cluster)
  est <- matrix(NA, ncol = 2, nrow = N)
  for(i in 1:N){
    #sample clusters
    selectCl <- sample(cl, size = n, replace = T)
    df.bs <- sapply(selectCl, function(x) which(dat$cluster == x))
    df.bs <- dat[unlist(df.bs),]
    #rename selected clusters
    df.bs$cluster <- rep(paste(selectCl, 1:length(selectCl), sep="-"), count[selectCl])
    df.bs$cluster <- as.factor(df.bs$cluster)
    ##########one-stage bootstrap 
    est1 <- est.Wcdf(df.bs$x, df.bs$cluster, ri=ri, std = F)
    if(subboot){
      #sample level 1
      l1idx <- tapply(seq(nrow(df.bs)), df.bs[,'cluster'], I)
      idxselect1 <- lapply(l1idx, function(x) sample(x, size = length(x), replace = T))
      df.bs <- df.bs[unlist(idxselect1),]
      ##########two-stage bootstrap 
      est2 <- est.Wcdf(x=df.bs$x, cluster = df.bs$cluster, ri=ri, std = F)
      est[i,] <- c(est1,est2)
    }
    else{
      est[i,] <- c(est1)
    }
  }
  return(est)
}

#CI after fisher transformation 
fisherCIdelta <- function(x, s){
  tr <- log((1+x)/(1-x))/2 
  ts <- s/((1+x)*(1-x))
  l <- tr - 1.96 * ts
  u <- tr + 1.96 * ts
  l <- (exp(2*l)-1)/(exp(2*l)+1) #tanh()
  u <- (exp(2*u)-1)/(exp(2*u)+1)
  cbind(l,u)
}

#Coverage of CI
CI_coverage <- function(r, x, s, fisher=T, tol = 1e-7){
  l <- x - 1.96 * s 
  u <- x + 1.96 * s
  p <- mean((l - tol  <= r) & (r <= u + tol))
  if(fisher){
    ci <- fisherCIdelta(x, s)
    p <- c(mean((l - tol  <= r) & (r <= u + tol)),
           mean((ci[,1] - tol <= r) & (r <= ci[,2] + tol)))
  }
  return(p)
}

# Scenario I fixed n 
simulation <- function(iter, n.cll, h){
  ri <- seq(0, 1, 0.1)
  us <- rep(1, length(ri))
  rs <- (1-ri)*us/ri
  us[ri == 0] <- 0
  rs[ri == 0] <- 10
  sim_num <- 5
  n <- 100

  for(j in seq_along(us)){
    print(paste("ncll",n.cll[1],"-", n.cll[2], "iter", iter, "idx", j))

    u <- us[j]; r <- rs[j]
    ans <- list()
    trial <- (iter - 1) * sim_num
    for(i in 1:sim_num){
      dat <- generate.data.s1(seed=trial+i, n=n, cluster.limit=n.cll, mu=1, vu=u, vr=r, hlf=h)
      if(n.cll[1] == n.cll[2]){
        ec <- est.Wcdf(dat$x, dat$cluster, ri=1, std=T)
        ans[[i]] <- c(ec$est, ec$se)
      }
      else{
        ei <- est.Wcdf(dat$x, dat$cluster, ri=0, std=T)
        ec <- est.Wcdf(dat$x, dat$cluster, ri=1, std=T)
        other <- est.iter(dat$x, dat$cluster, ri=1, std=T)
        ans[[i]] <- list("est" = c(ei$est, ec$est, other$est),
                         "se" = c(ei$se, ec$se, other$se),
                         "Niter" = other$Niter)
      }
      save(ans, file = paste('output/ncll', n.cll[1], '-', n.cll[2], '-', 'idx', j-1, '-', iter,'.RData', sep=""))
    }
  }
}

for(i in 1:200){
  simulation(i, c(30,30), h=F)
  simulation(i, c(2,50), h=F)
  simulation(i, c(2,30), h=T)
}

# Scenario I fixed cluster sizes
simulation <- function(iter, n){

  ri <- seq(0, 1, 0.1)
  us <- rep(1, length(ri))
  rs <- (1-ri)*us/ri
  us[ri == 0] <- 0
  rs[ri == 0] <- 10

  sim_num <- 5
  n.cll <- c(30,30); h <- F

  for(j in seq_along(us)){
    u <- us[j]; r <- rs[j]; idx <- ri[j]*10
    ans <- list()
    trial <- (iter - 1) * sim_num
    for(i in 1:sim_num){
      dat <- generate.data.s1(seed=trial+i, n=n, cluster.limit=n.cll, mu=1, vu=u, vr=r, hlf=h)
      ec <- est.Wcdf(dat$x, dat$cluster, ri=1, std=T)
      ans[[i]] <- c(ec$est, ec$se)
      save(ans, file = paste('output/n',n, '-', 'idx', idx, '-', iter,'.RData', sep=""))
    }
  }
}

for(j in c(25, 50, 100, 200, 500, 1000)){
  for(i in 1:200){
    simulation(i, j)
  }
}


# Scenario III fixed cluster sizes (check generate.data()) 
simulation <- function(iter, n){

  mu <- 1
  ri <- seq(0,1,0.1); us <- rep(1, length(ri))
  rs <- us*(1 - ri)/ri
  us[ri==0] <- 0
  rs[ri==0] <- 10
  us <- log(1/2+sqrt(us*exp(-2*mu) + 1/4))

  sim_num <- 5
  n.cll <- c(30,30); h <- F

  for(j in seq_along(us)){
    u <- us[j]; r <- rs[j]; idx <- ri[j]*10
    ans <- list()
    trial <- (iter - 1) * sim_num
    for(i in 1:sim_num){
      dat <- generate.data.s3(seed=trial+i, n=n, cluster.limit=n.cll, mu=mu, vu=u, vr=r, hlf=h)
      ec <- est.Wcdf(dat$x, dat$cluster, ri=1, std=T)
      ans[[i]] <- c(ec$est, ec$se)
      save(ans, file = paste('output/n',n, '-', 'idx', idx, '-', iter,'.RData', sep=""))
    }
  }
}

for(j in c(25, 50, 100, 200, 500, 1000)){
  for(i in 1:200){
    simulation(i, j)
  }
}

# 
# # Scenario III fixed n (check generate.data()) 
simulation <- function(iter, n.cll, h){

  mu <- 1
  ri <- seq(0,1,0.1); us <- rep(1, length(ri))
  rs <- us*(1 - ri)/ri
  us[ri==0] <- 0
  rs[ri==0] <- 10
  us <- log(1/2+sqrt(us*exp(-2*mu) + 1/4))

  sim_num <- 5
  n <- 200

  for(j in seq_along(us)){
    print(paste("ncll",n.cll[1],"-", n.cll[2], "iter", iter, "idx", j))

    u <- us[j]; r <- rs[j]
    ans <- list()
    trial <- (iter - 1) * sim_num
    for(i in 1:sim_num){
      dat <- generate.data.s3(seed=trial+i, n=n, cluster.limit=n.cll, mu=mu, vu=u, vr=r, hlf=h)
      if(n.cll[1] == n.cll[2]){
        ec <- est.Wcdf(dat$x, dat$cluster, ri=1, std=T)
        ans[[i]] <- c(ec$est, ec$se)
      }
      else{
        ei <- est.Wcdf(dat$x, dat$cluster, ri=0, std=T)
        ec <- est.Wcdf(dat$x, dat$cluster, ri=1, std=T)
        other <- est.iter(dat$x, dat$cluster, ri=1, std=T)
        ans[[i]] <- list("est" = c(ei$est, ec$est, other$est),
                         "se" = c(ei$se, ec$se, other$se),
                         "Niter" = other$Niter)
      }
      save(ans, file = paste('output/ncll', n.cll[1], '-', n.cll[2], '-', 'idx', j-1, '-', iter,'.RData', sep=""))
    }
  }
}

for(i in 1:200){
  simulation(i, c(30,30), h=F)
  simulation(i, c(2,50), h=F)
  simulation(i, c(2,30), h=T)
}


# 
# #############Scenario 4, unequal within-variance, fixed number of clusters 
simulation <- function(iter, n.cll, k="same"){
  ri <- seq(0, 1, 0.1)
  us <- rep(1, length(ri))
  if(k=="diff"){
    k1 <- n.cll[1]/sum(n.cll); a1 <- 0.5; a2 <- 1.5
    rs <- (us/ri-us)/(k1*a1+a2-k1*a2)
  }
  if(k=="same"){
    k1 <- 1/2; a1 <- 0.5; a2 <- 1.5
    rs <- (us/ri-us)/(k1*a1+a2-k1*a2)
  }
  us[ri == 0] <- 0
  rs[ri == 0] <- 10

  sim_num <- 5
  n <- 200

  for(j in seq_along(us)){
    print(paste("ncll",n.cll[1],"-", n.cll[2], "iter", iter, "idx", j))

    u <- us[j]; r <- rs[j]
    ans <- list()
    trial <- (iter - 1) * sim_num
    for(i in 1:sim_num){
      dat <- generate.data.s4(seed=trial+i, n=n, cluster.limit=n.cll, mu=1, vu=u, vr1=a1*r,vr2=a2*r)

      ei <- est.Wcdf(dat$x, dat$cluster, ri=0, std=T)
      ec <- est.Wcdf(dat$x, dat$cluster, ri=1, std=T)
      other <- est.iter(dat$x, dat$cluster, ri=1, std=T)
      ans[[i]] <- list("est" = c(ei$est, ec$est, other$est),
                       "se" = c(ei$se, ec$se, other$se),
                       "Niter" = other$Niter)

      save(ans, file = paste('output/',k,'-n.cll', n.cll[1], '-', n.cll[2], '-', 'idx', j-1, '-', iter,'.RData', sep=""))
    }
  }
}

for(i in 1:200){
  simulation(i, c(2,30), k="same")
  simulation(i, c(2,30), k="diff")
}


##Under scenario 4, populations favor iterative methods 
simulation <- function(iter, n.cll, k1, lb){
  ri <- seq(0, 1, 0.1)
  us <- rep(1, length(ri))
  a1 <- 0.5; a2 <- 1.5
  rs <- (us/ri-us)/(k1*a1+a2-k1*a2)
  us[ri == 0] <- 0
  rs[ri == 0] <- 10
  sim_num <- 5
  n <- 200
  for(j in seq_along(us)){
    print(paste("ncll",n.cll[1],"-", n.cll[2], "iter", iter, "idx", j))

    u <- us[j]; r <- rs[j]
    ans <- list()
    trial <- (iter - 1) * sim_num
    for(i in 1:sim_num){
      dat <- generate.data.s4(seed=trial+i, n=n, cluster.limit=n.cll, mu=1, vu=u, vr1=a1*r,vr2=a2*r)

      ei <- est.Wcdf(dat$x, dat$cluster, ri=0, std=T)
      ec <- est.Wcdf(dat$x, dat$cluster, ri=1, std=T)
      other <- est.iter(dat$x, dat$cluster, ri=ec$est, std=T)
      ans[[i]] <- list("est" = c(ei$est, ec$est, other$est),
                       "se" = c(ei$se, ec$se, other$se),
                       "Niter" = other$Niter)

      save(ans, file = paste('output/', lb,'-n.cll', n.cll[1], '-', n.cll[2], '-', 'idx', j-1, '-', iter,'.RData', sep=""))
    }
  }
}

for(i in 1:200){
  simulation(i, c(2,30), 9/32, "betw")
}


simulation <- function(iter, n.cll, lb){
  ri <- seq(0, 1, 0.1)
  us <- rep(1, length(ri))
  sim_num <- 5
  n <- 200
  if(lb == "combo") k1 <- (1 + 7 * ri) / 16
  else if(lb == "ess"){
    w1 <- 2/(n/(1+ri)+15*n/(1+29*ri))/(1+ri)
    w2 <- 30/(n/(1+ri)+15*n/(1+29*ri))/(1+29*ri)
    k1 <- w1/(w1+w2)
  }
  a1 <- 0.5; a2 <- 1.5
  rs <- (us/ri-us)/(k1*a1+a2-k1*a2)
  us[ri == 0] <- 0
  rs[ri == 0] <- 10
  for(j in seq_along(us)){
    print(paste("ncll",n.cll[1],"-", n.cll[2], "iter", iter, "idx", j))

    u <- us[j]; r <- rs[j]
    ans <- list()
    trial <- (iter - 1) * sim_num
    for(i in 1:sim_num){
      dat <- generate.data.s4(seed=trial+i, n=n, cluster.limit=n.cll, mu=1, vu=u, vr1=a1*r,vr2=a2*r)

      ei <- est.Wcdf(dat$x, dat$cluster, ri=0, std=T)
      ec <- est.Wcdf(dat$x, dat$cluster, ri=1, std=T)
      other <- est.iter(dat$x, dat$cluster, ri=ec$est, std=T)
      ans[[i]] <- list("est" = c(ei$est, ec$est, other$est),
                       "se" = c(ei$se, ec$se, other$se),
                       "Niter" = other$Niter)

      save(ans, file = paste('output/', lb,'-n.cll', n.cll[1], '-', n.cll[2], '-', 'idx', j-1, '-', iter,'.RData', sep=""))
    }
  }
}

for(i in 1:200){
  simulation(i, c(2,30), 9/32, "ess")
}


#########################Simulations Scenario 4 and finite populations
ri <- seq(0, 1, 0.1)
us <- rep(1, length(ri))
k1 <- 1/2; k2 <- 1 - w1
a1 <- 1; a2 <- 1
ri1 <- ri - 0.05
rs <- (us*(1-ri))/(k1*a1+k2*a2)/(ri-ri1)

# ####################unequal cluster sizes, used for comparison with bootstrap 
simulation <- function(iter, n){
  mu <- 1
  ri <- 0.5
  u <- 1
  r <- u*(1 - ri)/ri
  sim_num <- 5
  n.cll <- c(2,50); h <- F
  trial <- (iter - 1) * sim_num
  load(paste('output/ncll-2-50-n',n, '-', iter,'-val05.RData', sep=""))
  for(i in 1:sim_num){
      print(trial+i)
      dat <- generate.data.s1(seed=trial+i, n=n, cluster.limit=n.cll, mu=mu, vu=u, vr=r, hlf=h)
      ei <- est.Wcdf(dat$x, dat$cluster, ri=0, std=T)
      ans[[i]] <- c(ans[[i]], ei$est, ei$se)
      save(ans, file = paste('output/ncll-2-50-n',n, '-', iter,'-val05.RData', sep=""))
  }
}

for(j in c(25, 50, 100, 200, 500, 1000)){
  for(i in 1:200){
    simulation(i, j)
  }
}

####################Bootstrap performance 
####True values at the boundary
simulation <- function(iter, n, n.cll, bootN=200, h=F){
  sim_num <- 5
  ri <- c(0.1, 0.9)
  us <- c(1,1)
  rs <- us*(1 - ri)/ri
  for(j in seq_along(us)){
    u <- us[j];r <- rs[j]
    ans <- list()
    trial <- (iter - 1) * sim_num
    print(paste("j is ", j, sep=""))
    for(i in 1:sim_num){
      print(trial+i)
      dat <- generate.data.s1(seed=trial+i, n=n, cluster.limit=n.cll, mu=1, vu=u, vr=r, hlf=h)
      est <- est.Wcdf.boot(seed=trial+i, dat, N=bootN, ri=1, subboot=T)
      ans[[i]] <- est
      save(ans, file = paste('output/n',n, '-', 'idx', j, '-', iter,'.RData', sep=""))
    }
  }
}

for(i in 1:200){
  simulation(i, n=100, n.cll=c(30,30))
}


####Unequal sizes
simulation <- function(iter, n, n.cll, bootN=200, h=F){
  sim_num <- 5
  ri <- 0.5
  u <- 1
  r <- u*(1 - ri)/ri
  ans <- list()
  trial <- (iter - 1) * sim_num
  for(i in 1:sim_num){
      print(trial+i)
      dat <- generate.data.s1(seed=trial+i, n=n, cluster.limit=n.cll, mu=1, vu=u, vr=r, hlf=h)
      ec <- est.Wcdf.boot(seed=trial+i, dat, N=bootN, ri=1, subboot=T)
      ei <- est.Wcdf.boot(seed=trial+i, dat, N=bootN, ri=0, subboot=T)
      ans[[i]] <- cbind(ec, ei)
      save(ans, file = paste('output/Boot-ncll-2-50-n',n,'-', iter,'.RData', sep=""))
  }
}


for(i in 1:200){
  simulation(i, n=100, n.cll=c(2,50))
}

#####################Finite population 
generate.data <- function(seed=1234, n=30, cluster.limit=c(10,10), mu=1, vu=1,vr=1, hlf=F){
  set.seed(seed)
  #cluster size
  if(hlf) size.cluster <- rep(cluster.limit, each = n/2)
  else{
    if(cluster.limit[1] == cluster.limit[2]){
      size.cluster <- rep(cluster.limit[1], n)
    }else size.cluster <- replicate(n, sample(cluster.limit[1]:cluster.limit[2], 1))
  }
  #generate data for each cluster
  x <- list()
  cluster <- list()
  id <- list()
  #generate cluster means
  u <- rnorm(n, mean = mu, sd = sqrt(vu))
  for(i in 1:n){
    n.cluster <- size.cluster[i]
    rw <-  1 / (1 - n.cluster)
    s <- matrix(rw * vr, nrow = n.cluster - 1, ncol = n.cluster - 1)
    diag(s) <- vr
    r <- rep(NA, n.cluster)
    #r[1:(n.cluster-1)] <- mvtnorm::rmvnorm(1, mean = rep(0, n.cluster - 1), sigma = s)
    r[1:(n.cluster-1)] <- MASS::mvrnorm(1, mu = rep(0, n.cluster - 1), Sigma = s)
    r[n.cluster] <- - sum(r[1:(n.cluster-1)])
    x[[i]] <- r+u[i]
    cluster[[i]] <- rep(i, n.cluster)
    id[[i]] <- 1:n.cluster
  }
  dat <- data.frame("x"= unlist(x),
                    "cluster"=as.factor(unlist(cluster)),
                    "id"=as.factor(unlist(id)))
  return(dat)
}

simulation <- function(iter, n.cll, h=F){
  rw <- mean(- 1 / (n.cll - 1))
  ri <- seq(-1, 1, 0.1)
  us <- rep(1, length(ri))
  rs <- (1-ri)*us/(ri-rw)
  us[ri == -1] <- 0
  rs[ri == -1] <- 20
  sim_num <- 5
  
  n <- 200
  for(j in seq_along(us)){
    print(paste("ncll",n.cll[1],"-", n.cll[2], "iter", iter, "idx", j))
    u <- us[j]; r <- rs[j]
    ans <- list()
    trial <- (iter - 1) * sim_num
    for(i in 1:sim_num){
      dat <- generate.data(seed=trial+i, n=n, cluster.limit=n.cll, mu=1, vu=u, vr=r, hlf=h)
      eicc <- est.icc(rank(dat$x), dat$cluster)
      possibleError <- tryCatch(
        estlmer <- est.lmer(rank(dat$x), dat$cluster),
        error=function(e) e
      )
      if(inherits(possibleError, "error")) elmer <- NA
      else elmer <- estlmer$icc
      if(n.cll[1] == n.cll[2]){
        ec <- est.Wcdf(dat$x, dat$cluster, ri=1, std=T)
        ans[[i]] <- c(eicc$icc, elmer, ec$est, ec$se)
      }
      else{
        ei <- est.Wcdf(dat$x, dat$cluster, ri=0, std=T)
        ec <- est.Wcdf(dat$x, dat$cluster, ri=1, std=T)
        other <- est.iter(dat$x, dat$cluster, ri=1, std=T)
        ans[[i]] <- list("est" = c(eicc$icc,  elmer, ei$est, ec$est, other$est),
                         "se" = c(ei$se, ec$se, other$se),
                         "Niter" = other$Niter)
      }
      save(ans, file = paste('output/finite-allRho-ncll', n.cll[1], '-', n.cll[2], '-', 'idx', j-1, '-', iter,'.RData', sep=""))
    }
  }
}


for(i in 1:200){
    simulation(i, c(2, 2))
}


#####################Scenario 5 discrete 
simulation <- function(iter, n, ls){
  
  # ri <- seq(0, 1, 0.1)
  ri <- c(0.925, 0.95)
  us <- rep(1, length(ri))
  rs <- (1-ri)*us/ri
  us[ri == 0] <- 0
  rs[ri == 0] <- 10
  sim_num <- 5
  n.cll <- c(30,30); h <- F
  
  for(j in seq_along(us)){
    u <- us[j]; r <- rs[j]; idx <- ri[j]*10
    ans <- list()
    trial <- (iter - 1) * sim_num
    for(i in 1:sim_num){
      dat <- generate.data.s1(seed=trial+i, n=n, cluster.limit=n.cll, mu=1, vu=u, vr=r, hlf=h)
      est <- list()
      for(li in seq_along(ls)){
        l <- ls[li]
        bi <- c(-Inf, quantile(dat$x, seq(1/l, 1-1/l, 1/l)), Inf)
        xi <- as.numeric(cut(dat$x, breaks = bi, labels = 1:l))  
        ec <- est.Wcdf(xi, dat$cluster, ri=1, std=T)
        est[[li]] <- c(ec$est, ec$se)
      }
      est <- unlist(est)
      names(est) <- rep(ls, each = 2)
      ans[[i]] <- est
      save(ans, file = paste('output/Discrete-n',n, '-', 'idx', idx, '-', iter,'.RData', sep=""))
    }
  }
}

for(j in c(25, 50, 100, 200, 500, 1000)){
  for(i in 1:200){
    simulation(i, n=j, c(30, 30))
  }
}

