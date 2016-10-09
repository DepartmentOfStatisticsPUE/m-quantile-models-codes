# Model-Based Simulations
#rm(list=ls(all=TRUE))
library(MASS)
library(SemiPar)
library(nlme)

my.default.knots.2D <- function (x1, x2, num.knots)
{
  require("cluster")
  if (missing(num.knots))
    num.knots <- max(10, min(50, round(length(x1) / 4)))
  X <- cbind(x1, x2)
  dup.inds <- (1:nrow(X))[dup.matrix(X) == T]
  if (length(dup.inds) > 0)
    X <- X[-dup.inds,]
  knots <- clara(X, num.knots)$medoids
  return(knots)
}
#
#
#
Z.matrix <- function(lon, lat, knots) {
  K <- nrow(knots)
  dist.knot <- matrix(0, K, K)
  dist.knot[lower.tri(dist.knot)] <- dist(knots)
  dist.knot <- dist.knot + t(dist.knot)
  Omega <- tps.cov(dist.knot)
  dist.lon <- outer(lon, knots[, 1], "-")
  dist.lat <- outer(lat, knots[, 2], "-")
  dist.x <- sqrt(dist.lon ^ 2 + dist.lat ^ 2)
  svd.Omega <- svd(Omega)
  sqrt.Omega <- t(svd.Omega$v %*% (t(svd.Omega$u) * sqrt(svd.Omega$d)))
  Z <- t(solve(sqrt.Omega, t(tps.cov(dist.x))))
  return(Z)
}

npqrlm <-
  function(Y.n,
           X,
           Z.spline,
           quantile,
           tol = 0.001,
           maxit = 100,
           theta = 2,
           kk = 1.345)
  {
    # prototype function for panel data fitting of MQ models
    # the matrix X is assumed to contain an intercept
    # the vector s is a strata indicator assumed (so far) to be a one-way layout
    #theta : GCV parameter, defaulted to 2
    
    require(SemiPar)
    require(MASS)
    require(splines)
    assign("tol", tol, pos = 1)
    assign("maxit", maxit, pos = 1)
    assign("Y.n", Y.n, pos = 1)
    assign("X", X, pos = 1)
    assign("kk", kk, pos = 1)
    assign("Z.spline", Z.spline, pos = 1)
    assign("theta", theta, pos = 1)
    assign("quantile", quantile, pos = 1)
    
    n <- length(Y.n)
    X <- as.matrix(X)
    p1 <- ncol(X)
    p2 <- ncol(Z.spline)
    X.n <- cbind(X, Z.spline)
    p = p1 + p2
    
    b = rep(1, p)
    
    my.psi.q <- function(u, q, c) {
      s <- median(abs(u)) / 0.6745
      w <- psi.huber((u / s), c)
      ww <- 2 * (1 - q) * w
      ww[u > 0] <- 2 * q * w[u > 0]
      w <- ww
      w * u
    }
    
    my.b <- function(X, Y, W, lambda)
    {
      G <- as.matrix(diag(c(rep(0, p1), rep(1, p2)), p))
      solve(t(X) %*% W %*% X + G * lambda) %*% t(X) %*% W %*% Y
    }
    
    
    stima <- function(l) {
      # l : coefficiente di penalizzazione nell'IRPLS
      n <- nrow(X.n)
      diff <- 1
      iter <- 0
      while (diff > tol)
      {
        #inizia procedura di stima
        iter <- iter + 1
        res <- Y.n - X.n %*% b #calcolo residui
        W.n <-
          as.matrix(diag(c(
            my.psi.q(as.matrix(res), qtl, kk) / as.matrix(res)
          ), n))
        assign("W.n", W.n, pos = 1)
        b.ott <- my.b(X.n, Y.n, W.n, l)
        diff <- sum((as.matrix(b) - as.matrix(b.ott)) ^ 2)
        b <- b.ott
        if (iter > maxit)
        {
          warning(paste("failed to converge in", maxit, "steps at q = ", qtl))
          break
        }
      }
      y.hat = X.n %*% b
      list(
        fitted.values = as.matrix(y.hat),
        coef = as.matrix(b),
        we = as.matrix(diag(W.n))
      )
    }
    
    my.GCV <- function(l)
    {
      G <- as.matrix(diag(c(rep(0, p1), rep(1, p2)), p))
      tmp <- stima(l)
      y.hat <- tmp$fitted.values
      S <- (X.n) %*% solve(t(X.n) %*% W.n %*% X.n + G * l) %*% t(X.n) %*% W.n
      sum((Y.n - y.hat) ^ 2) / ((1 - theta * sum(diag(S)) / n) ^ 2)
    }
    
    length.q <- length(quantile)
    y.fit <- matrix(0, n, length.q)
    y.coef <- matrix(0, p, length.q)
    y.weight <- matrix(0, n, length.q)
    lambda.ott <- NULL
    for (k in 1:length.q)
    {
      qtl <- quantile[k]
      qtl <- assign("qtl", qtl, pos = 1)
      tmp <- optimize(my.GCV, c(0, 50))
      l.ott <- tmp$minimum
      lambda.ott[k] <- l.ott
      y.stim <- stima(l.ott)
      y.fit[, k] <- y.stim$fitted.values
      y.coef[, k] = y.stim$coef
      y.weight[, k] = y.stim$we
    }
    list(
      hat.values = y.fit,
      b.stim = y.coef,
      q.weights = y.weight,
      lambda.q = lambda.ott
    )
  }

# COMPUTE THE QUANTILE ORDER

# COMPUTING OF THE QUANTILE-ORDERS
"zerovalinter" <- function(y, x)
{
  if (min(y) > 0) {
    xmin <- x[y == min(y)]
    if (length(xmin) > 0)
      xmin <- xmin[length(xmin)]
    xzero <- xmin
  }
  
  
  
  else {
    if (max(y) < 0) {
      xmin <- x[y == max(y)]
      if (length(xmin) > 0)
        xmin <- xmin[1]
      xzero <- xmin
    }
    else {
      y1 <- min(y[y > 0])
      if (length(y1) > 0)
        y1 <- y1[length(y1)]
      y2 <- max(y[y < 0])
      if (length(y2) > 0)
        y2 <- y2[1]
      x1 <- x[y == y1]
      if (length(x1) > 0)
        x1 <- x1[length(x1)]
      x2 <- x[y == y2]
      if (length(x2) > 0)
        x2 <- x2[1]
      xzero <- (x2 * y1 - x1 * y2) / (y1 - y2)
      xmin <- x1
      if (abs(y2) < y1)
        xmin <- x2
    }
  }
  resu <-  xzero
  resu
}

# Function for Finding the Quantile Orders by Linear Interpolation
# Assumes that "zerovalinter" function has been already loaded

"gridfitinter" <- function(y, expectile, Q)
  # computing of the expectile-order of each observation of y by interpolation
{
  nq <- length(Q)
  diff <- y %*% t(as.matrix(rep(1, nq))) - expectile
  vectordest <- apply(diff, 1, zerovalinter, Q)
  
  #print(vectordest)
  #qord<-list(ord=c(vectordest))
  #qord
}


#y: study variable
#x: set of covariates without the intercept for sampled units
#regioncode.s: area code for sampled units
#x.r: set of covariates for out of sample units
#regioncode.r: area code for out of sample units
#p size of x +1 (intercept)
#z.spline: the spline matrix
sae.npmq = function(y,
                    x,
                    z.spline,
                    z.spline.r,
                    regioncode.s,
                    m,
                    p,
                    x.outs,
                    regioncode.r,
                    tol.value = 0.0001,
                    maxit.value = 100,
                    k.value = 1.345)
{
  MQE <- c(rep(0, m))
  MQNAIVE <- c(rep(0, m))
  
  
  datanew = cbind(y, x, regioncode.s)
  
  ni = as.numeric(table(regioncode.s))
  
  sample.sizer <- as.numeric(table(regioncode.r))
  
  Ni = sample.sizer + ni
  
  
  N <- sum(Ni)
  n <- sum(ni)
  
  x = matrix(x, n, p - 1)
  x.r = matrix(x.outs, (N - n), p - 1)
  x.t = rbind(x, x.r)
  x.c <- rep(1, n)
  x.design <- cbind(x.c, x)
  p = ncol(x.design)
  p2 <- ncol(z.spline)
  G.matrix = diag(c(rep(0, p), rep(1, p2)))
  
  ob <-
    npqrlm(
      y,
      x.design,
      z.spline,
      quantile = sort(c(
        seq(0.006, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98
      )),
      kk = k.value,
      maxit = maxit.value,
      tol = tol.value
    )
  
  q.values <-
    sort(c(seq(0.006, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98))
  
  qo <- matrix(c(gridfitinter(y, ob$hat.values, q.values)), nrow = n, ncol =
                 1)
  
  qmat <- matrix(c(qo, regioncode.s), nrow = sum(ni), ncol = 2)
  
  mqo <- aggregate(qmat[, 1], list(d2 = qmat[, 2]), mean)[, 2]
  
  uar <- sort(unique(regioncode.s))
  saq <- matrix(c(mqo, uar), nrow = m, ncol = 2)
  
  saq <- rbind(saq, c(0.5, 9999))
  
  beta.stored = matrix(0, m, 2)
  res.s = NULL
  tttmp1 <- NULL
  ci <- array(rep(0, n * m), dim = c(n, m, 1))
  ci1 <- array(rep(0, n * m), dim = c(n, m, 1))
  prs <- NULL
  prr <- NULL
  wi <- matrix(0, n, m)
  
  for (i in 1:m) {
    ob1 <-
      npqrlm(
        y,
        x.design,
        z.spline,
        quantile = mqo[i],
        kk = k.value,
        maxit = maxit.value,
        tol = tol.value
      )
    
    
    wd <- diag(c(ob1$q.weights))
    
    # Regional parameters from multiquantile model
    
    coef <-
      matrix(c(t(ob1$b.stim)), nrow = 1, ncol = p + p2) # need to be ordered by area
    
    coef <- t(coef)
    
    meat <-
      wd %*% cbind(x.design, z.spline) %*% solve(t(cbind(x.design, z.spline)) %*%
                                                   wd %*% cbind(x.design, z.spline) + ob1$lambda.q * G.matrix)
    
    x1 <- c(rep(1, (Ni[i] - ni[i])))
    
    ir <- rep(0, n)
    
    ir[regioncode.s == uar[i]] <- 1
    
    rj1 <- sample.sizer[i]
    
    r = NULL
    
    for (kk in 1:(p - 1))
    {
      r <- c(r, sum(x.r[, kk][regioncode.r == uar[i]]))
    }
    
    r = c(rj1, r, c(as.numeric(apply(
      z.spline.r[regioncode.r == uar[i], ], 2, sum
    ))))
    
    
    sj1 <- sum(rep(1, ni[i]))
    
    tss = NULL
    
    for (kk in 1:(p - 1))
    {
      tss <- c(tss, sum(x[, kk][regioncode.s == uar[i]]))
    }
    
    
    tss <-
      c(sj1, tss, c(as.numeric(apply(
        z.spline[regioncode.s == uar[i], ], 2, sum
      ))))
    
    w.welsh <- ((Ni[i]) / (ni[i])) * ir + meat %*% (r - ((Ni[i] - ni[i]) / ni[i]) *
                                                      tss)
    
    
    MQE[i] <- sum(w.welsh * y) / sum(w.welsh)
    y.i <- y[regioncode.s == uar[i]]
    
    y.pred.s <-
      cbind(1, x[regioncode.s == uar[i], ], z.spline[regioncode.s == uar[i], ]) %*%
      coef
    residual <- y.i - y.pred.s
    
    tttmp1[i] <- (sample.sizer[i] / ni[i]) * sum(residual)
    
    prs <- c(prs, y.pred.s)
    res.s <- c(res.s, residual)
    y.pred <-
      cbind(1, x.r[regioncode.r == uar[i], ], z.spline.r[regioncode.r == uar[i], ]) %*%
      coef
    prr <- c(prr, y.pred)
    
    data <- cbind(regioncode.s, w.welsh)
    
    for (kk in 1:n) {
      if (data[kk, 1] == uar[i])
        ci[kk, i, 1] <- data[kk, 2] - 1
      else if (data[kk, 1] != uar[i])
        ci[kk, i, 1] <- data[kk, 2]
    }
    
    MQNAIVE[i] <- (1 / Ni[i]) * as.real(sum(y.i) + sum(y.pred))
    
    #f<-(sample.sizer[i])/ni[i]
    ai <- r
    bi <-
      ai %*% solve(t(cbind(x.design, z.spline)) %*% wd %*% cbind(x.design, z.spline)) %*%
      t(cbind(x.design, z.spline)) %*% wd
    bi <- c(bi)
    wi[, i] <- c(ir + bi)
    
    datanaive <- cbind(regioncode.s, wi[, i])
    for (kk in 1:n) {
      if (datanaive[kk, 1] == uar[i])
        ci1[kk, i, 1] <- datanaive[kk, 2] - 1
      else if (datanaive[kk, 1] != uar[i])
        ci1[kk, i, 1] <- datanaive[kk, 2]
    }
    
    
    
  }
  
  res.d = res.s ^ 2
  
  res.d <- cbind(res.d, regioncode.s, ci[, , 1])
  
  v <- NULL
  for (oo in 1:m) {
    v[oo] <-
      1 / Ni[oo] ^ 2 * (sum((res.d[, (oo + 2)][res.d[, 2] == uar[oo]] ^ 2 + (sample.sizer[oo]) /
                               n) * res.d[, 1][res.d[, 2] == uar[oo]]) + sum(res.d[, (oo + 2)][res.d[, 2] !=
                                                                                                 uar[oo]] ^ 2 * res.d[, 1][res.d[, 2] != uar[oo]]))
    
  }
  
  res.d1 <- cbind(res.s ^ 2, regioncode.s, ci1[, , 1])
  v1 <- NULL
  bias <- NULL
  mse <- NULL
  for (oo in 1:m) {
    v1[oo] <-
      1 / Ni[oo] ^ 2 * (sum((res.d1[, (oo + 2)][res.d1[, 2] == uar[oo]] ^ 2 +
                               (sample.sizer[oo]) / n) * res.d1[, 1][res.d1[, 2] == uar[oo]]) + sum(res.d1[, (oo +
                                                                                                                2)][res.d1[, 2] != uar[oo]] ^ 2 * res.d1[, 1][res.d1[, 2] != uar[oo]]))
    
    bias[oo] <-
      (1 / Ni[oo]) * (sum(wi[, oo] * prs) - sum(c(prs[regioncode.s == uar[oo]], prr[regioncode.r ==
                                                                                      uar[oo]])))
    mse[oo] <- v1[oo] + (bias[oo]) ^ 2
  }
  
  
  list(
    npmq.cd = MQE,
    npmq.naive = MQNAIVE,
    mse.cd = v,
    mse.naive = mse,
    code.area = uar
  )
}
