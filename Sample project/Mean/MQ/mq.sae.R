# Model-Based Simulations
rm(list = ls(all = TRUE))
library(MASS)

QRLM <-
  function (x,
            y,
            case.weights = rep(1, nrow(x)),
            var.weights = rep(1, nrow(x)),
            ...,
            w = rep(1, nrow(x)),
            init = "ls",
            psi = psi.huber,
            scale.est = c("MAD", "Huber", "proposal 2"),
            k2 = 1.345,
            method = c("M", "MM"),
            maxit = 20,
            acc = 1e-04,
            test.vec = "resid",
            q = 0.5)
  {
    irls.delta <-
      function(old, new)
        sqrt(sum((old - new) ^ 2) / max(1e-20, sum(old ^ 2)))
    irls.rrxwr <- function(x, w, r) {
      w <- sqrt(w)
      max(abs((matrix(
        r * w, 1, length(r)
      ) %*% x) / sqrt(matrix(w, 1, length(
        r
      )) %*% (x ^ 2)))) / sqrt(sum(w * r ^ 2))
    }
    method <- match.arg(method)
    nmx <- deparse(substitute(x))
    if (is.null(dim(x))) {
      x <- as.matrix(x)
      colnames(x) <- nmx
    }
    else
      x <- as.matrix(x)
    if (is.null(colnames(x)))
      colnames(x) <- paste("X", seq(ncol(x)), sep = "")
    if (qr(x)$rank < ncol(x))
      stop("x is singular: singular fits are not implemented in rlm")
    if (!(any(test.vec == c("resid", "coef", "w", "NULL")) ||
          is.null(test.vec)))
      stop("invalid testvec")
    if (length(var.weights) != nrow(x))
      stop("Length of var.weights must equal number of observations")
    if (any(var.weights < 0))
      stop("Negative var.weights value")
    if (length(case.weights) != nrow(x))
      stop("Length of case.weights must equal number of observations")
    w <- (w * case.weights) / var.weights
    if (method == "M") {
      scale.est <- match.arg(scale.est)
      if (!is.function(psi))
        psi <- get(psi, mode = "function")
      arguments <- list(...)
      if (length(arguments)) {
        pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0)
        if (any(pm == 0))
          warning(paste("some of ... do not match"))
        pm <- names(arguments)[pm > 0]
        formals(psi)[pm] <- unlist(arguments[pm])
      }
      if (is.character(init)) {
        if (init == "ls")
          temp <- lm.wfit(x, y, w, method = "qr")
        else if (init == "lts")
          temp <- lqs.default(x, y, intercept = FALSE, nsamp = 200)
        else
          stop("init method is unknown")
        coef <- temp$coef
        resid <- temp$resid
      }
      else {
        if (is.list(init))
          coef <- init$coef
        else
          coef <- init
        resid <- y - x %*% coef
      }
    }
    else if (method == "MM") {
      scale.est <- "MM"
      temp <-
        lqs.default(x,
                    y,
                    intercept = FALSE,
                    method = "S",
                    k0 = 1.548)
      coef <- temp$coef
      resid <- temp$resid
      psi <- psi.bisquare
      if (length(arguments <- list(...)))
        if (match("c", names(arguments), nomatch = FALSE)) {
          c0 <- arguments$c
          if (c0 > 1.548) {
            psi$c <- c0
          }
          else
            warning("c must be at least 1.548 and has been ignored")
        }
      scale <- temp$scale
    }
    else
      stop("method is unknown")
    done <- FALSE
    conv <- NULL
    n1 <- nrow(x) - ncol(x)
    if (scale.est != "MM")
      scale <- mad(resid / sqrt(var.weights), 0)
    theta <- 2 * pnorm(k2) - 1
    gamma <- theta + k2 ^ 2 * (1 - theta) - 2 * k2 * dnorm(k2)
    qest <- matrix(0, nrow = ncol(x), ncol = length(q))
    qwt <- matrix(0, nrow = nrow(x), ncol = length(q))
    qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
    qres <- matrix(0, nrow = nrow(x), ncol = length(q))
    for (i in 1:length(q)) {
      for (iiter in 1:maxit) {
        if (!is.null(test.vec))
          testpv <- get(test.vec)
        if (scale.est != "MM") {
          if (scale.est == "MAD")
            scale <- median(abs(resid / sqrt(var.weights))) / 0.6745
          else
            scale <-
              sqrt(sum(pmin(
                resid ^ 2 / var.weights, (k2 * scale) ^ 2
              )) / (n1 * gamma))
          if (scale == 0) {
            done <- TRUE
            break
          }
        }
        w <- psi(resid / (scale * sqrt(var.weights))) * case.weights
        ww <- 2 * (1 - q[i]) * w
        ww[resid > 0] <- 2 * q[i] * w[resid > 0]
        w <- ww
        temp <- lm.wfit(x, y, w, method = "qr")
        coef <- temp$coef
        resid <- temp$residuals
        if (!is.null(test.vec))
          convi <- irls.delta(testpv, get(test.vec))
        else
          convi <- irls.rrxwr(x, wmod, resid)
        conv <- c(conv, convi)
        done <- (convi <= acc)
        if (done)
          break
      }
      if (!done)
        warning(paste("rlm failed to converge in", maxit, "steps at q = ", q[i]))
      qest[, i] <- coef
      qwt[, i] <- w
      qfit[, i] <- temp$fitted.values
      qres[, i] <- resid
    }
    list(
      fitted.values = qfit,
      residuals = qres,
      q.values = q,
      q.weights = qwt,
      coefficients = qest
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
mq_function = function(y,
                       x,
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

  ob <-
    QRLM(
      x.design,
      y,
      q = sort(c(
        seq(0.006, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98
      )),
      k = k.value,
      maxit = maxit.value,
      acc = tol.value
    )

  qo <-
    matrix(c(gridfitinter(y,
                          ob$fitted.values,
                          ob$q.values)),
           nrow = n,
           ncol = 1)

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
      QRLM(
        x.design,
        y,
        q = mqo[i],
        psi = psi.huber,
        k = k.value,
        maxit = maxit.value,
        acc = tol.value
      )

    wd <- diag(c(ob1$q.weights))

    # Regional parameters from multiquantile model

    coef <-
      matrix(c(t(ob1$coefficients)),
             nrow = 1,
             ncol = p) # need to be ordered by area

    coef <- t(coef)

    meat <-
      wd %*% x.design %*% solve(t(x.design) %*% wd %*% x.design)

    x1 <- c(rep(1, (Ni[i] - ni[i])))

    ir <- rep(0, n)

    ir[regioncode.s == uar[i]] <- 1

    rj1 <- sample.sizer[i]

    r = NULL

    for (kk in 1:(p - 1))
    {
      r <- c(r, sum(x.r[, kk][regioncode.r == uar[i]]))
    }

    r = c(rj1, r)


    sj1 <- sum(rep(1, ni[i]))

    tss = NULL

    for (kk in 1:(p - 1))
    {
      tss <- c(tss, sum(x[, kk][regioncode.s == uar[i]]))
    }


    tss <- c(sj1, tss)

    w.welsh <-
      ((Ni[i]) / (ni[i])) * ir + meat %*% (r - ((Ni[i] - ni[i]) / ni[i]) *
                                             tss)


    MQE[i] <- sum(w.welsh * y) / sum(w.welsh)
    y.i <- y[regioncode.s == uar[i]]

    y.pred.s <- cbind(1, x[regioncode.s == uar[i],]) %*% coef
    residual <- y.i - y.pred.s

    tttmp1[i] <- (sample.sizer[i] / ni[i]) * sum(residual)

    prs <- c(prs, y.pred.s)
    res.s <- c(res.s, residual)
    y.pred <- cbind(1, x.r[regioncode.r == uar[i],]) %*% coef
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
      ai %*% solve(t(x.design) %*% wd %*% x.design) %*% t(x.design) %*% wd
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
                               ni[oo]) * res.d[, 1][res.d[, 2] == uar[oo]]) + sum(res.d[, (oo + 2)][res.d[, 2] !=
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
    mq.cd = MQE,
    mq.naive = MQNAIVE,
    mse.cd = v,
    mse.naive = mse,
    code.area = uar
  )
}
