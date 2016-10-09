library(MASS)
library(nlme)
library(sp)
library(spgwr)

#M-estimator with GWR weights

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
            q = 0.5,
            w1)
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
        w <- ww * diag(w1)
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
#
zerovalinter <- function(y, x)
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
#
# Function for Finding the Quantile Orders by Linear Interpolation
#
gridfitinter <- function(y, expectile, Q)
  # computing of the expectile-order of each observation of y by interpolation
{
  nq <- length(Q)
  diff <- y %*% t(as.matrix(rep(1, nq))) - expectile
  vectordest <- apply(diff, 1, zerovalinter, Q)
}



##M-estimator original

QRLM1 <-
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


mqgwr.sae = function(x,
                     y,
                     m,
                     area,
                     lon,
                     lat,
                     x.r,
                     area.r,
                     lon.r,
                     lat.r,
                     method = "mqgwr",
                     k.value = 1.345,
                     mqgwrweight = TRUE)
{
  id.area <- sort(unique(area))
  m <- length(id.area)
  
  id.area.r = sort(unique(area.r))
  m.r <- length(id.area.r)
  
  tmp.cont = rep(0, m.r)
  for (i in 1:m.r)
  {
    for (j in 1:m)
    {
      if (id.area.r[i] == id.area[j])
        tmp.cont[i] = 1
    }
    
  }
  
  tmp0 = which(tmp.cont == 0)
  id.area.out = id.area.r[tmp0]
  id.area.in = id.area
  
  ni <- rep(0, m)
  for (i in 1:m)
    ni[i] <- sum(area == id.area.in[i])
  n <- sum(ni)
  
  ri <- rep(0, m)
  for (i in 1:m)
    ri[i] <- sum(area.r == id.area.in[i])
  r <- sum(ri)
  
  Ni = ri + ni
  
  RI.MQGWR_CD.Mean = rep(0, m)
  if (m.r > m)
    RI.MQGWR_CD.Mean.out = rep(0, (m.r - m))
  if (m.r == m)
    RI.MQGWR_CD.Mean.out = NULL
  mse = rep(0, m)
  
  # Compute the Distance Matrix
  
  eu_dist = as.matrix(dist(cbind(as.vector(lon), as.vector(lat))))
  x.design = cbind(1, x)
  p = ncol(x.design)
  q.value = sort(c(seq(0.002, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98))
  
  if (method == "mqgwr")
  {
    ob.trad <- QRLM1(x.design,
                     y,
                     maxit = 100,
                     q = q.value,
                     k = k.value)
    qo.trad <-
      matrix(c(gridfitinter(y, ob.trad$fitted.values, q.value)), nrow = n, ncol =
               1)
    
    if (mqgwrweight == TRUE)
      band = gwr.sel(y ~ 1 + x, coords = cbind(lon, lat), gweight = gwr.gauss)
    
    if (mqgwrweight == FALSE)
      band = gwr.sel(y ~ 1 + x, coords = cbind(lon, lat), gweight = gwr.bisquare)
    
    if (mqgwrweight == TRUE)
      w.sp <- gwr.gauss((eu_dist ^ 2), band)
    if (mqgwrweight == FALSE)
      w.sp <- gwr.bisquare((eu_dist ^ 2), band)
    q.new = 0
    for (ii in 1:n) {
      w.new <- diag(w.sp[, ii])
      ob <-
        QRLM(
          x.design,
          y,
          maxit = 100,
          q = sort(c(
            seq(0.002, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98
          )),
          w1 = w.new,
          k = k.value
        )
      qo <-
        matrix(c(gridfitinter(y, ob$fitted.values, ob$q.values)), nrow = n, ncol =
                 1)
      q.new[ii] = as.real(qo[ii, 1])
      if (is.na(q.new[ii]))
        q.new[ii] = as.real(qo.trad[ii, 1])
    }
    
    
    qmat1 <- matrix(c(q.new, area), nrow = n, ncol = 2)
    
    mqo1 = tapply(qmat1[, 1], qmat1[, 2], mean)
    
    saq <- matrix(0, nrow = m, ncol = 2)
    
    saq[, 1] = mqo1
    
    saq[, 2] = sort(unique(qmat1[, 2]))
    
    ci = array(rep(0, n * m), dim = c(n, m, 1))
    
    res.s = NULL
    
    
    for (i in 1:m)
    {
      pred.medr = 0
      x.r.area <-
        matrix(cbind(1, x.r)[area.r == id.area.in[i]], ri[i], p)
      lon.gwr = lon.r[area.r == id.area.in[i]]
      lat.gwr = lat.r[area.r == id.area.in[i]]
      
      tmp = matrix(0, n, 1)
      tmp1 = matrix(0, n, 1)
      
      for (j in 1:ri[i]) {
        dbase = as.matrix(rbind(
          cbind(lon.gwr[j], lat.gwr[j]),
          cbind(as.vector(lon), as.vector(lat))
        ))
        dist.r = (as.matrix(dist(dbase))[-1, 1])
        if (mqgwrweight == TRUE)
          w.new = gwr.gauss((dist.r) ^ 2, band)
        if (mqgwrweight == FALSE)
          w.new = gwr.bisquare((dist.r) ^ 2, band)
        w.new = diag(w.new)
        ob1 = QRLM(
          x.design,
          y,
          maxit = 100,
          q = c(saq[i, 1]),
          w1 = w.new,
          k = k.value
        )
        coef <-
          matrix(c(t(ob1$coef)), nrow = 1, ncol = p) # need to be ordered by area
        coef <- t(coef)
        W_star = diag(c(ob1$q.weight), n, n)
        S = W_star %*% x.design %*% solve(t(x.design) %*% W_star %*% x.design)
        xir = x.r.area[j, ]
        tmp = tmp + S %*% xir
        pred.medr[j] <- (x.r.area[j, ] %*% coef[, 1])
      }
      
      pred.meds = 0
      sj <- matrix(x.design[area == id.area.in[i]], ni[i], p)
      
      lon.gwr = (lon)[area == id.area.in[i]]
      lat.gwr = (lat)[area == id.area.in[i]]
      for (j in 1:ni[i]) {
        dbase = as.matrix(rbind(
          cbind(lon.gwr[j], lat.gwr[j]),
          cbind(as.vector(lon), as.vector(lat))
        ))
        dist.r = (as.matrix(dist(dbase))[-1, 1])
        if (mqgwrweight == TRUE)
          w.new = gwr.gauss((dist.r) ^ 2, band)
        if (mqgwrweight == FALSE)
          w.new = gwr.bisquare((dist.r) ^ 2, band)
        w.new = diag(w.new)
        ob1 = QRLM(
          x.design,
          y,
          maxit = 100,
          q = c(saq[i, 1]),
          w1 = w.new,
          k = k.value
        )
        coef <-
          matrix(c(t(ob1$coef)), nrow = 1, ncol = p) # need to be ordered by area
        coef <- t(coef)
        W_star = diag(c(ob1$q.weight), n, n)
        S = W_star %*% x.design %*% solve(t(x.design) %*% W_star %*% x.design)
        xis = sj[j, ]
        tmp1 = tmp1 + S %*% xis
        pred.meds[j] <- (sj[j, ] %*% coef[, 1])
      }
      
      
      f1 <- y[area == id.area.in[i]]
      res.s <- c(res.s, (f1 - pred.meds))
      
      ir = rep(0, n)
      ir[area == id.area.in[i]] <- 1
      welsh.cd = ir + ir * ((ri[i] + ni[i]) / ni[i]) + tmp - ((ri[i] +
                                                                 ni[i]) / ni[i]) * tmp1
      data <- cbind(as.vector(area), welsh.cd)
      
      for (kk in 1:n)
      {
        if (data[kk, 1] == id.area.in[i])
          ci[kk, i, 1] = data[kk, 2] - 1
        else if (data[kk, 1] != id.area.in[i])
          ci[kk, i, 1] = data[kk, 2]
      }
      RI.MQGWR_CD.Mean[i] = as.real(1 / (ri[i] + ni[i]) * (t(welsh.cd) %*%
                                                             as.vector(y)))
    }
    res.s = res.s ^ 2
    res.d = cbind(res.s, as.vector(area), ci[, , 1])
    
    for (oo in 1:m)
    {
      mse[oo] = (1 / (ni[oo] + ri[oo]) ^ 2) * (sum((res.d[, (oo + 2)][res.d[, 2] ==
                                                                        oo] ^ 2 + (ri[oo]) / n) * res.d[, 1][res.d[, 2] == oo]) + sum(res.d[, (oo +
                                                                                                                                                 2)][res.d[, 2] != oo] ^ 2 * res.d[, 1][res.d[, 2] != oo]))
    }
  }
  
  if (method == "mqgwr-li")
    
  {
    if (mqgwrweight == TRUE)
      band = gwr.sel(y ~ 1 + x, coords = cbind(lon, lat), gweight = gwr.gauss)
    
    if (mqgwrweight == FALSE)
      band = gwr.sel(y ~ 1 + x, coords = cbind(lon, lat), gweight = gwr.bisquare)
    
    if (mqgwrweight == TRUE)
      w.sp <- gwr.gauss((eu_dist ^ 2), band)
    if (mqgwrweight == FALSE)
      w.sp <- gwr.bisquare((eu_dist ^ 2), band)
    q.new = 0
    n.q = length(q.value)
    fitted = matrix(0, n, n.q)
    for (qj in 1:n.q)
    {
      ob.trad <- QRLM1(x.design, y, maxit = 100, q = q.value[qj])
      for (ii in 1:n) {
        w.new <- (w.sp[, ii])
        err <-
          sum(w.new * ob.trad$q.weights * ob.trad$residuals) / sum(w.new * ob.trad$q.weights)
        fitted[ii, qj] = ob.trad$fitted.values[ii] + err
      }
    }
    q.new <-
      matrix(c(gridfitinter(y, fitted, q.value)), nrow = n, ncol = 1)
    
    qmat1 <- matrix(c(q.new, area), nrow = n, ncol = 2)
    
    mqo1 = tapply(qmat1[, 1], qmat1[, 2], mean)
    
    saq <- matrix(0, nrow = m, ncol = 2)
    
    saq[, 1] = mqo1
    
    saq[, 2] = sort(unique(qmat1[, 2]))
    
    ci = array(rep(0, n * m), dim = c(n, m, 1))
    
    res.s = NULL
    
    for (i in 1:m)
    {
      pred.medr = 0
      x.r.area <-
        matrix(cbind(1, x.r)[area.r == id.area.in[i]], ri[i], p)
      lon.gwr = lon.r[area.r == id.area.in[i]]
      lat.gwr = lat.r[area.r == id.area.in[i]]
      
      tmp = matrix(0, 1, n)
      tmp1 = matrix(0, 1, n)
      
      ob.trad <- QRLM1(x.design, y, maxit = 100, q = c(saq[i, 1]))
      coef <-
        matrix(c(t(ob.trad$coef)), nrow = 1, ncol = 2) # need to be ordered by area
      coef <- t(coef)
      wd <- diag(c(ob.trad$q.weights))
      meat <- wd %*% x.design %*% solve(t(x.design) %*% wd %*% x.design)
      meat1 = (diag(1, n, n) - x.design %*% solve(t(x.design) %*% wd %*%
                                                    x.design) %*% t(x.design) %*% wd)
      
      for (j in 1:ri[i]) {
        dbase = as.matrix(rbind(
          cbind(lon.gwr[j], lat.gwr[j]),
          cbind(as.vector(lon), as.vector(lat))
        ))
        dist.r = (as.matrix(dist(dbase))[-1, 1])
        if (mqgwrweight == TRUE)
          w.new = gwr.gauss((dist.r) ^ 2, band)
        if (mqgwrweight == FALSE)
          w.new = gwr.bisquare((dist.r) ^ 2, band)
        
        uno = matrix(1, n, 1)
        err1 = ((t(uno) %*% (diag(
          c(w.new), n, n
        ) %*% wd %*% meat1)) * as.real(solve(
          t(uno) %*% diag(c(w.new), n, n) %*% wd %*% uno
        )))
        tmp = tmp + err1
        pred.medr[j] <-
          (x.r.area[j, ] %*% coef[, 1] + as.real(err1 %*% matrix(y, n, 1)))
      }
      
      pred.meds = 0
      sj <- matrix(x.design[area == id.area.in[i]], ni[i], p)
      
      lon.gwr = (lon)[area == id.area.in[i]]
      lat.gwr = (lat)[area == id.area.in[i]]
      for (j in 1:ni[i]) {
        dbase = as.matrix(rbind(
          cbind(lon.gwr[j], lat.gwr[j]),
          cbind(as.vector(lon), as.vector(lat))
        ))
        dist.r = (as.matrix(dist(dbase))[-1, 1])
        if (mqgwrweight == TRUE)
          w.new = gwr.gauss((dist.r) ^ 2, band)
        if (mqgwrweight == FALSE)
          w.new = gwr.bisquare((dist.r) ^ 2, band)
        
        err2 = ((t(uno) %*% (diag(
          c(w.new), n, n
        ) %*% wd %*% meat1)) * as.real(solve(
          t(uno) %*% diag(c(w.new), n, n) %*% wd %*% uno
        )))
        tmp1 = tmp1 + err2
        pred.meds[j] <-
          (sj[j, ] %*% coef[, 1] + as.real(err2 %*% matrix(y, n, 1)))
      }
      
      
      f1 <- y[area == id.area.in[i]]
      res.s <- c(res.s, (f1 - pred.meds))
      
      sj.tot <- apply(sj, 2, sum)
      rj.tot <- apply(x.r.area, 2, sum)
      
      ir = rep(0, n)
      ir[area == id.area.in[i]] <- 1
      
      
      welsh.cd = ir * ((ri[i] + ni[i]) / ni[i]) + meat %*% (rj.tot - (ri[i] /
                                                                        ni[i]) * sj.tot) + t(tmp) - ((ri[i]) / ni[i]) * t(tmp1)
      data <- cbind(as.vector(area), welsh.cd)
      
      for (kk in 1:n)
      {
        if (data[kk, 1] == id.area.in[i])
          ci[kk, i, 1] = data[kk, 2] - 1
        else if (data[kk, 1] != id.area.in[i])
          ci[kk, i, 1] = data[kk, 2]
      }
      RI.MQGWR_CD.Mean[i] = as.real(1 / (ri[i] + ni[i]) * (t(welsh.cd) %*%
                                                             as.vector(y)))
    }
    res.s = res.s ^ 2
    res.d = cbind(res.s, as.vector(area), ci[, , 1])
    
    for (oo in 1:m)
    {
      mse[oo] = (1 / (ni[oo] + ri[oo]) ^ 2) * (sum((res.d[, (oo + 2)][res.d[, 2] ==
                                                                        oo] ^ 2 + (ri[oo]) / n) * res.d[, 1][res.d[, 2] == oo]) + sum(res.d[, (oo +
                                                                                                                                                 2)][res.d[, 2] != oo] ^ 2 * res.d[, 1][res.d[, 2] != oo]))
    }
    
  }
  #out of sample area
  if (m.r > m) {
    rr.sample = m.r - m
    for (i in 1:rr.sample)
    {
      Ri = sum(area.r == id.area.out[i])
      x.r.area <- matrix(cbind(1, x.r)[area.r == id.area.out[i]], Ri, p)
      pred.medr = 0
      lon.gwr = lon.r[area.r == id.area.out[i]]
      lat.gwr = lat.r[area.r == id.area.out[i]]
      for (j in 1:Ri) {
        dbase = as.matrix(rbind(
          cbind(lon.gwr[j], lat.gwr[j]),
          cbind(as.vector(lon), as.vector(lat))
        ))
        dist.r = (as.matrix(dist(dbase))[-1, 1])
        if (mqgwrweight == TRUE)
          w.new = gwr.gauss((dist.r) ^ 2, band)
        if (mqgwrweight == FALSE)
          w.new = gwr.bisquare((dist.r) ^ 2, band)
        w.new = diag(w.new)
        ob1 = QRLM(x.design,
                   y,
                   maxit = 100,
                   q = 0.5,
                   w1 = w.new)
        coef <-
          matrix(c(t(ob1$coef)), nrow = 1, ncol = 2) # need to be ordered by area
        coef <- t(coef)
        pred.medr[j] <- (x.r.area[j, ] %*% coef[, 1])
      }
      
      
      RI.MQGWR_CD.Mean.out[i] <- 1 / (Ri) * (sum(pred.medr))
    }
  }
  
  list(
    Area.code.in = id.area.in,
    Area.code.out = id.area.out,
    Est.Mean.in = RI.MQGWR_CD.Mean,
    Est.Mean.out = RI.MQGWR_CD.Mean.out,
    Est.mse.in = mse
  )
  
}
