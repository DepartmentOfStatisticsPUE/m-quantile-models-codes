#############################################
# M-quantiles models				        #
# for poverty indicators: hcr, poverty-gap  #
#############################################


rm(list = ls(all.names = TRUE))


set.seed(1977)


###########LIBRARY
library(nlme)
library(MASS)
library(np)


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



mq.coef <- function(myx, myy, myregioncode, maxiter = 100) {
  #This function estimate the m-quantile regression coefficients
  #myx<- x sample matrix of auxiliary variables
  #myy<- y vector
  #mynumauxvar<- number of auxiliary variables (include constant)
  #myregioncode<- area code for y and x units, data must be ordered by area code
  #maxiter<- OPTIONAL, number of maximum iteration for ob algorithm

  myar <- unique(myregioncode)
  myareas <- length(myar)
  mysamplesize <- sum(as.numeric(table(myregioncode)))
  mynumauxvar <- dim(myx)[2]

  ob <-
    QRLM(myx,
         myy,
         maxit = maxiter,
         q = sort(c(
           seq(0.006, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98
         )),
         k = 1.345)
  qo <-
    matrix(c(gridfitinter(myy, ob$fitted.values, ob$q.values)), nrow = mysamplesize, ncol =
             1)
  qmat <- matrix(c(qo, myregioncode),
                 nrow = length(myregioncode),
                 ncol = 2)
  mqo <- aggregate(qmat[, 1], list(d2 = qmat[, 2]), mean)[, 2]
  saq <- matrix(c(mqo, myar), nrow = myareas, ncol = 2)
  saq <- rbind(saq, c(0.5, 9999))
  ob1 <- QRLM(myx,
              myy,
              maxit = maxiter,
              q = c(mqo[1:myareas]),
              k = 1.345)
  mycoef <-
    matrix(c(t(ob1$coefficients)), nrow = myareas, ncol = mynumauxvar) # need to be ordered by area
  mycoef <- t(mycoef)
  list(q.mean = mycoef, q.unit = qmat)
}



compute.hcr.pg <-
  function(my.ys,
           my.x.s,
           my.X.pop,
           myregioncode,
           myregioncodepop,
           L,
           areas,
           ar,
           pop.size,
           sample.size,
           myz) {
    #f.EB.0<-array(0,dim=c(areas))
    #f.EB.1<-array(0,dim=c(areas))
    f.MQ.0 <- array(0, dim = c(areas))
    f.MQ.1 <- array(0, dim = c(areas))

    res.mq <- NULL

    #Fit the model M-quantile
    tmp = mq.coef(my.x.s, my.ys, myregioncode)
    beta.mq <- tmp$q.mean
    for (i in 1:areas) {
      #MQ
      ysd <- my.ys[myregioncode == ar[i]]
      x.sd <- my.x.s[myregioncode == ar[i], ]
      Eds.mq <- x.sd %*% beta.mq[, i]
      res.d.mq <- ysd - Eds.mq
      res.mq <- c(res.mq, res.d.mq)
    }


    for (i in 1:areas) {
      #	ysd<-my.ys[myregioncode==ar[i]]
      #	x.sd<-my.x.s[myregioncode==ar[i],]
      x.rd <- my.X.pop[myregioncodepop == ar[i], ]
      #	I.s<-ysd<z
      #	Ed.s<-ysd

      #Monte Carlo approximation to the best predictor of yi
      F.0.hl.mq <- matrix(0, L, 1)
      F.1.hl.mq <- matrix(0, L, 1)

      for (l in 1:L) {
        #MQ
        Ehl.mc <- x.rd %*% beta.mq[, i] + sample(res.mq, pop.size[i], replace =
                                                   TRUE)
        I.mc.mq <- Ehl.mc < myz
        F.0.hl.mq[l, 1] <- sum(I.mc.mq) / pop.size[i]
        Ehl.mc[Ehl.mc < 0] <- 0
        F.1.hl.mq[l, 1] <- (1 / pop.size[i]) * sum((1 - Ehl.mc / myz) * I.mc.mq)
      }
      #	f.EB.0[i]<-mean(F.0.hl[,1])
      #	f.EB.1[i]<-mean(F.1.hl[,1])
      f.MQ.0[i] <- mean(F.0.hl.mq[, 1])
      f.MQ.1[i] <- mean(F.1.hl.mq[, 1])

    }#i ends here

    f.MQ.0[which(f.MQ.0 > 1)] <- 1
    f.MQ.1[which(f.MQ.1 > 1)] <- 1

    res <- list(
      HCR.MQ = f.MQ.0,
      PG.MQ = f.MQ.1,
      res.mq = res.mq,
      mycoef = beta.mq
    )
    res

  }#Function compute.hcr.pg ends here



MQ.SAE.poverty <-
  function(my.ys,
           my.x.s,
           my.X.pop,
           myregioncode,
           myregioncodepop,
           L = 50,
           myMSE = TRUE,
           myB = 1,
           myR = 400,
           method = "eu",
           pov.l = NULL) {
    #This function compute the HCR and the PovertyGap statistics for Small Area


    areas <- length(unique(myregioncode))
    ar <- unique(myregioncode)
    pop.size <- table(myregioncodepop)
    sample.size <- table(myregioncode)
    myid <- 1:sum(pop.size)
    if (is.null(pov.l))
      z <- 0.6 * median(my.ys)
    if (!is.null(pov.l))
      z <- povl

    myq.true.boot <- array(0, dim = c(2, areas, myB))
    myarea.q.boot.r <- array(0, dim = c(2, areas, myB, myR))
    myarea.q.boot <- array(0, dim = c(2, areas, myB))
    myq.true.boot.m <- array(0, dim = c(2, areas))
    myarea.q.boot.m <- array(0, dim = c(2, areas))
    BIAS.boot <- matrix(0, 2, areas)
    VAR.boot <- matrix(0, 2, areas)
    MSE.boot <- matrix(0, 2, areas)
    CI.boot.hcr <- matrix(0, areas, 2)
    CI.boot.pg <- matrix(0, areas, 2)

    estimate <-
      compute.hcr.pg(
        my.ys,
        my.x.s,
        my.X.pop,
        myregioncode,
        myregioncodepop,
        L,
        areas,
        ar,
        pop.size,
        sample.size,
        z
      )
    res.mq <- estimate$res.mq


    if (myMSE) {
      #Generate B bootstrap Population (size N)

      if (method == "sc") {
        #Centering residuals in each areas (use this for area conditioned approach)
        res.s.centered <- NULL
        for (i in 1:areas) {
          res.s.centered[i] <-
            list(res.mq[myregioncode == ar[i]] - mean(res.mq[myregioncode == ar[i]]))
        }

        #smoothed density of residuals areas conditioned
        Fhat.ord <- NULL
        res.ord <- NULL
        for (i in 1:areas) {
          bw <- npudensbw( ~ res.s.centered[[i]], ckertype = "epanechnikov")
          Fhat <- fitted(npudist(bws = bw))
          res.ord[i] <- list(sort(res.s.centered[[i]]))
          Fhat.ord[i] <- list(sort(Fhat))
        }
      }

      if (method == "su") {
        #Centering residuals for the whole sample (use this for area unconditioned approach)
        res.s.centered <- sort(res.mq - mean(res.mq))

        #smoothed density of residuals areas unconditioned
        Fhat.ord <- NULL
        bw <- npudensbw( ~ res.s.centered, ckertype = "epanechnikov")
        Fhat <- fitted(npudist(bws = bw))
        Fhat.ord <- sort(Fhat)
      }

      if (method == "ec") {
        #Centering residuals in each areas (use this for area conditioned approach)
        res.s.centered <- NULL
        for (i in 1:areas) {
          res.s.centered[i] <-
            list(res.mq[myregioncode == ar[i]] - mean(res.mq[myregioncode == ar[i]]))
        }
      }

      if (method == "eu") {
        #Centering residuals for the whole sample (use this for area unconditioned approach)
        res.s.centered <- sort(res.mq - mean(res.mq))
      }

      for (b in 1:myB) {
        if (method == "sc") {
          #Sample from kernel density areas conditioned
          samp.boot <- NULL
          for (i in 1:areas) {
            s.boot <- NULL
            for (g in 1:pop.size[i]) {
              s.boot[g] <-
                which(Fhat.ord[[i]] == quantile(Fhat.ord[[i]], prob = runif(1), type = 3))
            }
            samp.boot[i] <- list(s.boot)
          }
          #Population smoothed density of residuals area conditioned
          y.boot <- NULL
          y.boot.i <- NULL
          for (i in 1:areas) {
            y.boot.i <-
              my.X.pop[myregioncodepop == ar[i], ] %*% estimate$mycoef[, i] + res.ord[[i]][samp.boot[[i]]]
            y.boot <- c(y.boot, y.boot.i)
          }
        }

        if (method == "su") {
          #Sample from kernel density areas unconditioned
          samp.boot <- NULL
          for (i in 1:areas) {
            s.boot <- NULL
            for (g in 1:pop.size[i]) {
              s.boot[g] <- which(Fhat.ord == quantile(Fhat.ord, prob = runif(1), type =
                                                        3))
            }
            samp.boot[i] <- list(s.boot)
          }
          #Population smoothed density of residuals areas unconditioned
          y.boot <- NULL
          y.boot.i <- NULL
          for (i in 1:areas) {
            y.boot.i <-
              my.X.pop[myregioncodepop == ar[i], ] %*% estimate$mycoef[, i] + res.s.centered[samp.boot[[i]]]
            y.boot <- c(y.boot, y.boot.i)
          }
        }

        if (method == "ec") {
          #Population empirical density of residuals area conditioned
          y.boot <- NULL
          y.boot.i <- NULL
          for (i in 1:areas) {
            y.boot.i <-
              my.X.pop[myregioncodepop == ar[i], ] %*% estimate$mycoef[, i] + sample(res.s.centered[[i]], pop.size[i], replace =
                                                                                       TRUE)
            y.boot <- c(y.boot, y.boot.i)
          }
        }

        if (method == "eu") {
          #Population empirical density of residuals area unconditioned
          y.boot <- NULL
          y.boot.i <- NULL
          for (i in 1:areas) {
            y.boot.i <-
              my.X.pop[myregioncodepop == ar[i], ] %*% estimate$mycoef[, i] + sample(res.s.centered, pop.size[i], replace =
                                                                                       TRUE)
            y.boot <- c(y.boot, y.boot.i)
          }
        }


        for (ii in 1:areas) {
          y.d.boot <- y.boot[myregioncodepop == ar[ii]]
          myq.true.boot[1, ii, b] <- sum(y.d.boot < z) / pop.size[ii]
          y.d.boot[y.d.boot < 0] <- 0
          myq.true.boot[2, ii, b] <-
            (1 / pop.size[ii]) * sum((1 - y.d.boot / z) * (y.d.boot < z))
        }


        for (rr in 1:myR) {
          mysboot <- NULL
          s.boot.i <- NULL
          for (ii in 1:areas) {
            s.boot.i <- sample(myid[myregioncodepop == ar[ii]], sample.size[ii])
            mysboot <- c(mysboot, s.boot.i)
          }
          ys.boot <- y.boot[mysboot]
          x.s.boot <- my.X.pop[mysboot, ]
          estimate.boot <-
            compute.hcr.pg(
              ys.boot,
              x.s.boot,
              my.X.pop,
              myregioncode,
              myregioncodepop,
              L,
              areas,
              ar,
              pop.size,
              sample.size,
              z
            )
          myarea.q.boot.r[1, , b, rr] <- estimate.boot$HCR.MQ
          myarea.q.boot.r[2, , b, rr] <- estimate.boot$PG.MQ
        }


        for (ii in 1:areas) {
          myarea.q.boot[1, ii, b] <- mean(myarea.q.boot.r[1, ii, b, ])
          myarea.q.boot[2, ii, b] <- mean(myarea.q.boot.r[2, ii, b, ])
        }

      }#B ends here

      for (i in 1:areas) {
        myq.true.boot.m[1, i] <- mean(myq.true.boot[1, i, ])
        myq.true.boot.m[2, i] <- mean(myq.true.boot[2, i, ])
        myarea.q.boot.m[1, i] <- mean(myarea.q.boot[1, i, ])
        myarea.q.boot.m[2, i] <- mean(myarea.q.boot[2, i, ])
      }

      for (i in 1:areas) {
        BIAS.boot[1, i] <- myarea.q.boot.m[1, i] - myq.true.boot.m[1, i]
        BIAS.boot[2, i] <- myarea.q.boot.m[2, i] - myq.true.boot.m[2, i]
        aux.0 <- matrix(0, myB, 1)
        aux.1 <- matrix(0, myB, 1)
        for (b in 1:myB) {
          aux.0[b, 1] <-
            (1 / myR) * sum((myarea.q.boot.r[1, i, b, ] - myarea.q.boot[1, i, b]) ^
                              2)
          aux.1[b, 1] <-
            (1 / myR) * sum((myarea.q.boot.r[2, i, b, ] - myarea.q.boot[2, i, b]) ^
                              2)
        }
        VAR.boot[1, i] <- (1 / myB) * sum(aux.0[, 1])
        VAR.boot[2, i] <- (1 / myB) * sum(aux.1[, 1])
      }

      for (i in 1:areas) {
        MSE.boot[1, i] <- ((BIAS.boot[1, i]) ^ 2) + VAR.boot[1, i]
        MSE.boot[2, i] <- ((BIAS.boot[2, i]) ^ 2) + VAR.boot[2, i]
        CI.boot.hcr[i, ] <-
          quantile(c(myarea.q.boot.r[1, i, , ]), prob = c(0.025, 0.975))
        CI.boot.pg[i, ] <-
          quantile(c(myarea.q.boot.r[2, i, , ]), prob = c(0.025, 0.975))
      }

      rmse.hcr.mq <- sqrt(MSE.boot[1, ])
      rmse.pg.mq <- sqrt(MSE.boot[2, ])

    }#end if myMSE

    if (myMSE == FALSE) {
      rmse.hcr.mq <- NULL
      rmse.pg.mq <- NULL
    }

    res <-
      list(
        HCR.MQ = estimate$HCR.MQ,
        PG.MQ = estimate$PG.MQ,
        RMSE.HCR.MQ = rmse.hcr.mq,
        RMSE.PG.MQ = rmse.pg.mq,
        Area.Code = unique(myregioncode),
        Pov.Line = z
      )#,CI.HCR=CI.boot.hcr,CI.PG=CI.boot.pg


  }#Function ends here
