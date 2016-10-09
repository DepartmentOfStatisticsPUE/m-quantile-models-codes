rm(list = ls(all.names = TRUE))


set.seed(1977)

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
    QRLM(myx, myy, maxit = maxiter, q = sort(c(
      seq(0.006, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98
    )))
  qo <-
    matrix(c(gridfitinter(myy, ob$fitted.values, ob$q.values)), nrow = mysamplesize, ncol =
             1)
  qmat <- matrix(c(qo, myregioncode),
                 nrow = length(myregioncode),
                 ncol = 2)
  mqo <- aggregate(qmat[, 1], list(d2 = qmat[, 2]), mean)[, 2]
  saq <- matrix(c(mqo, myar), nrow = myareas, ncol = 2)
  saq <- rbind(saq, c(0.5, 9999))
  ob1 <- QRLM(myx, myy, maxit = maxiter, q = c(mqo[1:myareas]))
  mycoef <-
    matrix(c(t(ob1$coefficients)), nrow = myareas, ncol = mynumauxvar) # need to be ordered by area
  mycoef <- t(mycoef)
  mycoef
}


intsolver <-
  function(myqest,
           myyboot,
           myX,
           myregioncodepop,
           mypopsize,
           myar,
           myareas,
           adjseed,
           mysboot,
           mymaxit = 100) {
    myres <- array(0, dim = c(myareas))

    myy <- myyboot[mysboot]
    myx <- myX[mysboot, ]
    myregioncode <- myregioncodepop[mysboot]
    myregioncoder <- myregioncodepop[-mysboot]
    mysamplesizer <- as.numeric(table(myregioncoder))
    myX.r <- myX[-mysboot, ]

    # M-quantiles
    coef.boot <- mq.coef(myx, myy, myregioncode)

    # Quantile Estimation Using Chambers Dunstan Estimator
    for (i in 1:myareas) {
      f1 <- myy[myregioncode == myar[i]]
      X.aux.i <- as.matrix(myX.r[myregioncoder == myar[i], ])
      pred.medr <- (X.aux.i %*% coef.boot[, i])
      x.design.i <- as.matrix(myx[myregioncode == myar[i], ])
      pred.meds <- (x.design.i %*% coef.boot[, i])
      res.s <- f1 - pred.meds
      #z<-matrix(rep(res.s,sample.sizer[i]),nrow=sample.sizer[i],ncol=sample.size[i])
      z <- sample(res.s, mysamplesizer[i], replace = TRUE)
      z <- z + pred.medr
      #	z<-pmax(0,z)
      comb <- c(f1, pred.medr)
      start0 <- quantile(comb, prob = c(myqest))

      sameside <- T
      myiter <- 0
      while (sameside & myiter < mymaxit) {
        ff2 <- sum(c(z) <= start0)
        ff1 <- sum(c(f1) <= start0)
        f0 <- 1 / (mypopsize[i]) * (ff1 + ff2)
        if (f0 <= myqest)
          start1 <- start0 + adjseed
        if (f0 > myqest)
          start1 <- start0 - adjseed
        ff2 <- sum(c(z) <= start1)
        ff1 <- sum(c(f1) <= start1)
        f.new <- 1 / (mypopsize[i]) * (ff1 + ff2)
        start.bef <- start0
        start.aft <- start1
        if (f0 <= myqest & f.new >= myqest)
          sameside <- F
        if (f0 >= myqest & f.new <= myqest)
          sameside <- F
        start0 <- start1
        myiter <- myiter + 1
      }
      if (myiter >= 100)
        warning("intsolver sameside did not converge in ",
                mymaxit,
                " iteration")

      ff2 <- sum(c(z) <= start.bef)
      ff1 <- sum(c(f1) <= start.bef)
      f.bef <- 1 / (mypopsize[i]) * (ff1 + ff2)
      ff2 <- sum(c(z) <= start.aft)
      ff1 <- sum(c(f1) <= start.aft)
      f.aft <- 1 / (mypopsize[i]) * (ff1 + ff2)
      fdif <- abs(f.bef - f.aft)
      if (fdif >= 0.01) {
        start.med <- (start.bef + start.aft) / 2

        eps <- 0.001
        tol <- 50
        myiter <- 0
        while (abs(tol) >= 0.05 & myiter < mymaxit) {
          ff2 <- sum(c(z) <= start.med)
          ff1 <- sum(c(f1) <= start.med)
          fmed <- 1 / (mypopsize[i]) * (ff1 + ff2)
          tol <- (fmed - myqest)
          fmedl <- 1 / (mypopsize[i]) * (ff1 + ff2) - eps
          fmedu <- 1 / (mypopsize[i]) * (ff1 + ff2) + eps
          if (fmed < myqest & fmedl < fmedu)
            start.bef <- start.med
          if (fmed < myqest & fmedl > fmedu)
            start.aft <- start.med
          if (fmed > myqest & fmedl < fmedu)
            start.aft <- start.med
          if (fmed > myqest & fmedl > fmedu)
            start.bef <- start.med
          #if (fmedl<fmedu) start.bef<-start.med
          #if (fmedl>fmedu) start.aft<-start.med
          #if (fmedl<fmedu) start.aft<-start.med
          #if (fmedl>fmedu) start.bef<-start.med
          start.med <- (start.bef + start.aft) / 2
          #print(tol)
          myiter <- myiter + 1
        }
        if (myiter >= 100)
          warning("intsolver fdif>0.01 tol did not converge in ",
                  mymaxit,
                  " iteration")
      }

      if (fdif < 0.01) {
        start.med <- (start.bef + start.aft) / 2

        eps <- 0.001
        tol <- 50
        myiter <- 0
        while (abs(tol) > 0.05 & myiter < mymaxit) {
          ff2 <- sum(c(z) <= start.med)
          ff1 <- sum(c(f1) <= start.med)
          fmed <- 1 / (mypopsize[i]) * (ff1 + ff2)
          tol <- (fmed - myqest)
          fmedl <- 1 / (mypopsize[i]) * (ff1 + ff2) - eps
          fmedu <- 1 / (mypopsize[i]) * (ff1 + ff2) + eps
          if (fmed < myqest & fmedl < fmedu)
            start.bef <- start.med
          if (fmed < myqest & fmedl > fmedu)
            start.aft <- start.med
          if (fmed > myqest & fmedl < fmedu)
            start.aft <- start.med
          if (fmed > myqest & fmedl > fmedu)
            start.bef <- start.med
          #if (fmedl<fmedu) start.bef<-start.med
          #if (fmedl>fmedu) start.aft<-start.med
          #if (fmedl<fmedu) start.aft<-start.med
          #if (fmedl>fmedu) start.bef<-start.med
          start.med <- (start.bef + start.aft) / 2
          #print(tol)
          myiter <- myiter + 1
        }
        if (myiter >= 100)
          warning("intsolver fdif<0.01 tol did not converge in ",
                  mymaxit,
                  " iteration")
      }

      myres[i] <- start.med
      #print(i)
    }#Iteration i in bootstrap ends here
    myres
  }


boot.CD.R <-
  function(myqest,
           myyboot,
           myX,
           myregioncodepop,
           mypopsize,
           mysamplesize,
           myar,
           myareas,
           myadjseed,
           myR,
           myid) {
    myproc.a <- array(0, dim = c(myareas, myR))

    #Sampling from bootstrap population
    for (r in 1:myR) {
      mysboot <- NULL
      s.boot.i <- NULL
      for (i in 1:myareas) {
        s.boot.i <- sample(myid[myregioncodepop == myar[i]], mysamplesize[i])
        mysboot <- c(mysboot, s.boot.i)
      }

      myproc.a[, r] <-
        intsolver(
          myqest,
          myyboot,
          myX,
          myregioncodepop,
          mypopsize,
          myar,
          myareas,
          myadjseed,
          mysboot
        )
    }#R ends here
    myproc.a
  }


MQ.SAE.quant <-
  function(myqgrid,
           myy,
           myx,
           myX,
           myregioncode,
           myregioncodepop,
           adjseed = max(0.15, mean(myy) / 500),
           myMSE = FALSE,
           B = 1,
           R = 400,
           method = "su",
           mymaxit = 100) {
    #This function estimate quantiles via CD estimator when n/N -> p with p very small
    #myqgrid<- quantiles order to be estimated (i.e. 0.25,0.50,0.75)
    #myy<- y vector
    #myx<- x sample matrix of auxiliary variables
    #myX<- X population matrix of auxiliary variables
    #myregioncode<- area code for y and x units, data must be ordered by area code
    #myregioncodepop<- area code for X units (population), data must be ordered by area code
    #adjseed<- OPTIONAL, tune the value used to find two good starting point to solve the integral
    #myMSE<-TRUE compute the MSE of the CD quantiles estimate via bootstrap, FALSE does not compute any MSE
    #B<- number of bootstrap population
    #R<- number of bootstrap samples
    #method<- which method to be used to estimate residuals distribution, choice= "su" (smooth unconditional),"eu" (emprirical unconditional),"sc" (smooth conditional),"ec" (empirical coditional)
    #mymaxit: maximum iteration allowed in the while routines

    myar <- unique(myregioncode)
    myareas <- length(myar)
    mypopsize <- as.numeric(table(myregioncodepop))
    mysamplesize <- as.numeric(table(myregioncode))
    myquantnum <- length(myqgrid)
    id <- seq(1:sum(mypopsize))

    myarea.q <- matrix(0, myquantnum, myareas)
    myq.true.boot <- array(0, dim = c(myquantnum, myareas, B))
    myarea.q.boot.r <- array(0, dim = c(myquantnum, myareas, B, R))
    myarea.q.boot <- array(0, dim = c(myquantnum, myareas, B))
    myq.true.boot.m <- array(0, dim = c(myquantnum, myareas))
    myarea.q.boot.m <- array(0, dim = c(myquantnum, myareas))
    BIAS.boot <- matrix(0, myquantnum, myareas)
    VAR.boot <- matrix(0, myquantnum, myareas)
    MSE.boot <- matrix(0, myquantnum, myareas)
    kk <- 0

    mycoef <- mq.coef(myx, myy, myregioncode)

    for (qq in myqgrid) {
      qest <- qq
      kk <- kk + 1
      myres <- NULL

      for (i in 1:myareas) {
        f1 <- myy[myregioncode == myar[i]]
        X.aux.i <- as.matrix(myX[myregioncodepop == myar[i], ])
        pred.medtot <- (X.aux.i %*% mycoef[, i])
        x.design.i <- as.matrix(myx[myregioncode == myar[i], ])
        pred.meds <- (x.design.i %*% mycoef[, i])
        res.s <- f1 - pred.meds
        myres[i] <- list(res.s)
        #			z<-matrix(rep(res.s,pop.size[i]),nrow=pop.size[i],ncol=sample.size[i])
        z <- sample(res.s, mypopsize[i], replace = TRUE)
        z <- z + pred.medtot
        #			z<-pmax(0,z)
        comb <- c(pred.medtot)
        start0 <- quantile(comb, prob = c(qest))

        sameside <- T
        myiter <- 0
        while (sameside & myiter < mymaxit) {
          ff2 <- sum(c(z) <= start0)
          #				ff1<-sum(c(f1)<=start0)
          f0 <- 1 / (mypopsize[i]) * (ff2)
          if (f0 <= qest)
            start1 <- start0 + adjseed
          if (f0 > qest)
            start1 <- start0 - adjseed
          ff2 <- sum(c(z) <= start1)
          #				ff1<-sum(c(f1)<=start1)
          f.new <- 1 / (mypopsize[i]) * (ff2)
          start.bef <- start0
          start.aft <- start1
          if (f0 <= qest & f.new >= qest)
            sameside <- F
          if (f0 >= qest & f.new <= qest)
            sameside <- F
          start0 <- start1
          myiter <- myiter + 1
        }
        if (myiter >= 100)
          warning("CD.quant sameside did not converge in ",
                  mymaxit,
                  " iteration")

        ff2 <- sum(c(z) <= start.bef)
        #			ff1<-sum(c(f1)<=start.bef)
        f.bef <- 1 / (mypopsize[i]) * (ff2)
        ff2 <- sum(c(z) <= start.aft)
        #			ff1<-sum(c(f1)<=start.aft)
        f.aft <- 1 / (mypopsize[i]) * (ff2)
        fdif <- abs(f.bef - f.aft)
        if (fdif >= 0.01) {
          start.med <- (start.bef + start.aft) / 2

          eps <- 0.001
          tol <- 50
          myiter <- 0
          while (abs(tol) >= 0.05 & myiter < mymaxit) {
            ff2 <- sum(c(z) <= start.med)
            #					ff1<-sum(c(f1)<=start.med)
            fmed <- 1 / (mypopsize[i]) * (ff2)
            tol <- (fmed - qest)
            fmedl <- 1 / (mypopsize[i]) * (ff2) - eps
            fmedu <- 1 / (mypopsize[i]) * (ff2) + eps
            if (fmed < qest & fmedl < fmedu)
              start.bef <- start.med
            if (fmed < qest & fmedl > fmedu)
              start.aft <- start.med
            if (fmed > qest & fmedl < fmedu)
              start.aft <- start.med
            if (fmed > qest & fmedl > fmedu)
              start.bef <- start.med
            #if (fmedl<fmedu) start.bef<-start.med
            #if (fmedl>fmedu) start.aft<-start.med
            #if (fmedl<fmedu) start.aft<-start.med
            #if (fmedl>fmedu) start.bef<-start.med
            start.med <- (start.bef + start.aft) / 2
            #print(tol)
            myiter <- myiter + 1
          }
          if (myiter >= 100)
            warning("CD.quant fdif>0.01 tol did not converge in ",
                    mymaxit,
                    " iteration")
        }

        if (fdif < 0.01) {
          start.med <- (start.bef + start.aft) / 2

          eps <- 0.001
          tol <- 50
          myiter <- 0
          while (abs(tol) > 0.05 & myiter < mymaxit) {
            ff2 <- sum(c(z) <= start.med)
            #					ff1<-sum(c(f1)<=start.med)
            fmed <- 1 / (mypopsize[i]) * (ff2)
            tol <- (fmed - qest)
            fmedl <- 1 / (mypopsize[i]) * (ff2) - eps
            fmedu <- 1 / (mypopsize[i]) * (ff2) + eps
            if (fmed < qest & fmedl < fmedu)
              start.bef <- start.med
            if (fmed < qest & fmedl > fmedu)
              start.aft <- start.med
            if (fmed > qest & fmedl < fmedu)
              start.aft <- start.med
            if (fmed > qest & fmedl > fmedu)
              start.bef <- start.med
            #if (fmedl<fmedu) start.bef<-start.med
            #if (fmedl>fmedu) start.aft<-start.med
            #if (fmedl<fmedu) start.aft<-start.med
            #if (fmedl>fmedu) start.bef<-start.med
            start.med <- (start.bef + start.aft) / 2
            #print(tol)
            myiter <- myiter + 1
          }
          if (myiter >= 100)
            warning("CD.quant fdif<0.01 tol did not converge in ",
                    mymaxit,
                    " iteration")
        }
        #			print(i)
        myarea.q[kk, i] <- start.med
      }#Iteration i ends here

      if (myMSE) {
        ######Bootstrap######

        #Generate B bootstrap Population (size N)

        if (method == "sc") {
          #Centering residuals in each areas (use this for area conditioned approach)
          res.s.centered <- NULL
          for (i in 1:myareas) {
            res.s.centered[i] <- list(myres[[i]] - mean(myres[[i]]))
          }

          #smoothed density of residuals areas conditioned
          Fhat.ord <- NULL
          res.ord <- NULL
          for (i in 1:myareas) {
            bw <- npudensbw( ~ res.s.centered[[i]], ckertype = "epanechnikov")
            Fhat <- fitted(npudist(bws = bw))
            res.ord[i] <- list(sort(res.s.centered[[i]]))
            Fhat.ord[i] <- list(sort(Fhat))
          }
        }

        if (method == "su") {
          #Centering residuals for the whole sample (use this for area unconditioned approach)
          res.s.centered <- NULL
          for (i in 1:myareas) {
            res.s.centered <- c(res.s.centered, myres[[i]])
          }
          res.s.centered <- sort(res.s.centered - mean(res.s.centered))

          #smoothed density of residuals areas unconditioned
          Fhat.ord <- NULL
          bw <- npudensbw( ~ res.s.centered, ckertype = "epanechnikov")
          Fhat <- fitted(npudist(bws = bw))
          Fhat.ord <- sort(Fhat)
        }

        if (method == "ec") {
          #Centering residuals in each areas (use this for area conditioned approach)
          res.s.centered <- NULL
          for (i in 1:myareas) {
            res.s.centered[i] <- list(myres[[i]] - mean(myres[[i]]))
          }
        }

        if (method == "eu") {
          #Centering residuals for the whole sample (use this for area unconditioned approach)
          res.s.centered <- NULL
          for (i in 1:myareas) {
            res.s.centered <- c(res.s.centered, myres[[i]])
          }
          res.s.centered <- sort(res.s.centered - mean(res.s.centered))
        }

        for (b in 1:B) {
          #				cat(date(),"Generating Bootstrap Population ",b,"\n")

          if (method == "sc") {
            #Sample from kernel density areas conditioned
            samp.boot <- NULL
            for (i in 1:myareas) {
              s.boot <- NULL
              for (g in 1:mypopsize[i]) {
                s.boot[g] <-
                  which(Fhat.ord[[i]] == quantile(Fhat.ord[[i]], prob = runif(1), type = 3))
              }
              samp.boot[i] <- list(s.boot)
            }
            #Population smoothed density of residuals area conditioned
            y.boot <- NULL
            y.boot.i <- NULL
            for (i in 1:myareas) {
              y.boot.i <-
                myX[myregioncodepop == myar[i], ] %*% mycoef[, i] + res.ord[[i]][samp.boot[[i]]]
              y.boot <- c(y.boot, y.boot.i)
            }
          }

          if (method == "su") {
            #Sample from kernel density areas unconditioned
            samp.boot <- NULL
            for (i in 1:myareas) {
              s.boot <- NULL
              for (g in 1:mypopsize[i]) {
                s.boot[g] <- which(Fhat.ord == quantile(Fhat.ord, prob = runif(1), type =
                                                          3))
              }
              samp.boot[i] <- list(s.boot)
            }
            #Population smoothed density of residuals areas unconditioned
            y.boot <- NULL
            y.boot.i <- NULL
            for (i in 1:myareas) {
              y.boot.i <-
                myX[myregioncodepop == myar[i], ] %*% mycoef[, i] + res.s.centered[samp.boot[[i]]]
              y.boot <- c(y.boot, y.boot.i)
            }
          }

          if (method == "ec") {
            #Population empirical density of residuals area conditioned
            y.boot <- NULL
            y.boot.i <- NULL
            for (i in 1:myareas) {
              y.boot.i <-
                myX[myregioncodepop == myar[i], ] %*% mycoef[, i] + sample(res.s.centered[[i]], mypopsize[i], replace =
                                                                             TRUE)
              y.boot <- c(y.boot, y.boot.i)
            }
          }

          if (method == "eu") {
            #Population empirical density of residuals area unconditioned
            y.boot <- NULL
            y.boot.i <- NULL
            for (i in 1:myareas) {
              y.boot.i <-
                myX[myregioncodepop == myar[i], ] %*% mycoef[, i] + sample(res.s.centered, mypopsize[i], replace =
                                                                             TRUE)
              y.boot <- c(y.boot, y.boot.i)
            }
          }


          for (ii in 1:myareas) {
            myq.true.boot[kk, ii, b] <-
              quantile(y.boot[myregioncodepop == myar[ii]], prob = c(qest))
          }

          myarea.q.boot.r[kk, , b, 1:R] <-
            boot.CD.R(
              qest,
              y.boot,
              myX,
              myregioncodepop,
              mypopsize,
              mysamplesize,
              myar,
              myareas,
              adjseed,
              R,
              id
            )
          #print(myarea.q.boot.r[kk,,b,])
          #myarea.q.boot.r[kk,,b,(R/2+1):R]<-mycollect[[2]]
          #print(myarea.q.boot.r[kk,,b,])

          for (i in 1:myareas) {
            myarea.q.boot[kk, i, b] <- mean(myarea.q.boot.r[kk, i, b, ])
          }

        }#B ends here

        for (i in 1:myareas) {
          myq.true.boot.m[kk, i] <- mean(myq.true.boot[kk, i, ])
          myarea.q.boot.m[kk, i] <- mean(myarea.q.boot[kk, i, ])
        }

        for (i in 1:myareas) {
          BIAS.boot[kk, i] <- myarea.q.boot.m[kk, i] - myq.true.boot.m[kk, i]
          aux <- matrix(0, B, 1)
          for (b in 1:B) {
            aux[b, 1] <-
              (1 / R) * sum((myarea.q.boot.r[kk, i, b, ] - myarea.q.boot[kk, i, b]) ^
                              2)
          }
          VAR.boot[kk, i] <- (1 / B) * sum(aux[, 1])
        }

        for (i in 1:myareas) {
          MSE.boot[kk, i] <- ((BIAS.boot[kk, i]) ^ 2) + VAR.boot[kk, i]
        }

      }#end if myMSE
    }#Iteration k ends here
    RMSE.CD <- as.data.frame(sqrt(MSE.boot))
    names(RMSE.CD) <- myar
    row.names(RMSE.CD) <- myqgrid
    quantiles <- as.data.frame(myarea.q)
    names(quantiles) <- myar
    row.names(quantiles) <- myqgrid
    rest <-
      list(quantiles = quantiles,
           rmse = RMSE.CD,
           Area.Code = myar)#,var=VAR.boot,bias=BIAS.boot,qtrueboot=myq.true.boot
  }
