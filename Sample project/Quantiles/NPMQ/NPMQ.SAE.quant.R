rm(list = ls(all.names = TRUE))


set.seed(1977)

library(MASS)
library(np)
library(SemiPar)
library(splines)


npqrlm <-
  function(Y.n,
           X,
           Z.spline,
           quantile,
           tol = 0.001,
           maxit = 100,
           theta = 2,
           kk = 1.345) {
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
    Z.spline <- as.matrix(Z.spline)
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



mq.coef <-
  function(myx,
           myy,
           myzspline,
           mykvalue = 1.345,
           myregioncode,
           maxiter = 100) {
    #This function estimate the m-quantile regression coefficients
    #myx<- x sample matrix of auxiliary variables
    #myy<- y vector
    #mynumauxvar<- number of auxiliary variables (include constant)
    #myregioncode<- area code for y and x units, data must be ordered by area code
    #maxiter<- OPTIONAL, number of maximum iteration for ob algorithm

    myar <- sort(unique(myregioncode))
    myzspline <- as.matrix(myzspline)
    myareas <- length(myar)
    mysamplesize <- sum(as.numeric(table(myregioncode)))
    mynumauxvar <- dim(myx)[2] + ncol(myzspline)

    ob <-
      npqrlm(myy, myx, myzspline, q = sort(c(
        seq(0.006, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98
      )), kk = mykvalue)
    q.values <-
      sort(c(seq(0.006, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98))
    qo <-
      matrix(c(gridfitinter(y, ob$hat.values, q.values)), nrow = mysamplesize, ncol =
               1)
    qmat <- matrix(c(qo, myregioncode),
                   nrow = length(myregioncode),
                   ncol = 2)
    mqo <- aggregate(qmat[, 1], list(d2 = qmat[, 2]), mean)[, 2]
    saq <- matrix(c(mqo, myar), nrow = myareas, ncol = 2)
    saq <- rbind(saq, c(0.5, 9999))
    ob1 <-
      npqrlm(myy,
             myx,
             myzspline,
             quantile = c(mqo[1:myareas]),
             kk = mykvalue)
    mycoef <-
      matrix(c(t(ob1$b.stim)), nrow = myareas, ncol = mynumauxvar) # need to be ordered by area
    mycoef <- t(mycoef)
    mycoef
  }


intsolver <-
  function(myqest,
           myyboot,
           myX,
           myzspline,
           myzsplinepop,
           myregioncodepop,
           mypopsize,
           myar,
           myareas,
           adjseed,
           mysboot) {
    myres <- array(0, dim = c(myareas))

    myy <- myyboot[mysboot]
    myx <- myX[mysboot, ]
    myregioncode <- myregioncodepop[mysboot]
    myregioncoder <- myregioncodepop[-mysboot]
    mysamplesizer <- as.numeric(table(myregioncoder))
    myX.r <- myX[-mysboot, ]
    myzspline <- myzsplinepop[mysboot, ]
    myzspliner <- myzsplinepop[-mysboot, ]

    # M-quantiles
    coef.boot <- mq.coef(myx, myy, myregioncode)

    # Quantile Estimation Using Chambers Dunstan Estimator
    for (i in 1:myareas) {
      f1 <- myy[myregioncode == myar[i]]
      pred.medr <-
        cbind(myX.r[myregioncoder == myar[i], ], myzspliner[myregioncoder == myar[i], ]) %*%
        coef.boot[, i]
      pred.meds <-
        cbind(myx[myregioncode == myar[i], ], myzspline[myregioncode == myar[i], ]) %*%
        coef.boot[, i]
      res.s <- f1 - pred.meds
      #z<-matrix(rep(res.s,sample.sizer[i]),nrow=sample.sizer[i],ncol=sample.size[i])
      z <- sample(res.s, mysamplesizer[i], replace = TRUE)
      z <- z + pred.medr
      comb <- c(f1, pred.medr)
      start0 <- quantile(comb, prob = c(myqest))
      sameside <- T

      while (sameside) {
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
      }

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

        while (abs(tol) >= 0.5) {
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
        }
      }

      if (fdif < 0.01) {
        start.med <- (start.bef + start.aft) / 2
        eps <- 0.001
        tol <- 50

        while (abs(tol) > 0.5) {
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
        }
      }

      myres[i] <- start.med
      #print(i)
    }#Iteration i in bootstrap ends here
    myres
  }



NPMQ.SAE.quant <-
  function(myqgrid,
           myy,
           myx,
           myX,
           myzspline,
           myzsplinepop,
           myregioncode,
           myregioncodepop,
           adjseed = max(0.15, mean(myy) / 500)) {
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
    #method<- which method to be used to estimate residuals distribution, choice= "su" (smooth unconditional),"eu" (emprirical unconditional),"sc" (smooth conditional),"ec" (empirical uncoditional)


    myar <- unique(myregioncode)
    myareas <- length(myar)
    mypopsize <- as.numeric(table(myregioncodepop))
    mysamplesize <- as.numeric(table(myregioncode))
    myquantnum <- length(myqgrid)
    id <- seq(1:sum(mypopsize))

    myarea.q <- matrix(0, myquantnum, myareas)
    myq.true.boot.m <- array(0, dim = c(myquantnum, myareas))
    myarea.q.boot.m <- array(0, dim = c(myquantnum, myareas))
    kk <- 0

    mycoef <- mq.coef(myx, myy, myzspline, mykvalue = 1.345, myregioncode)

    for (qq in myqgrid) {
      qest <- qq
      kk <- kk + 1
      myres <- NULL

      for (i in 1:myareas) {
        f1 <- myy[myregioncode == myar[i]]
        pred.medtot <-
          cbind(myX[myregioncodepop == myar[i], ], myzsplinepop[myregioncodepop ==
                                                                  myar[i], ]) %*% mycoef[, i]
        x.design.i <- as.matrix(myx[myregioncode == myar[i], ])
        pred.meds <-
          cbind(myx[myregioncode == myar[i], ], myzspline[myregioncode == myar[i], ]) %*%
          mycoef[, i]
        res.s <- f1 - pred.meds
        myres[i] <- list(res.s)
        #			z<-matrix(rep(res.s,pop.size[i]),nrow=pop.size[i],ncol=sample.size[i])
        z <- sample(res.s, mypopsize[i], replace = TRUE)
        z <- z + pred.medtot
        comb <- c(pred.medtot)
        start0 <- quantile(comb, prob = c(qest))
        sameside <- T

        while (sameside) {
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
        }

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

          while (abs(tol) >= 0.5) {
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
          }
        }

        if (fdif < 0.01) {
          start.med <- (start.bef + start.aft) / 2
          eps <- 0.001
          tol <- 50

          while (abs(tol) > 0.5) {
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
          }
        }
        #print(i)
        myarea.q[kk, i] <- start.med
      }#Iteration i ends here

    }
    quantiles <- myarea.q
    for (i in 1:myareas) {
      check <- sort(quantiles[, i]) - quantiles[, i]
      check <- sum(abs(check))
      if (check != 0) {
        warning("Quantile crossing produced in area ", i)
      }
    }
    rest <- list(quantiles = myarea.q, Area.Code = myar)

  }
