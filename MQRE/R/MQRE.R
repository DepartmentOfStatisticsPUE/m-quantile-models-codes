# require(MASS)
# require(nlme)
# require(Matrix)
# require(Rcpp)
# require(RcppArmadillo)
# require(RcppEigen)
#' The function provides the estimates based on M-quantile random effects regression.
#'
#' @examples
#' # Data generation
#' library(saeSim)
#' set.seed(1)
#' groups<-30 # Define the number of groups
#' sample_size<-sample(x=5:15,size=groups,replace=TRUE) # Define the sample sizes
#' setup<-sim_base(data = base_id(nDomains=groups,nUnits=sample_size)) %>%
#'  sim_gen_x(mean=0,sd=10)       %>%
#'  sim_gen_v(mean=0,sd=300)         %>%
#'  sim_gen_e(mean=0,sd=800)         %>%
#'  sim_resp_eq(y = 4500 - 400*x + v + e)
#'  # Generate 10 samples
#' Pop<-sim(setup,R=10)
#' sample_data<-Pop[[1]] # Select the first sample for the test
#'
#' # Fit the MQRE model
#' fit<-MQRE_function(formula=y~x,saind=sample_data$idD,data=sample_data,qtl=0.4)
#' @export
#' @useDynLib MQRE
#' @param formula An object of class "formula"
#' @param saind Small area indicator (specification of the grouping)
#' @param qtl Define the M-quantile
#' @param data Data frame containing the variables in the model
#' @param tol Tolerance
#' @param maxit Maximum number of iterations
#' @param k Tuning constant
#' @param k_v Tuning constant for the estimation of the random effects
#' @return Returns a list including:
#' \itemize{
#' \item coefficients - Regression parameters
#' \item stde.beta - Standard error of the regression parameters
#' \item sigma2u - Varianz of the random effects (Level 2)
#' \item std.u - Standard error of the varianz of the random effects (Level 2)
#' \item sigma2e - Varianz of the error term (Level 1)
#' \item std.e - Standard error of the varianz of the error term (Level 1)
#' \item quantile - Estimated M-quantile
#' \item samplesize - Sample size in each group
#' \item iterations - Number of iterations in the fixed point algorithm
#' \item Kvalue - Value of the K matrix in the estimation equations
#' \item rand.eff - Vector of the random effects
#' }
#' @references
#' Tzavidis, N.; Salvati, N.; Schmid, T.; Flouri, E. and Midouhas, E.: Longitudinal analysis of the Strengths and Difficulties Questionnaire scores of the Millennium Cohort Study children in England using M-quantile random effects regression, Journal of the Royal Statistical Society: Series A, accepted.
#'

MQRE_function <-
  function(formula,
           saind,
           qtl,
           data,
           tol = 1e-04,
           maxit = 100,
           k = 1.345,
           k_v = 1.345)

  {
    #sourceCpp("C++\\ersterVersuch.cpp")

    makeXY <- function(formula, data) {
      mf <- model.frame(formula = formula, data = data)
      x <- model.matrix(attr(mf, "terms"), data = mf)
      y <- model.response(mf)

      list(y = y,
           x = x)
    }
    daten <- makeXY(formula = formula, data = data)

    y <- daten$y
    group <- saind
    x <- daten$x
    n = length(y)
    ni = table(group)
    m = length(ni)
    areanumber = m
    p = ncol(x)
    gg = unique(group)


    z = matrix(0, n, m)
    kk = 0
    for (j in 1:m) {
      for (i in 1:ni[j]) {
        kk = kk + 1
        z[kk, j] = 1
      }
    }

    const <-
      4 * k ^ 2 * (1 - pnorm(k)) * ((1 - qtl) ^ 2 + qtl ^ 2) - 4 * k * dnorm(k) *
      ((1 - qtl) ^ 2 + qtl ^ 2) + 4 * (1 - qtl) ^ 2 * (pnorm(0) - (1 - pnorm(k))) + 4 *
      qtl ^ 2 * (pnorm(k) - pnorm(0))
    K2_value <- const
    K2 <- diag(K2_value, n)


    #STEP 1 (Definition of the starting values)
    w = cbind(x[, -1], z)
    fit.H = lm(y ~ w)
    e.1.H = residuals(fit.H)
    sigma.e.hat.H = sum(e.1.H ^ 2) / (n - (qr(w)$rank))

    xx = x
    fit2.H = lm(y ~ xx[, -1])
    e.2.H = residuals(fit2.H)
    A.H = sum(e.2.H ^ 2)
    xx = as.matrix(xx)
    B.H = sigma.e.hat.H * (n - qr(xx)$rank)
    o.H = diag(n) - (xx %*% (ginv(t(xx) %*% xx) %*% t(xx)))
    C.H = sum(diag(t(z) %*% o.H %*% z))
    sigma.v.hat.H = (A.H - B.H) / C.H
    sigmasq0 = sigma.e.hat.H #initial values of sigma_sq_e#
    sigmasq0.v = sigma.v.hat.H #initial values sigma_sq_V#

    if (dim(x)[2] == 1) {
      ls <- lm(y ~ 1)
    } else{
      ls <- lm(y ~ x[, -1])
    }
    beta1 <- c(ls$coefficients)
    estsigma2u <- sigmasq0.v
    estsigma2e <- sigmasq0
    sigma1 <- abs(c(sigmasq0, sigmasq0.v))

    # Estimation of beta and sigma

    ZZ <- z %*% t(z)

    E.psi.der      <- pnorm(k) - pnorm(-k)

    estimators1 <-
      EstFunc_ML2(
        y = y,
        x = x,
        Z = z,
        sigma_v = sigma1[2],
        sigma_e = sigma1[1],
        m = m,
        beta = beta1,
        tol = tol,
        maxit = maxit,
        k = k,
        K2_val = K2_value,
        qtl = qtl,
        k_val = E.psi.der
      )




    beta1 <- estimators1$beta
    sigma1 <- estimators1$sigma

    Vmatrix <- makeVMat(
      Z = z,
      sigma_v = sigma1[2],
      sigma_e = sigma1[1],
      m = m
    )
    V <- Vmatrix$V
    V.inv <- Vmatrix$Vinv
    U <- Diagonal(x = diag(V))
    U.inv <- chol2inv(chol(U))
    xbeta <- c(x %*% beta1)
    r <- c(sqrt(U.inv) %*% (y - xbeta))
    psi.r <- psi_Q(u = as.vector(r), q = qtl, k = k)
    der.psi.r <- psi_Q_der(u = as.vector(r), q = qtl,  k = k)



    si = array(c(rep(0, (m * 2))), dim = c(2, 1, m))
    Hi = array(c(rep(0, (m * 2))), dim = c(2, 2, m))

    for (i in 1:m)
    {
      zi = matrix(rep(1, ni[i]), ni[i], 1)
      ZZi <- zi %*% t(zi)
      s <- matrix(0, 2, 1)
      Hessian = matrix(0, 2, 2)

      Vh <- sigma1[1] * diag(1, ni[i], ni[i]) + sigma1[2] * zi %*% t(zi)
      Vih <- chol2inv(chol(Vh))
      Uh <- diag(diag(Vh), ni[i], ni[i])
      Uih <- chol2inv(chol(Uh))

      res1 <- r[group == i]
      res2 <- psi.r[group == i]
      res3 <- der.psi.r[group == i]

      K2.i <- diag(K2_value, ni[i])
      Part1 <- Vih %*% sqrt(Uh) %*% matrix(c(res2), ni[i], 1)
      Part2 <- Vih %*% ZZi

      s[1, 1] <-
        as.numeric(0.5 * t(Part1) %*% ZZi %*% Part1 - 0.5 * sum(diag(K2.i %*% Part2)))
      s[2, 1] <-
        as.numeric(0.5 * t(Part1) %*% Part1 - 0.5 * sum(diag(K2.i %*% Vih)))
      si[, 1, i]  <- s

      Part3 <- t(diag(res3, ni[i]) %*% Uih %*% res1) %*% sqrt(Uh) %*% Vih
      Part4 <- t(res2) %*% sqrt(Uih) %*% Vih
      Part5 <- t(res2) %*% sqrt(Uh) %*% Vih
      Hessian[1, 1] = ((-0.5 * Part3 + 0.5 * Part4 - Part5 %*% ZZi %*% Vih) %*%
                         ZZi %*% t(Part5) + 0.5 * sum(diag(K2.i %*% Part2 %*% Part2)))[1, 1]

      Hessian[1, 2] = ((-0.5 * Part3 + 0.5 * Part4 - Part5 %*% Vih) %*%
                         ZZi %*% t(Part5) + 0.5 * sum(diag(K2.i %*% Vih %*% Part2)))[1, 1]

      Hessian[2, 1] = ((-0.5 * Part3 + 0.5 * Part4 - Part5 %*% ZZi %*% Vih) %*%
                         t(Part5) + 0.5 * sum(diag(K2.i %*% Part2 %*% Vih)))[1, 1]

      Hessian[2, 2] = ((-0.5 * Part3 + 0.5 * Part4 - Part5 %*% Vih) %*%
                         t(Part5) + 0.5 * sum(diag(K2.i %*% Vih %*% Vih)))[1, 1]
      Hi[, , i] <- Hessian

    }

    B2 = matrix(0, 2, 2)

    B2[1, 1] = sum(si[1, 1, ] * si[1, 1, ])
    B2[1, 2] = sum(si[1, 1, ] * si[2, 1, ])
    B2[2, 1] = sum(si[2, 1, ] * si[1, 1, ])
    B2[2, 2] = sum(si[2, 1, ] * si[2, 1, ])

    C2 = matrix(0, 2, 2)
    C2[1, 1] = sum(Hi[1, 1, ])
    C2[1, 2] = sum(Hi[1, 2, ])
    C2[2, 1] = sum(Hi[2, 1, ])
    C2[2, 2] = sum(Hi[2, 2, ])

    I.matrix <- solve(C2) %*% B2 %*% solve(C2)
    stde.u = sqrt(I.matrix[1, 1])
    stde.e = sqrt(I.matrix[2, 2])

    #STEP 3

    si = array(c(rep(0, (m * p))), dim = c(p, p, m))
    Hi = array(c(rep(0, (m * p))), dim = c(p, p, m))

    for (i in 1:m)
    {
      zi = matrix(rep(1, ni[i]), ni[i], 1)
      Vh <- sigma1[1] * diag(1, ni[i], ni[i]) + sigma1[2] * zi %*% t(zi)
      Vih <- chol2inv(chol(Vh))
      Uh <- diag(diag(Vh), ni[i])
      Uih <- chol2inv(chol(Uh))

      xi <- matrix(x[group == i, ], nrow = ni[i], ncol = p)

      res1 <- r[group == i]
      res2 <- psi.r[group == i]
      res3 <- der.psi.r[group == i]

      si[, , i] <-
        as.matrix(t(xi) %*% Vih %*% sqrt(Uh) %*% res2 %*% t(res2) %*% sqrt(Uh) %*%
                    Vih %*% xi)
      Hi[, , i] <- as.matrix(t(xi) %*% Vih %*% diag(res3, ncol = ni[i]) %*% xi)
    }

    B1 = matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        B1[i, j] = sum(si[i, j, ])
      }
    }

    C1 = matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        C1[i, j] = sum(Hi[i, j, ])
      }
    }

    var.beta = solve(C1) %*% B1 %*% t(solve(C1))

    sqrt.beta = NULL
    for (oo in 1:p)
    {
      sqrt.beta[oo] = sqrt(var.beta[oo, oo])
    }

    zd <- beta1 / (sqrt.beta)
    sig <- 0

    for (oo in 1:p) {
      {
        if (zd[oo] > 0)
          sig[oo] <- (1 - pnorm(zd[oo])) * 2
        else
          sig[oo] <- (pnorm(zd[oo])) * 2
      }
    }


    Beta = beta1
    Std.Error = (sqrt.beta)
    z_value = (zd)
    p_value = (sig)
    risultati <- cbind(Beta[, 1], Std.Error, z_value[, 1], p_value)

    nnn = NULL
    if (p > 1) {
      for (oo in 2:p)
      {
        tmp.n = paste("X", (oo - 1))
        nnn = c(nnn, tmp.n)
      }
    }

    if (p == 1)
      res = data.frame(round(risultati, 5), row.names = c("Const"))
    if (p > 1)
      res = data.frame(round(risultati, 5), row.names = c("Const", nnn))
    #print(res)

    Value = (c(sigma1[1], sigma1[2]))
    Std.Error = (c(stde.e, stde.u))
    I.low.e = sigma1[1] * exp(qnorm(0.025, 0, 1) * (1 / sigma1[1]) * stde.e)
    I.low.u = sigma1[2] * exp(qnorm(0.025, 0, 1) * (1 / sigma1[2]) * stde.u)
    I.up.e = sigma1[1] * exp(qnorm(0.975, 0, 1) * (1 / sigma1[1]) * stde.e)
    I.up.u = sigma1[2] * exp(qnorm(0.975, 0, 1) * (1 / sigma1[2]) * stde.u)
    low <- c(I.low.e, I.low.u)
    up <- c(I.up.e, I.up.u)
    risultati <- cbind(Value, Std.Error, low, up)

    res = data.frame(round(risultati, 5),
                     row.names = c("Variance level 1", "Variance level 2"))
    #print(res)

    # Prediction Random Effects
    E.psi.der      <- pnorm(k_v) - pnorm(-k_v)


    rob.v       <-
      RandEff(
        y = y,
        x = x,
        Z = z,
        sigma_v = sigma1[2],
        sigma_e = sigma1[1],
        m = m,
        beta = beta1[, 1],
        tol = tol,
        maxit = maxit,
        k = k,
        k_v = k_v,
        qtl = qtl,
        k_val = E.psi.der
      )

    # List of Results
    list(
      coefficients = beta1,
      stde.beta = sqrt.beta,
      sigma2u = (sigma1[2]),
      std.u = stde.u,
      sigma2e = (sigma1[1]),
      std.e = stde.e,
      quantile = qtl,
      iterations = estimators1$Iteration,
      Kvalue = (K2[1, 1]),
      rand.eff = rob.v
    )

  }
