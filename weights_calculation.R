library("irlba")
#########################################

eval_obj_grad <- function(w, d1, d2, ind, kap) {

  z <- matrix(-1, d1, d2)
  z[ind == 1] <- w - 1
  B = z

  Bsing = irlba(B, 1 , maxit = 3000)

  maxsing <- Bsing$d[1]

  u1 = Bsing$u[, 1]
  v1 = Bsing$v[, 1]

  wmat = matrix(0, d1, d2)
  wmat[ind == 1] = w

  full_grad = ind * u1 %*% t(v1) + 2 * kap * wmat

  mygrad = full_grad[ind == 1]

  V <- sum(w^2)

  return(list(objective = maxsing + kap * V, gradient = as.vector(mygrad)))  # modified
}


eval_obj <- function(w, d1, d2, ind, kap) {

  z <- matrix(-1, d1, d2)
  z[ind == 1] <- w - 1
  B = z

  Bsing = irlba(B, 1, work = 20, maxit = 2000)

  maxsing <- Bsing$d[1]

  V <- sum(w^2)
  return(maxsing + kap * V)
}


eval_grad <- function(w, d1, d2, ind, kap) {

  z <- matrix(-1, d1, d2)
  z[ind == 1] <- w - 1
  B = z

  Bsing = irlba(B, 1, work = 20, maxit = 2000)


  u1 = Bsing$u[, 1]
  v1 = Bsing$v[, 1]

  wmat = matrix(0, nrow = d1, ncol = d2)
  wmat[ind == 1] = w

  #full_grad = ind %*% v1 %*% t(u1) + 2 * kap * wmat
  full_grad = ind * u1 %*% t(v1) + 2 * kap * wmat

  mygrad = full_grad[ind == 1]

  return(as.vector(mygrad))  # modified
}

#########################################

#' MCB core function

MCB.core <- function(ind, kaps, lower = 1, upper = Inf, traceit = TRUE, w0 = NULL,
  maxit = 2000, xtol_rel = 1e-08, full = FALSE, algorithm = "lbfgsb3c") {

  N <- length(ind)
  d1 <- dim(ind)[1]
  d2 <- dim(ind)[2]

  # construct a vector
  tind <- as.logical(ind)
  n1 <- sum(tind)
  n2 <- sum(!tind)

  # lower and upper bound
  if (length(lower) == 1)
    lower <- rep(lower, n1)
  if (length(upper) == 1)
    upper <- rep(upper, n1)

  if (is.null(w0)) {
    w0 <- lower + 1
  }

  nkaps <- length(kaps)
  outlist <- list()
  ws <- matrix(1, nr = nkaps, nc = n1)
  SNs <- array(dim = nkaps)

  rmseTo1s <- array(dim = nkaps)
  outlist <- NULL
  obj1 <- array(dim = nkaps)

  for (i in (1:nkaps)) {

    if (traceit)
      cat("####", i, ":", "\tkappa =", kaps[i], "\n")
    if (algorithm == "lbfgsb3c") {
      res <- lbfgsb3c::lbfgsb3c(par = w0, fn = eval_obj, gr = eval_grad, lower = lower,
        upper = upper, d1 = d1, d2 = d2, ind = ind, kap = kaps[i], control = list(factr = xtol_rel,
          maxit = maxit))
      res$par[res$par < 1] = 1

      obj1[i] <- res$value
      ws[i, ] <- res$par

      cat(paste(res$message, ' '))
      cat(res$counts)
      cat('\n')

    } else {
      res <- nloptr(x0 = w0, eval_f = eval_obj_grad, lb = lower, ub = upper,
        d1 = d1, d2 = d2, ind = ind, kap = kaps[i], opts = list(algorithm = "NLOPT_LD_LBFGS",
          xtol_rel = xtol_rel, print_level = 0, maxeval = maxit, check_derivatives = F))
      res$par[res$par < 1] = 1
      obj1[i] <- res$objective
      ws[i, ] <- res$sol
      cat(res$message)
      cat(res$counts)
      cat('\n')
    }
  ############ trying ###################
  #res <- lbfgsb3c::lbfgsb3c(par = w0, fn = eval_obj, lower = lower,
  #  upper = upper, d1 = d1, d2 = d2, ind = ind, kap = kaps[i], control = list(reltol = xtol_rel,
  #    maxit = maxit))
  #######################################


    if (full)
      outlist[[i]] <- res



    temp = matrix(-1, d1, d2)
    temp[ind == 1] = ws[i, ] - 1

    # ee <- eigen(temp %*% t(temp))
    ee <- irlba(temp, 1, maxit = 3000)
    # compute SN
    SNs[i] <- ee$d[1]

    # unorms[i] <- sqrt(sum(nq1/(-kaps[i]) * ee$vectors[, 1]^2))

    # fittedus[i, ] <- as.vector(nP1 %*% ee$vectors[, 1] * (N))

    rmseTo1s[i] <- sqrt(mean((ws[i, ] - 1)^2))

    w0 <- ws[i, ]
  }
  temp = matrix(-1, d1, d2)
  temp[ind == 1] = 0
  SNmax = irlba(temp, 1, maxit = 3000)$d[1]

  return(list(outlist = outlist, ws = ws, SNs = SNs, rmseTo1s = rmseTo1s, kaps = kaps,
    obj1 = obj1, SNmax = SNmax))
}


#########################################

#' MCB with kappa selection


MCB.SN <- function(ind, kaps, lower = 1, upper = Inf, traceit = TRUE, w0 = NULL,
  maxit = 2000, xtol_rel = 1e-08, method = 1, full = FALSE, algorithm = "lbfgsb3c") {
kaps <- sort(kaps)
  ores <- MCB.core(ind = ind, kaps = kaps, lower = lower, upper = upper, traceit = traceit,
    w0 = w0, maxit = maxit, xtol_rel = xtol_rel, full = full)

  # find which SN is the smallest
  if(length(kaps) == 1){kap_ind = 1}else {
     if (method == 1) {

    kap_ind <- which.min(ores$SNs)[1]

  } else if (method == 2) {
diffsignrev = rev(diff(ores$SNs)<0)*1
    diffsignrevcum = cumsum(diffsignrev)
    idx = which(diffsignrevcum>0)[1]
    idx2 = length(kaps)-idx + 1
    kap_ind <- which((diff(ores$SNs)/diff(ores$kaps) < (1e3)))
    kap_ind <- which(kap_ind %in% c(idx2:length(kaps)))[1]

  }

  if(is.na(kap_ind)){
    kap_ind <- which.min(ores$SNs)[1]
  }
  }
  

  return(list(w = ores$ws[kap_ind, ], kap_ind = kap_ind, warns = c(kap_ind == 1,
    kap_ind == length(kaps)), ores = ores, kaps = kaps))
}
