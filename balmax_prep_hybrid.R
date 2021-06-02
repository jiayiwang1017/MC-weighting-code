library(foreach)
library(doParallel)
library(abind)
library(Matrix)

############ balmaxfit cpp function ########
#### mu_seq is not scaled, and is not multiplicated by R0 ####
balmaxfit_cv = function(Atrue, Aobs, omega, what, alpha0, nfolds = 5, R0_seq, mu_seq, rho = 0.1,
  thresh_admm = 1e-04, maxit_admm = 200, traceit = F, CV_types = c("uni", "weighted",
    "true"), Aini = NULL) {


  n1 = dim(Aobs)[1]
  n2 = dim(Aobs)[2]
  if (is.null(Aini)) {
    Wini = rbind(cbind(diag(1, n1), Aobs), cbind(t(Aobs), diag(1, n2)))

    Xini = 1.0 * Wini

    Zini = diag(1, n1 + n2)

  } else {
    Wini = rbind(cbind(diag(1, n1), Aini), cbind(t(Aini), diag(1, n2)))

    Xini = 1.0 * Wini

    Zini = diag(1, n1 + n2)
  }


  iobs = which(omega == 1)
  folds <- cut(sample(seq(1, length(iobs))), breaks = nfolds, labels = FALSE)
  validateIndexes = list()
  for (i in 1:nfolds) {
    temp_idx <- which(folds == i, arr.ind = TRUE)
    validateIndexes[[i]] <- iobs[temp_idx]
  }

  #registerDoParallel(nfolds)
  #acomb <- function(...) abind(..., along = 4)
  ### just do validation, no cross validation 
  #cvbalmax = foreach(fold_num = 1:nfolds, .combine = "acomb", .multicombine = TRUE) %dopar% {
  fold_num = 1
  W = 1 * Wini
  X = 1 * Xini
  Z = 1 * Zini

  omegavalidate = rep(0, n1 * n2)
  # iomegavalidate1[ivalidate1] = 1 omegavalidate1 = matrix(iomegavalidate1, n1,
  # n2)
  omegavalidate[validateIndexes[[fold_num]]] = 1
  omegatrain = as.vector(omega) - omegavalidate

  Aobstrain_vec = Aobs[omegatrain == 1]
  Aobsvalidate_vec = Aobs[omegavalidate == 1]
  whattrain = what[omegatrain == 1]
  whatvalidate = what[omegavalidate == 1]
  Temp = Aobs
  Temp[omegavalidate == 1] = 0
  sc = svd(Temp)$d[1]
  rm(Temp)
  cvbalmax_values = balmaxE_cv_cpp(Atrue, Aobstrain_vec, Aobsvalidate_vec, omegatrain,
        omegavalidate, n1, n2, whattrain, whatvalidate, R0_seq, mu_seq * sc, alpha0, thresh_admm,
        maxit_admm, traceit, rho, W, X, Z)



  #cvbalmax_values = apply(cvbalmax, c(1, 2, 3), sum)
  #colnames(cvbalmax_values) = CV_types
  cvbalmax_grid = NULL
  select_tuning = list()

  for (k in 1:length(CV_types)) {
    minerror = min(cvbalmax_values[,, k][cvbalmax_values[,, k] > 0])
    idx = which(cvbalmax_values[,, k] == minerror)[1]
    cvbalmax_grid = c(cvbalmax_grid, idx)
    tun_par = list()
    if (idx %% length(R0_seq) == 0) { tun_par$R0 = R0_seq[length(R0_seq)] } else {
      tun_par$R0 = R0_seq[idx %% length(R0_seq)]
    }
    tun_par$mu = mu_seq[idx %/% length(R0_seq) + 1] * tun_par$R0
    select_tuning[[k]] = tun_par

  }



  return(list(cvbalmax_values = cvbalmax_values, cvbalmax_grid = cvbalmax_grid, select_tuning = select_tuning))

}

############################################################# replace missing values by column means ##########

A_initial <- function(A, omega) {
  A[omega == 0] = NA
  A = apply(A, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })
  return(A)
}


balmaxfit_true <- function(Atrue, Aobs, omega, what, alpha0, R0_seq, mu_seq, rho = 0.1, thresh_admm = 1e-04,
  maxit_admm = 200, traceit = F) {

  n1 = dim(Aobs)[1]
  n2 = dim(Aobs)[2]


  Wini = rbind(cbind(diag(1, n1), Aobs), cbind(t(Aobs), diag(1, n2)))

  Xini = 1.0 * Wini

  Zini = diag(1, n1 + n2)

  mu_seq = mu_seq * svd(Aobs)$d[1]

  Aobs_vec = Aobs[omega == 1]
  what_vec = what[omega == 1]
  OUT = balmaxE_true_cpp(Atrue, Aobs_vec, as.vector(omega), n1, n2, what_vec, R0_seq, mu_seq, alpha0, thresh_admm, maxit_admm, traceit, rho, Wini, Xini, Zini)

  return(OUT)
}
