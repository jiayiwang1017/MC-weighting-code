library(foreach)
library(doParallel)
library(abind)
library(Matrix)


shrink_rownorm <- function(A, t){
  rownorms <- apply(A, 1, function(x) sqrt(sum(x^2)))
  idx <- rownorms> t
  A[idx, ] <-  A[idx,]/rownorms[idx]*t
  return(A)
}


A_initial <- function(A, omega) {
  A[omega == 0] = NA
  A = apply(A, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })
  return(A)
}



##### initialize with convex solution ####



maxnorm_nonconvex_hybrid <- function(Aobs, omega, what, Aini = NULL, L0 = NULL, R0 = NULL, beta, mu,tau = 1,
  r = NULL, maxiter = 3000, retol = 1e-06) {
  ### beta: maxnorm tau: step-size r: rank (can restrict to small number when
  ### dimension is very large)

  n1 <- nrow(Aobs)
  n2 <- ncol(Aobs)
  if (is.null(r)) {
    r <- min(n1, n2)
  }
  if (is.null(L0) & is.null(R0)) {
    if(is.null(Aini)) {Aini = A_initial(Aobs, omega)}
    SVDA <- svd(Aini)
    L0 <- SVDA$u[, 1:r] %*% diag(sqrt(SVDA$d[1:r]))
    R0 <- SVDA$v[, 1:r] %*% diag(sqrt(SVDA$d[1:r]))
  }
  retol_L <- rep(0, maxiter)
  retol_R <- rep(0, maxiter)
  obj <- rep(0, maxiter)
  obj[1] = 1/2 * sum(((L0 %*% t(R0) - Aobs) * omega)^2 * what)+ 0.5 * mu * (norm(L0,"F")^2 + norm(R0,"F")^2)

  for (i in 1:maxiter) {
    grad_L <-((L0 %*% t(R0) - Aobs) * omega * what) %*% R0 + mu * L0
    L <- L0 - tau * grad_L
    L <- shrink_rownorm(L, sqrt(beta))
    grad_R <- ((R0 %*% t(L) - t(Aobs)) * t(omega) * t(what)) %*% L + mu * R0

    R <- R0 - tau * grad_R
    R <- shrink_rownorm(R, sqrt(beta))

    #retol_L[i] <- norm(L-L0,"F")/norm(L0,"F")
    #retol_R[i] <- norm(R-R0,"F")/norm(R0,"F")
    obj[i+1] <- 1/2 * sum(((L %*% t(R) - Aobs) * omega)^2 * what) + 0.5 * mu * (norm(L,"F")^2 + norm(R,"F")^2)
    if((i>10) && (obj[i] - obj[i+1])/obj[i] < retol && (obj[i+1] <= obj[i]) #&& max(retol_L,retol_R)<retol
    ){
      break
    }
    if(( obj[i+1] > obj[i])){
    tau = 0.5 * tau
    #print(tau)
    L0 = L0
    R0 = R0
    obj[i+1] = obj[i]
    }
    else{
      L0 = L
      R0 = R
    }
    

  }

  return(list(L = L, R = R, Ahat = L %*% t(R), iterations = i, errors = cbind(retol_L[1:i],retol_R[1:i]), obj = obj[1:i],tau = tau))
}





### mu_seq not scaled and not multiply by beta  ###
nonconvex_hybrid_validation <- function(Aobs, omega, what, Aobsvalidate, omegavalidate, Aobstest, omegatest, Aini = NULL, L0 = NULL, R0 = NULL, beta_seq, mu_seq, tau = 10,
  r = NULL, maxiter = 2500, retol = 1e-06){
      n1 = dim(Aobs)[1]
      n2 = dim(Aobs)[2]

      if (is.null(r)) {
        r <- min(n1, n2)
      }

      
      if (is.null(L0) & is.null(R0)) {
        if(is.null(Aini)) {Aini = A_initial(Aobs, omega)}
        SVDA <- svd(Aini)
        L0 <- SVDA$u[, 1:r] %*% diag(sqrt(SVDA$d[1:r]))
        R0 <- SVDA$v[, 1:r] %*% diag(sqrt(SVDA$d[1:r]))
      }

            sc = svd(Aobs)$d[1]
            #OUTtemp = NULL
            cv_errors = matrix(0, nrow = length(beta_seq), ncol = length(mu_seq))
            test_errors = matrix(0, nrow = length(beta_seq), ncol = length(mu_seq))
            Ranks = matrix(0, nrow = length(beta_seq), ncol = length(mu_seq))
            #taus = rep(0,length(R0_seq))
            for(i in 1:length(beta_seq)){
              for(j in 1:length(mu_seq)){
                cat(c(i,j))
                cat("###########")
                  if((i==1) && (j==1)){

                    #W = rbind(cbind(diag(1, n1), Aobs), cbind(t(Aobs), diag(1, n2)))
                    #X = 1 * W
                    #Z = diag(1, n1 + n2)
                    #OUT = balmaxE_cpp_wraper(Aobs = Aobs, n1, n2, omega, what, beta_seq[i], mu_seq[j]*sc * beta_seq[i] ,alpha0, rho = 0.1, W, X, Z, thresh_admm = 1e-5, maxit_admm = 400, traceit = F)
                    OUT <- maxnorm_nonconvex_hybrid(Aobs = Aobs, omega = omega, what = what, Aini = Aini, L0 =L0, R0 = R0, beta = beta_seq[i], mu = mu_seq[j]* sc  * beta_seq[i], tau = tau,r = r, maxiter = maxiter + 200, retol = retol)
                    #OUTtemp = OUT
                  }else {
                    #if(j==1){r = Ranks[i-1,1]}else {
                    #   r = Ranks[i,j-1]
                    #}
                     if(i==1 && (j==2)){
                       Aini = OUT$Ahat
                        #L0 = OUT$L
                        #R0 = OUT$R
                       #Aini = NULL
                       OUT <- maxnorm_nonconvex_hybrid(Aobs = Aobs, omega = omega, what = what, Aini = Aini, L0 = NULL, R0 = NULL, beta = beta_seq[i], mu = mu_seq[j]* sc  * beta_seq[i], tau = tau,
                r = r, maxiter = maxiter, retol = retol)
                     }else
                     {
                       if(j==1){OUT = OUTtemp} 
                       Aini = OUT$Ahat
              #L0 = OUT$L
            #R0 = OUT$R
            #Aini = NULL
                       OUT <- maxnorm_nonconvex_hybrid(Aobs = Aobs, omega = omega, what = what, Aini = Aini, L0 = NULL, R0 = NULL, beta = beta_seq[i], mu = mu_seq[j]* sc  * beta_seq[i], tau = tau,
                r = r, maxiter = maxiter, retol = retol)
               
                     }
                  }
                temp = svd(OUT$Ahat)
                 ds = temp$d
                 rank = sum(ds/ds[1] > 1e-04 )
                 #rank = which(cumsum(ds)/sum(ds)>=0.9)[1]
                 Ranks[i,j] = rank
                cat(paste0("rank:", rank,"########"))
                sds = sum(ds[1:rank])
                cat(paste0("sds:", sds,"########"))
                if(i==1 && j==1){sds0 = sds}
                 if(rank == 1){
                   OUT$Ahat = temp$d[1] * (temp$u[,1:rank] %*% t(temp$v[,1:rank]))
                 }else{
                   OUT$Ahat =  temp$u[,1:rank] %*% diag(temp$d[1:rank]) %*%  t(temp$v[,1:rank])
                 }
                # error_out(OUT$Ahat, Aobstest, omegatest)
                cv_errors[i,j] =  norm((OUT$Ahat -
                Aobsvalidate)*omegavalidate, type = "F")/sqrt(sum(omegavalidate))
                 test_errors[i,j] =  norm((OUT$Ahat -
                Aobstest)*omegatest, type = "F")/sqrt(sum(omegatest))
              cat(paste0("val_error:", cv_errors[i,j],"########"))
              cat("\n")
              if(j>1 && (cv_errors[i,j] >= cv_errors[i,j-1]) && sds <= sds0){break}

              if(i== 1 && j==1){minerror = cv_errors[i,j]
               OUTfinal = OUT}else{
                if(cv_errors[i,j] < minerror){
                  OUTfinal = OUT 
                  minerror = cv_errors[i,j]
                }
              }
                 if(j==1){OUTtemp = OUT} 
                 sds0 = sds
              }
             
            }
             cvbalmax_grid = NULL
  select_tuning = list()

 
    minerror = min(cv_errors[cv_errors > 0])
    idx = which(cv_errors == minerror)[1]
    cvbalmax_grid = c(cvbalmax_grid, idx)
    tun_par = list()
    tun_par$beta = beta_seq[(idx -1 ) %% length(beta_seq) + 1]
    tun_par$mu = mu_seq[(idx -1) %/% length(beta_seq) + 1]  * tun_par$beta * sc
    select_tuning = tun_par
 
          
         return(list(cvbalmax_values = cv_errors, cvbalmax_grid = cvbalmax_grid, select_tuning = select_tuning, test_errors = test_errors, OUTfinal = OUTfinal, weights = what[omega==1], Ranks = Ranks))

  }
