


suppressWarnings(Rcpp::sourceCpp("MC_cpp_hybrid2.cpp", showOutput = FALSE, verbose = FALSE))

source(file = "balmax_prep_hybrid.R")
source(file = "nonconvex_hybrid_prep_testval2.R")


setwd(paste(wd, "/val_data/weights+hybrid+all3", sep = ""))
print(paste0("balancing_idx:",balancing_idx))



n1 = dim(Yobs)[1]
n2 = dim(Yobs)[2]

# R_balmax <- seq(0.2, 0.05, length = gridlength_balmax)
alpha0 <- max(abs(Yobs[Yobs != 0]))

#beta_seq <- alpha0 * exp(seq(log(0.5), log(2.5), length.out = 15))
beta_seq <- alpha0 * exp(seq(log(0.8), log(3), length.out = 15))
#beta_seq <- alpha0 * exp(seq(log(0.2), log(1), length.out = 15))
#R0_seq <- alpha0 * exp(seq(log(0.03), log(0.8), length.out = 15))
beta_seq <- sort(beta_seq, decreasing = TRUE)
cat("beta_seq:")
cat("\n")
print(beta_seq)

####  mu sequence ###
#mu_ind = 1/sqrt(sum(omega)*(sum(dim(omega))))
#mu_seq <-  exp(seq(log(5e-04), log(0.1), length.out = 15))
#mu_seq <-  exp(seq(log(1e-03), log(0.2), length.out = 15))
mu_seq <-  exp(seq(log(1e-03), log(0.2), length.out = 15))
mu_seq <- c(0, mu_seq)
#mu_seq <- sort( mu_seq, decreasing = TRUE)
cat("mu_seq:")
cat("\n")
print(mu_seq)


CV_OUT <- nonconvex_hybrid_validation(Aobs = Yobs, omega = omega, what = what, Aobsvalidate = Aobsvalidate, omegavalidate = omegavalidate, Aobstest = Aobstest, omegatest = omegatest,
  Aini = NULL, L0 = NULL, R0 = NULL, beta_seq, mu_seq, tau = 1e-03, r = NULL,
  maxiter = 500, retol = 1e-05)
save(CV_OUT, file = paste0(balancing_idx,"_cvout.RData"))

cv_errors <- CV_OUT$cvbalmax_values
test_errors <- CV_OUT$test_errors
cat("cv_errors:")
cat("\n")
print(cv_errors)
cat("\n")

cat("test_errors:")
cat("\n")
print(test_errors)
cat("\n")

idx = CV_OUT$cvbalmax_grid 
cat("final_test_errors:")
cat("\n")
print(test_errors[idx])
cat("\n")

SUM_W_OUT = CV_OUT$OUTfinal


error_wmax = error_out(Ahat = CV_OUT$OUTfinal$Ahat , Atest = Aobstest, rate = omegatest) 


print(error_wmax)


#save(SUM_W_OUT, file = paste0(balancing_idx,"_weights+hybrid_output.RData"))
setwd(wd)


# nohup Rscript ./nonconvex_hybrid.R </dev/null> ./combine_hybrid2/stdout 2>./combine_hybrid2/stderr &
