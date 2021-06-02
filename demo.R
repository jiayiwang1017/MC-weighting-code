#########################################
## Yobs: Observed matrix (zero for entries that are not observed); n1 by n2
## omega: observation indicator matrix; n1 by n2
## what: weight matrix; n1 by n2
Yobs = as.matrix(read.table("coat_train_data.txt"))
omega = 1 * (Yobs!=0)
n1 = nrow(Yobs)
n2 = ncol(Yobs)
what = matrix(1, n1, n2)

############ construct weights ##########
source("weights_calculation.R")

##########################################
##### hyperparameters: kappa, gamma ######
##########################################
## for example: 
kappa = 1e-08
gamma = Inf
fit <- MCB.SN(ind = omega, kaps = kappa, lower = 1, upper = gamma, traceit = FALSE, w0 = NULL) 
what[omega==1] = fit$w



#########################################
############# convex algorithm ##########
#########################################
Rcpp::sourceCpp("MC_cpp_hybrid2.cpp")
source(file = "balmax_prep_hybrid.R")
###################################################################################
## hyperparameters: R0 (max norm regularization), mu (nuclear norm regularization) #
####################################################################################
## for example: 
R0 = 10
mu = 30
alpha0 = max(Yobs)
Wini = rbind(cbind(diag(1, n1), Yobs), cbind(t(Yobs), diag(1, n2)))
Xini = 1.0 * Wini
Zini = diag(1, n1 + n2)
OUT = balmaxE_cpp_wraper(Aobs = Yobs, n1 = n1, n2= n2, omega = omega, what = what, R0 = R0, mu = mu, alpha0 = alpha0, rho = 0.1, W = Wini, X = Xini , Z = Zini, thresh_admm = 1e-05, maxit_admm = 500, traceit = F)
############# OUT$Ahat is the estimated target matrix ###########


#########################################
########### nonconvex algorithm #########
#########################################
source(file = "nonconvex_hybrid_prep_testval2.R")
###################################################################################
## hyperparameters: R0 (max norm regularization), mu (nuclear norm regularization) #
####################################################################################
## for example: 
R0 = 10
mu = 30
OUT = maxnorm_nonconvex_hybrid(Aobs = Yobs, omega = omega, what = what, Aini = NULL, L0 = NULL, R0 = NULL, beta = R0, mu = mu, tau = 1e-03, r = NULL, maxiter = 500, retol = 1e-05)
############# OUT$Ahat is the estimated target matrix ###########

