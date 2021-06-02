#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include <stdlib.h>
#include <math.h>
#include <R_ext/Lapack.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

void FLTZ_weighted_cpp(mat &Wnew, const vec &Aobs, const int n1, const int n2, const vec &omega, const vec &what, const mat &WC, const double R0, const double alpha0, const double rho)
{
        // Aobs: observed Aobs vector (by column)
        //R0: maxi-norm
        //alpha0: maximum value (natrually be bounded by R0)
        //rho: admm hyperparameter
        //omega: flattend matrix T (by column)
        //what: non zero weights vector (by column)

        vec W12 = zeros<vec>(n1 * n2);

        // W11 submatrix
        // .submat( first_row, first_col, last_row, last_col )
        mat W11_old = WC.submat(0, 0, n1 - 1, n1 - 1);
        Wnew.submat(0, 0, n1 - 1, n1 - 1) = min(max(W11_old, -1 * R0 * ones<mat>(n1, n1)), R0 * ones<mat>(n1, n1));
        Wnew.submat(0, 0, n1 - 1, n1 - 1).diag() = min(max(W11_old.diag(), zeros<vec>(n1)), R0 * ones<vec>(n1));

        // W22 submatrix
        mat W22_old = WC.submat(n1, n1, n1 + n2 - 1, n1 + n2 - 1);
        Wnew.submat(n1, n1, n1 + n2 - 1, n1 + n2 - 1) = min(max(W22_old, -1 * R0 * ones<mat>(n2, n2)), R0 * ones<mat>(n2, n2));
        Wnew.submat(n1, n1, n1 + n2 - 1, n1 + n2 - 1).diag() = min(max(W22_old.diag(), zeros<vec>(n2)), R0 * ones<vec>(n2));

        // W12 submatrix
        // vector version for easy indexing
        vec W12_old = vectorise(WC.submat(0, n1, n1 - 1, n1 + n2 - 1));
        uvec unob = find(omega == 0);
        uvec ob = find(omega == 1);
        W12.elem(unob) = min(max(W12_old.elem(unob), -alpha0 * ones<vec>(unob.n_elem)), alpha0 * ones<vec>(unob.n_elem));
        W12.elem(ob) = min(max((Aobs % what + rho * W12_old.elem(ob)) / (what + rho * ones<vec>(ob.n_elem)), -alpha0 * ones<vec>(ob.n_elem)), alpha0 * ones<vec>(ob.n_elem));

        // fill back to Wnew
        Wnew.submat(0, n1, n1 - 1, n1 + n2 - 1) = reshape(W12, n1, n2);
        Wnew.submat(n1, 0, n1 + n2 - 1, n1 - 1) = reshape(W12, n1, n2).t();
}

// [[Rcpp::export]]
mat FLTZ_weighted_cpp_wraper(const vec &Aobs, const int n1, const int n2, const vec &omega, const vec &what, const mat &WC, const double R0, const double alpha0, const double rho)
{
        mat Wnew = zeros(size(WC));
        FLTZ_weighted_cpp(Wnew, Aobs, n1, n2, omega, what, WC, R0, alpha0, rho);
        return (Wnew);
}

vec balmaxE_cpp(const vec &Aobs, const int n1, const int n2, const vec &omega, const vec &what, const double R0, const double mu, const double alpha0, double *rho, mat &W, mat &X, mat &Z, const double thresh_admm, const int maxit_admm, bool traceit)
{
        // X: ensure p.s.d
        // Z: dual
        int iter;

        double obj1, obj2 = 0;
        vec errors = zeros<vec>(maxit_admm);
        mat X0 = 1.0 * X;
        mat Z0 = 1.0 * Z;
        mat W0 = 1.0 * W;
        uvec ob = find(omega == 1);
        vec W12 = zeros<vec>(n1 * n2);
        uvec posind;
        double resp, resd;
        int k = 0;
        int n = n1 + n2;
        int il, iu, m, lda = n, ldz = n, info, lwork;
        double abstol, vu, vl;
        double wkopt;
        double *work;
        int iwork[5 * n], ifail[n];
        //double w[n], z[n * n];
        double *w, *z;
        w = (double *)malloc(n * sizeof(double));
        z = (double *)malloc(n * n * sizeof(double));
        abstol = 1e-08;
        vec eigval;
        mat eigvec;
        mat D = diagmat(ones<vec>(n1 + n2));
        lwork = 8 * n;
        work = (double *)malloc(lwork * sizeof(double));
        mat tempX(n, n);
        //vu = INFINITY;
        vl = 1e-06;
        int flag = 1;
        for (iter = 0; iter < maxit_admm; iter++)
        {
                //update X
                if (flag == 0)
                {
                        eig_sym(eigval, eigvec, W0 - (1.0 / (*rho)) * (Z0 + mu * D));

                        posind = find(eigval > vl);
                        //std::cout << posind.n_elem << '\n';
                        //std::cout << eigval.max() << eigval.min() << '\n';
                        X = eigvec.cols(posind) * diagmat(eigval.elem(posind)) * (eigvec.cols(posind).t());
                        if (posind.n_elem < (n1 + n2) / 2)
                        {
                                flag = 1;
                        }
                }
                else
                {
                        tempX = W0 - (1.0 / (*rho)) * (Z0 + mu * D);
                        //vl = 0.0;
                        vu = 1 * n * tempX.max();
                        //vl = -1.0 * n * tempX.max();
                        /* Solve eigenproblem, only search for negative eigenvalues*/
                        dsyevx_("V", "V", "U", &n, tempX.memptr(), &lda, &vl, &vu, &il, &iu,
                                &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);
                        //std::cout << m << '\n';

                        vec eigval(w, m, /*copy_aux_mem =*/false, /*strict = */ false);
                        mat eigvec(z, n, m, false, false);
                        X = eigvec * diagmat(eigval) * (eigvec.t());
                        if (m > (n1 + n2) / 2)
                        {
                                flag = 0;
                        }
                }

                //update W
                FLTZ_weighted_cpp(W, Aobs, n1, n2, omega, what, X + (1.0 / (*rho)) * Z0, R0, alpha0, (*rho));
                // update Z
                // why times 1.618?
                Z = Z0 + 1.618 * (*rho) * (X - W);

                // measurements
                resp = norm(X - W, "fro");
                resd = std::max(norm((*rho) * (W0 - W) + (1.618 - 1.0) * (*rho) * (X - W), "fro"),
                                norm((1.618 - 1.0) * (*rho) * (X - W), "fro"));
                if (resp < 0.5 * resd & k >= 10)
                {
                        (*rho) = 0.7 * (*rho);
                        k = 0;
                        //std::cout << rho << '\n';
                }
                if (resd < 0.5 * resp & k >= 10)
                {
                        (*rho) = 1.3 * (*rho);
                        k = 0;
                        //std::cout << rho << '\n';
                }
                //std::cout << error << '\n';
                errors(iter) = std::max(resp, resd);

                if (traceit)
                {
                        W12 = vectorise(W.submat(0, n1, n1 - 1, n1 + n2 - 1));
                        obj1 = sum(square(what % (W12.elem(ob) - Aobs))) / (n1 * n2);

                        //std::cout << omega.n_elem << what.n_elem << obW12vec.n_elem << Aobs.n_elem << '\n';
                        mat temp = (W - X) % (Z0.t());
                        obj2 = obj1 + 2 * sum(temp.diag()) + (*rho) * accu(square(W - X));
                        std::cout << "iter:" << iter << "obj1:" << obj1 << "obj2:" << obj2 << "error:" << errors[iter] << '\n';
                }
                //std::cout << errors(iter) << '\n';
                if (errors(iter) < thresh_admm)
                {
                        break;
                }
                W0 = 1.0 * W;
                Z0 = 1.0 * Z;
                X0 = 1.0 * X;
                k = k + 1;
        }

        errors = errors.head(iter + 1);
        //std::cout << iter << '\n';
        return (errors);
}

// [[Rcpp::export]]
List balmaxE_cpp_wraper(const mat &Aobs, const int n1, const int n2, const mat &omega, const mat &what, const double R0, const double mu, const double alpha0, double rho, mat &W, mat &X, mat &Z, const double thresh_admm, const int maxit_admm, bool traceit)
{
        // match arguments with original R function
        vec errors;
        //mat W = 1.0 * Wini;
        //mat Z = 1.0 * Zini;
        //mat X = 1.0 * Xini;
        uvec ob = find(vectorise(omega) == 1);

        errors = balmaxE_cpp(Aobs.elem(ob), n1, n2, vectorise(omega), what.elem(ob), R0, mu, alpha0, &rho, W, X, Z, thresh_admm, maxit_admm, traceit);
        List OUT;
        OUT["Ahat"] = W.submat(0, n1, n1 - 1, n1 + n2 - 1);
        OUT["errors"] = errors;
        OUT["W"] = W;
        OUT["Z"] = Z;
        OUT["X"] = X;
        OUT["rho"] = rho;
        return (OUT);
}

// [[Rcpp::export]]
cube balmaxE_cv_cpp(const mat &A0, const vec &Aobstrain, const vec &Aobsvalidate, const vec &omegatrain, const vec &omegavalidate, const int n1, const int n2, const vec &whattrain, const vec &whatvalidate, const vec &R0_seq, const vec &mu_seq, const int alpha0, const double thresh_admm, const int maxit_admm, bool traceit, double rho, mat &W, mat &X, mat &Z)
{
        int i, j;

        mat Ahat = zeros<mat>(n1, n2);
        vec errors;
        double R0;
        uvec obvalidate = find(omegavalidate == 1);
        vec Ahatvalidate = zeros<vec>(obvalidate.n_elem);

        cube cv_errors = zeros<cube>(R0_seq.n_elem, mu_seq.n_elem, 3);
        //mat cv_errors_uni = zeros<mat>(R0_seq.n_elem, mu_seq.n_elem);
        //mat cv_errors_weighted = zeros<mat>(R0_seq.n_elem, mu_seq.n_elem);
        //mat cv_errors_best = zeros<mat>(R0_seq.n_elem, mu_seq.n_elem);

        for (i = 0; i < R0_seq.n_elem; i++)
        {

                //R0 = alpha0 * sqrt(std::min(n1, n2)) * R0_ratio(i);
                // mu times R0
                for (j = 0; j < mu_seq.n_elem; j++)
                {
                        rho = 0.1;
                        //std::cout << i << j << '\n';
                        //mu_seq is scaled, but is not multiplicated by R0
                        if (alpha0 > R0_seq(i))
                        {

                                errors = balmaxE_cpp(Aobstrain, n1, n2, omegatrain, whattrain, R0_seq(i), mu_seq(j) * R0_seq(i), R0_seq(i), &rho, W, X, Z, thresh_admm, maxit_admm, traceit);
                        }
                        else
                        {
                                errors = balmaxE_cpp(Aobstrain, n1, n2, omegatrain, whattrain, R0_seq(i), mu_seq(j) * R0_seq(i), alpha0, &rho, W, X, Z, thresh_admm, maxit_admm, traceit);
                        }

                        //std::cout << errors.n_elem << '\n';

                        Ahat = W.submat(0, n1, n1 - 1, n1 + n2 - 1);

                        Ahatvalidate = Ahat.elem(obvalidate);

                        // uni
                        cv_errors(i, j, 0) = norm(Ahatvalidate - Aobsvalidate) / sqrt(sum(omegavalidate));
                        // weighted
                        cv_errors(i, j, 1) = norm((Ahatvalidate - Aobsvalidate) % sqrt(whatvalidate)) / sqrt(sum(omegavalidate));
                        // true
                        cv_errors(i, j, 2) = norm(Ahat - A0) / sqrt(n1 * n2);

                        if (j > 0 && (cv_errors(i, j, 0) > cv_errors(i, j - 1, 0)) && (cv_errors(i, j, 1) > cv_errors(i, j - 1, 1)))
                        {
                                break;
                        }
                }
        }
        return (cv_errors);
}

// [[Rcpp::export]]
mat balmaxE_true_cpp(const mat &Atrue, const vec &Aobs, const vec &omega, const int n1, const int n2, const vec &what, const vec &R0_seq, const vec &mu_seq, const int alpha0, const double thresh_admm, const int maxit_admm, bool traceit, double rho, const mat &Wini, const mat &Xini, const mat &Zini)
{
        int i, j;

        mat Ahat = zeros<mat>(n1, n2);
        vec errors;
        double R0;

        mat W = 1.0 * Wini;
        mat Z = 1.0 * Zini;
        mat X = 1.0 * Xini;

        mat true_errors = zeros<mat>(R0_seq.n_elem, mu_seq.n_elem);
        //mat cv_errors_uni = zeros<mat>(R0_seq.n_elem, mu_seq.n_elem);
        //mat cv_errors_weighted = zeros<mat>(R0_seq.n_elem, mu_seq.n_elem);
        //mat cv_errors_best = zeros<mat>(R0_seq.n_elem, mu_seq.n_elem);

        for (i = 0; i < R0_seq.n_elem; i++)
        {

                //R0 = alpha0 * sqrt(std::min(n1, n2)) * R0_ratio(i);
                for (j = 0; j < mu_seq.n_elem; j++)
                {
                        //W = 1.0 * Wini;
                        //Z = 1.0 * Zini;
                        //X = 1.0 * Xini;
                        std::cout << i << j << '\n';
                        rho = 0.1;
                        //std::cout << alpha0 << R0_seq(i) << '\n';
                        //mu_seq is scaled, but is not multiplicated by
                        if (alpha0 > R0_seq(i))
                        {
                                errors = balmaxE_cpp(Aobs, n1, n2, omega, what, R0_seq(i), mu_seq(j) * R0_seq(i), R0_seq(i), &rho, W, X, Z, thresh_admm, maxit_admm, traceit);
                        }
                        else
                        {
                                errors = balmaxE_cpp(Aobs, n1, n2, omega, what, R0_seq(i), mu_seq(j) * R0_seq(i), alpha0, &rho, W, X, Z, thresh_admm, maxit_admm, traceit);
                        }

                        std::cout << errors.n_elem << '\n';

                        Ahat = W.submat(0, n1, n1 - 1, n1 + n2 - 1);

                        true_errors(i, j) = norm(Ahat - Atrue, "fro") / sqrt(n1 * n2);

                        if (j > 0 && (true_errors(i, j) > true_errors(i, j - 1)))
                        {
                                break;
                        }
                }
        }
        return (true_errors);
}
