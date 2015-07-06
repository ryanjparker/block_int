#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}

// [[Rcpp::export]]
double sum_diag_mm_Rcpp(NumericMatrix A, NumericMatrix B) {
	int n = A.nrow();
	int i,j;
	double s = 0;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			s += A(i,j) * B(j,i);
		}
	}

	return(s);
}

// [[Rcpp::export]]
NumericMatrix ce_cov_Rcpp(NumericVector theta, NumericMatrix X) {
	int n      = X.nrow();
	int Ntheta = theta.size();

	NumericMatrix Sigma(n, n);

	int i,j,k;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			Sigma(i,j) = 0;
			for (k = 0; k < Ntheta; k++) {
				Sigma(i,j) += theta[k]*pow(X(i,k)-X(j,k), 2);
			}
			Sigma(i,j) = exp(-Sigma(i,j));
			Sigma(j,i) = Sigma(i,j);
		}
	}

	return(Sigma);
}

// [[Rcpp::export]]
NumericMatrix ce_cov_single_Rcpp(NumericVector theta, NumericMatrix X, NumericMatrix X2) {
	// compute Cov(X, X2)
	int n      = X.nrow();
	int Ntheta = theta.size();
	int n2     = X2.nrow();

	NumericMatrix Sigma(n, n2);

	int i,j,k;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n2; j++) {
			Sigma(i,j) = 0;
			for (k = 0; k < Ntheta; k++) {
				Sigma(i,j) += theta[k]*pow(X(i,k)-X2(j,k), 2);
			}
			Sigma(i,j) = exp(-Sigma(i,j));
		}
	}

	return(Sigma);
}

// [[Rcpp::export]]
NumericMatrix ce_partial_Rcpp(int e, NumericVector theta, NumericMatrix Sigma, NumericMatrix X) {
	int n      = X.nrow();
	int Ntheta = theta.size();

	NumericMatrix P(n, n);

	int i,j,k;
	k = e-1;

	double et = exp(theta[k]);

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			P(i,j) = -pow(X(i,k)-X(j,k), 2) * et * Sigma(i,j);
			P(j,i) = P(i,j);
		}
	}

	return(P);
}

// [[Rcpp::export]]
NumericMatrix ce_full_pred_X_cov(NumericMatrix X, NumericMatrix Xobs, NumericVector theta) {
	int Npred  = X.nrow();
	int Nobs   = Xobs.nrow();
	int Ntheta = theta.size();

	NumericMatrix Sigma(Npred, Nobs);

	int i,j,k;

	for (i = 0; i < Npred; i++) {
		for (j = 0; j < Nobs; j++) {
			Sigma(i,j) = 0;
			for (k = 0; k < Ntheta; k++) {
				Sigma(i,j) += theta[k]*pow(X(i,k)-Xobs(j,k), 2);
			}
			Sigma(i,j) = exp(-Sigma(i,j));
		}
	}

	return(Sigma);
}
