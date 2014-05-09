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
