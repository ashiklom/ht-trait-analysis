# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat rmvnorm(int n, const arma::vec& mu, const arma::mat& Sigma) {
    arma::mat L = arma::chol(Sigma, "lower");
    arma::mat Z = arma::randn(mu.n_rows, n);
    arma::mat LZ = L * Z;
    LZ.each_col() += mu;
    return LZ;
}

arma::mat sample_reflectance(const arma::colvec& refl,
                             const arma::colvec& unc,
                             const arma::mat& cormat,
                             int n) {
    arma::mat dunc = arma::diagmat(unc);
    arma::mat covmat = dunc * cormat * dunc;
    return rmvnorm(n, refl, covmat);
}

// [[Rcpp::export]]
arma::cube sample_reflectance_array(const arma::mat& refl,
                                    const arma::mat& unc,
                                    const arma::mat& corr,
                                    int n) {
    int nwl = refl.n_rows;
    int nobs = refl.n_cols;
    arma::cube result(nwl, n, nobs);
    for (int i = 0; i < nobs; i++) {
        result.slice(i) = sample_reflectance(refl.col(i), unc.col(i), corr, n);
    }
    return result;
}
