#include <Rcpp.h>
using namespace Rcpp;

template <int RTYPE>
inline Rcpp::Vector<RTYPE> 
anti_subset(const Rcpp::Vector<RTYPE>& x, Rcpp::IntegerVector idx) {
  Rcpp::IntegerVector xi = Rcpp::seq(0, x.size() - 1);
  return x[Rcpp::setdiff(xi, idx)];
}

// [[Rcpp::export]]
Rcpp::IntegerVector AntiSubset(Rcpp::IntegerVector x, Rcpp::IntegerVector idx) {
  return anti_subset(x, idx); // Write some commment here
}

// [[Rcpp::export]]
NumericMatrix mutateC(NumericMatrix mat, NumericMatrix boxbounds, double fParam) {
  NumericMatrix newmat(mat.nrow(), mat.ncol());
  int matsize = newmat.nrow();
  int parsize = newmat.ncol();
  for (int i = 0; i < matsize; i++) {
    int boundsok = 0;
    while (boundsok == 0) {
      IntegerVector rngg = seq(0, matsize - 1);
      rngg = AntiSubset(rngg, i);
      IntegerVector rs = sample(rngg, 3, false);
      NumericVector num = mat(rs[0], _ ) + fParam * (mat(rs[1], _ ) - mat(rs[2], _ ));
      IntegerVector tst(parsize);
      for (int j = 0; j < parsize; j++) {
        if (num[j] >= boxbounds(j, 0) && num[j] <= boxbounds(j, 1)) {
          tst[j] = 1;
        }
      }
      if (sum(tst) == parsize) {
        boundsok = 1;
        newmat(i, _) = num;
      }
    }
  }
  return newmat;
}