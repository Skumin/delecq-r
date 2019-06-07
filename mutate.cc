#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerVector sequenceVec(int n) {
  std::vector<int> ivec (n, 1);
  std::iota(ivec.begin(), ivec.end(), 0);
  return wrap(ivec);
}

// [[Rcpp::export]]
IntegerVector eraseInd(IntegerVector ivec, int eras) {
  IntegerVector newvec = clone(ivec);
  std::vector<int> newvec1 = as<std::vector<int> >(newvec);
  newvec1.erase (newvec1.begin()+eras);
  return wrap(newvec1);
}

// [[Rcpp::export]]
NumericMatrix mutateC(NumericMatrix mat, NumericMatrix boxbounds, double fParam) {
  NumericMatrix newmat(mat.nrow(), mat.ncol());
  int matsize = newmat.nrow();
  int parsize = newmat.ncol();
  IntegerVector rngg1 = sequenceVec(matsize);
  for (int i = 0; i < matsize; i++) {
    IntegerVector rngg = eraseInd(rngg1, i);
    int boundsok = 0;
    while (boundsok == 0) {
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
