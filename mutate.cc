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

// [[Rcpp::export]]
NumericMatrix project_populationC(NumericMatrix mat) {
  NumericMatrix Emat(1, mat.ncol());
  int Mmat = mat.ncol();
  NumericMatrix projmat(mat.nrow(), mat.ncol());
  for (int i = 0; i < mat.nrow(); i++) {
    NumericVector num = mat(i, _ );
    double z = sum(num) - 1;
    double u = z/Mmat;
    NumericMatrix v = transpose(Emat) * u;
    NumericVector v1 = v( _, 0);
    projmat(i, _ ) = num - v1;
  }
  return projmat;
}

// [[Rcpp::export]]
NumericMatrix gen_init_popC(int NP, NumericMatrix boxbounds) {
  int dm = boxbounds.nrow();
  NumericMatrix xi(NP, dm);
  for (int i = 0; i < NP; i++) {
    int boundsok = 0;
    while (boundsok == 0) {
      NumericVector nums(dm);
      for (int j = 0; j < dm; j++) {
        nums[j] = R::runif(0, 1) * (boxbounds(j, 1) - boxbounds(j, 0)) + boxbounds(j, 0);
      }
      nums = nums/sum(nums);
      IntegerVector ids(dm);
      for (int j = 0; j < dm; j++) {
        if (nums[j] >= boxbounds(j, 0) && nums[j] <= boxbounds(j, 1)) {
          ids[j] = 1;
        }
      }
      if (sum(ids) == ids.size()) {
        xi(i, _) = nums;
        boundsok = 1;
      }
    }
  }
  return xi;
}
