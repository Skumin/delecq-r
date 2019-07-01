library(Rcpp)

sourceCpp('cpp_functions.cc')
source('main.R')

fn <- function(x) {
  -(x[1] - 2)^2 - (x[2] - 1)^2
}

microbenchmark::microbenchmark(run_deleq(fn, NP = 60, boxbounds = matrix(c(0, 0, 1, 1), nrow = 2), Emat = matrix(c(1, 1), nrow = 1), constr = matrix(1), cr = 0.7, f.param = 0.9, maxgen = 200, show.progress = FALSE))