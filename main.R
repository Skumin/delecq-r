crossover_deleq <- function(mat, newmat, cr) {
  return((1 - cr) * mat + cr * newmat)
}

run_deleq <- function(fun, ..., boxbounds, Emat, constr, cr = 0.7, f.param = 0.9, maxgen = 100, NP, show.progress = TRUE) {
  gen <- 1
  mat <- gen_init_pop(NP = NP, boxbounds = boxbounds, Emat = Emat, constr = constr)
  funvals <- apply(X = mat, MARGIN = 1, FUN = fun, ...)
  rbest <- which.max(funvals)
  trugen <- maxgen
  while(gen <= trugen) {
    if(show.progress) {
      print(paste0('Generation ', gen))
    }
    pbest <- rbest
    pbest_ind <- mat[pbest, ]
    pbest_val <- fun(pbest_ind, ...)
    newmat <- mutate_deleq(mat = mat, boxbounds = boxbounds, fParam = f.param)
    newmat <- crossover_deleq(mat = mat, newmat = newmat, cr = cr)
    funvals1 <- apply(X = newmat, MARGIN = 1, FUN = fun, ...)
    mat[funvals1 > funvals, ] <- newmat[funvals1 > funvals, ]
    funvals <- ifelse(funvals1 > funvals, funvals1, funvals)
    rbest <- which.max(funvals)
    if(sum(mat[rbest, ]) != 1) { # If unfeasible
      if(show.progress) {
        print('* Unfeasible, projecting back.')
      }
      trugen <- trugen - 1
      if(gen <= trugen) {
        projmat <- project_population(mat = mat, Emat = Emat, constr = constr)
        funvals1 <- apply(X = projmat, MARGIN = 1, FUN = fun, ...)
        rbest <- which.max(funvals1)
        if(pbest_val > fun(projmat[rbest, ], ...)) {
          projmat[rbest, ] <- pbest_ind
          rbest <- pbest
        }
        mat <- projmat
      } else {
        mat[rbest, ] <- pbest_ind
      }
    }
    gen <- gen + 1
    if(show.progress) {
      print(paste0("** Value = ", funvals[rbest]))
    }
  }
  return(list(params = mat[rbest, ], value = fun(mat[rbest, ], ...), generations = trugen))
}
