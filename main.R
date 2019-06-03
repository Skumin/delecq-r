# R implementation of "An improved differential evolution algorithm for optimization 
# including linear equality constraints" by Barbosa, Bernardino, and Angelo
gen_init_pop <- function(E, const, NP, boxbounds) {
  if(!is.matrix(boxbounds)) {
    stop("'boxbounds' must be a matrix.")
  }
  if(ncol(boxbounds) != 2) {
    stop("'boxbounds' must have two columns.")
  }
  if(length(const) > 1 & !is.matrix(E)) {
    stop("'E' must be a matrix for this constraints vector.")
  }
  if(length(const) == 1 & !is.matrix(E)) {
    E <- t(as.matrix(E))
  }
  if(nrow(E) != length(const)) {
    stop("The number of constraints must be equal to the number of rows of 'E'.")
  }
  M <- E %*% t(E)
  if(length(const) == 1) {
    y <- const/M
    x0 <- t(E) %*% y
    dm <- ncol(E)
    xi <- matrix(0, ncol = dm, nrow = NP)
    for(i in seq_len(NP)) {
      boundsok <- FALSE
      while(!boundsok) {
        d <- as.matrix(runif(dm))
        z <- E %*% d
        u <- z/M
        v <- t(E) %*% u
        num <- x0 + d - v
        if(all(sapply(seq_along(num), function(i) num[i] >= boxbounds[i, 1] & num[i] <= boxbounds[i, 2]))) {
          xi[i, ] <- num
          boundsok <- TRUE
        }
      }
    }
  } else {
    stop("This is not yet implemented.")
  }
  return(xi)
}

gen_init_pop_simple <- function(NP, boxbounds) {
  dm <- nrow(boxbounds)
  if(is.null(NP)) {
    NP <- dm * 20
  }
  if(!is.matrix(boxbounds)) {
    stop("'boxbounds' must be a matrix.")
  }
  if(ncol(boxbounds) != 2) {
    stop("'boxbounds' must have two columns.")
  }
  xi <- matrix(0, ncol = dm, nrow = NP)
  for(i in seq_len(NP)) {
    boundsok <- FALSE
    while(!boundsok) {
      num <- runif(dm)
      num <- num/sum(num)
      if(all(sapply(seq_along(num), function(i) num[i] >= boxbounds[i, 1] & num[i] <= boxbounds[i, 2]))) {
        xi[i, ] <- num
        boundsok <- TRUE
      }
    }
  }
  return(xi)
}

mutate_delecq <- function(mat, boxbounds, f.param) {
  if(!is.matrix(mat)) {
    stop("'mat' must be a matrix.")
  }
  newmat <- matrix(0, ncol = ncol(mat), nrow = nrow(mat))
  for(i in seq_len(nrow(mat))) {
    boundsok <- FALSE
    while(!boundsok) {
      rs <- sample(seq_len(nrow(mat))[-i], size = 3, replace = FALSE)
      num <- mat[rs[1], ] + f.param * (mat[rs[2], ] - mat[rs[3], ])
      if(all(sapply(seq_along(num), function(i) num[i] >= boxbounds[i, 1] & num[i] <= boxbounds[i, 2]))) {
        newmat[i, ] <- num
        boundsok <- TRUE
      }
    }
  }
  return(newmat)
}

crossover_delecq <- function(mat, newmat, cr) {
  return((1 - cr) * mat + cr * newmat)
}

project_population_delecq <- function(mat) {
  E <- matrix(1, ncol = ncol(mat))
  M <- E %*% t(E)
  projmat <- matrix(0, ncol = ncol(mat), nrow = nrow(mat))
  for(i in seq_len(nrow(mat))) {
    num <- mat[i, ]
    z <- sum(num) - 1
    u <- z/M
    v <- t(E) %*% u
    projmat[i, ] <- as.vector(num - v)
  }
  return(projmat)
}

run_delecq <- function(fun, ..., boxbounds, cr = 0.5, f.param = 0.5, maxgen = 500, NP = NULL, show.progress = TRUE) {
  gen <- 1
  mat <- gen_init_pop_simple(NP = NP, boxbounds = boxbounds)
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
    newmat <- mutate_delecq(mat = mat, boxbounds = boxbounds, f.param = f.param)
    newmat <- crossover_delecq(mat = mat, newmat = newmat, cr = cr)
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
        projmat <- project_population_delecq(mat)
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
