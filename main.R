# R Delecq
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
  if(!is.matrix(boxbounds)) {
    stop("'boxbounds' must be a matrix.")
  }
  if(ncol(boxbounds) != 2) {
    stop("'boxbounds' must have two columns.")
  }
  dm <- nrow(boxbounds)
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