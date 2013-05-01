#' A dispersal function that only has 3 values, one for local, one for adjacent
#' cells, and zero for others.  Requires distances between cells to = 1
#' @export
adjacent.dispersal <- function(distance, t1, t2, local, adjacent) {
  ifelse(distance < t1, local,
         ifelse(distance < t2, adjacent,0)
  )
}

#' @export
normal.dispersal <- function(distance, var=1, base=1, normalize=1) {
  base*dnorm(x=distance,mean=0,sd=var,log=FALSE)
}

#' @export
exp.dispersal <- function(distance, decay = 1, base=1, normalize=1) {
  base*exp(-distance*decay)
}