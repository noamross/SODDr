#' A dispersal function that only has 3 values, one for local, one for adjacent
#' cells, and zero for others.  Requires distances between cells to = 1
#' @export
adjacent.dispersal <- function(distance, local, adjacent) {
  ifelse(distance < 0.5, local,
         ifelse(distance < 1.5, adjacent,0)
  )
}