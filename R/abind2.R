#'@export
#'@import abind
abind2 <- function(...) abind(..., rev.along=0, use.first.dimnames=TRUE)
