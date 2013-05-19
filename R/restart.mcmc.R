#' @export
restart.mcmc <- function(SPM) {
  n.samples <- nrow(SPM$p.beta.theta.samples)
  n.batch <- ncol(SPM$tuning)
  batch.length <- n.samples/n.batch
  tuningvec <- exp(SPM$tuning[,n.batch])
  tuningvec.w <- exp(SPM$tuning.w[,n.batch])
  startingvec <- colMeans(SPM$p.beta.theta.samples[(n.samples-batch.length+1):
                                                     n.samples,])
  starting.w <- rowMeans(SPM$p.w.samples[,(n.samples-batch.length+1):n.samples])
  Alocs <- grep("K\\[", names(startingvec))
  philocs <- grep("phi\\[", names(startingvec))
  betalocs <- grep("\\(Intercept\\)", names(startingvec))
  starting=list(beta=startingvec[betalocs], phi=startingvec[philocs],
                A=startingvec[Alocs], w=starting.w)
  tuning=list(beta=tuningvec[betalocs], phi=tuningvec[philocs],
              A=tuningvec[Alocs], w=tuningvec.w)
  return(list(starting=starting, tuning=tuning))
  }