rm(list=ls())
library(Biostrings)
setwd("C:/Users/hanic/Desktop/PRG/exercise_11")


viterbi <- function(obs, N, M, A, B, pi) {
  obs <- strsplit(as.character(obs), "")[[1]]
  T <- length(obs)
  n <- length(N)
  
  emit_idx <- match(obs, M)
  
  delta <- matrix(0, n, T)
  psi <- matrix(0, n, T)
  
  for (i in 1:n) {
    delta[i, 1] <- pi[i] + B[i, emit_idx[1]]
    psi[i, 1] <- 0
  }
  
  for (t in 2:T) {
    for (j in 1:n) {
      vals <- delta[, t-1] + A[, j]
      psi[j, t] <- which.max(vals)
      delta[j, t] <- max(vals) + B[j, emit_idx[t]]
    }
  }
  
  path <- integer(T)
  path[T] <- which.max(delta[, T])
  for (t in (T-1):1) {
    path[t] <- psi[path[t+1], t+1]
  }
  
  states <- N[path]
  
  return(list(delta = delta, path = states))
}



load("HMM1.RData")
N <- HMM1$N
M <- HMM1$M
A <- HMM1$A
B <- HMM1$B
pi <- HMM1$pi


str <- AAString("ACG")
res <- viterbi(str, N, M, A, B, pi)

res$delta
res$path
