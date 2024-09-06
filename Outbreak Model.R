
dS <- function(B, S, I){
  -B*S*I
}

dE <- function(B, S, I, a){
  -B*S*I - a*E
}

dI <- function(a, E, y, I){
  a*E - y*I
}

dR <- function(y, I){
  y*I
}

SEIR.model<-function (init, beta.s, gamma.e, gamma.i, times) {
  library(deSolve)
  seir <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS <- -beta.s * S * I
      dE <- beta.s * S * I - gamma.e * E
      dI <- gamma.e * E - gamma.i * I
      dR <- gamma.i * I
      return(list(c(dS, dE, dI, dR)))
    })
  }
  parameters <- c(beta.s = beta.s, gamma.e = gamma.e, gamma.i = gamma.i)
  out <- ode(y = init, times = times, func = seir, parms = parameters)
  out <- as.data.frame(out)
  out$time <- NULL
  return(out)
}

N <- 1000
parms <- c(beta=0.2,gamma.1=0.1,gamma.e=0.2)
times <- 1:550
xstart <- c(S=N-1, E=1, I=0, R=0)

p = 1/N
init <- c(S = 1-p, E = 0, I = p, R = 0)
mod <- SEIR.model(init = init, beta.s = parms[1], gamma.e = parms[3], gamma.i = parms[2], times = times)

SEIR_out <- round(mod*N, digits = 0)
