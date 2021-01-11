library(deSolve) # importation de la librarie diferential equation Solve 
library(phaseR)

rm(list = ls()) # Ré-initialisation de toutes les variables

SIRV <- function(t, x, parameters)
{
  m <- parameters[1] 
  v <- parameters[2]
  beta <- parameters[3]
  g <- parameters[4]
  mu <- parameters[5]
  
  S <- x[1]
  I <- x[2]
  R <- x[3]
  V <- x[4]
  
  dSdt <- (1-v)*m*N + (mu-beta*I)*V -m*S -beta*S*I # equation de S
  dIdt <- beta*S*I-m*I-g*I # equation de I
  dRdt <- g*I -m*R # equation de R
  dVdt <- v*m*N - m*V -(mu-beta*I)*V # equation de V
    
  list(c(dSdt, dIdt, dRdt, dVdt))
}

N = 1e06

m = 1/80 # esperance de vie 80ans
v = 0.8 # taux de vaccination 
g = 52/3 # durée de l'infection 3semaines
mu = 1/10

beta = 6.5*(m+g)/N

V0 = (v*m*N)/(m+mu)
S0 = ((1-v)*m*N + mu*V0)/m -1
I0 = 1
R0 = 0


param = c(m,v,beta,g,mu)
init = c(S0,I0,R0,V0)
times <- seq(0, 100, by = 0.001)

out <- ode(y = init, times = times, func = SIRV, parms = param)

S = out[,"1"]
I = out[,"2"]
R = out[,"3"]
V = out[,"4"]
t = out[,"time"]

par(mfrow=c(2,2))
plot(t,S, type = "l")
plot(t,I, type = "l")
plot(t,R, type = "l")
plot(t,V, type = "l")


par(mfrow=c(1,1))
Vvec=seq(from = 0, to = 1, by = 0.05)
Ivec = matrix(0,21,1)

for (k in 1:length(Vvec)){
  v0 <- Vvec[k]
  parameters <- c(m,v0,beta,g,mu)
  V0 = (v0*m*N)/(m+mu)
  S0 = ((1-v0)*m*N + mu*V0)/m -1
  I0 = 1
  R0 = 0
  initial.state <- c(S0, I0, R0, V0)
  times <- seq(0, 500, by = 0.01)
  out <- ode(y = initial.state, times = times, func = SIRV, parms = parameters, method = "ode45")
  I = out[,3]
  Ivec[k]= I[length(I)]
}

plot(Vvec,Ivec, type = "l", xlab = "v (%)", ylab = "I", sub = " Nombre d'infectés I à l'équilibre endémique en fonction de la fraction de population vaccinée v")
