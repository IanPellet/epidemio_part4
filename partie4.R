library(deSolve) # importation de la librarie diferential equation Solve 
library(phaseR)

rm(list = ls()) # Ré-initialisation de toutes les variables

# Parametres du modèles
N = 1000000 # popultation 
c0 = 1/18 # gama passage à l'âge adulte 18ans 
m0 = 1/80 # taux de mortalité/natalité, essperance de vie de 80ans
v0 = 0.8 # couverture vaccinale de 80%
mu0 = 1/10 # durée de l'immunité vaccinale 10ans

# Modèle avec Enfants-Adultes
SIRage <-function(t, state, parameters) { 
  with(as.list(c(state, parameters)),{
    dSEdt <- (1 - v)* b * N - m * SE - beta * SE * (IE + IA) - c * SE 
    dSAdt <- c * SE - m * SA - beta * SA * (IE + IA) 
    dIEdt <- beta * SE * (IE + IA) - m * IE - g * IE - c * IE 
    dIAdt <- beta * SA * (IE + IA) + c * IE - m * IA - g * IA 
    dREdt <- v * b * N + g * IE - m * RE - c * RE 
    dRAdt <- g * IA + c * RE - m * RA 
    list(c(dSEdt, dSAdt, dIEdt, dIAdt, dREdt, dRAdt)) })
}


# Modèle avec immunité vaccinale
SIRvimm <-function(t, state, parameters) { 
  with(as.list(c(state, parameters)),{
    
    dVEdt <- v*m*N - m*VE - c*VE - (mu - beta*(IE+IA))*VE
    dVAdt <- c*VE - m*VA - (mu - beta*(IE+IA))*VA
    dSEdt <- (1-v)*m*N + (mu - beta*(IE+IA))*VE - m*SE - c*SE  - beta*SE*(IE + IA)
    dSAdt <- c*SE + (mu - beta*(IE+IA))*VA - m*SA - beta*SA*(IE + IA) 
    dIEdt <- beta*SE*(IE+IA) - m*IE - c*IE  - g*IE 
    dIAdt <- c*IE + beta*SA*(IE+IA) - m*IA - g*IA 
    dREdt <- g*IE - m*RE - c*RE 
    dRAdt <- g*IA + c*RE - m*RA 
    list(c(dVEdt, dVAdt, dSEdt, dSAdt, dIEdt, dIAdt, dREdt, dRAdt)) })
}


# Paramètres d'entrée du modèle 
# SIRage
parameters_age <- c(beta = 6.5*(52/3)/N, b = m0, m =m0, g=52/3-1/80, v = v0, c = c0)
# SIRvimm
parameters_Vimm <- c(beta = 6.5*(52/3)/N, b = m0, m =m0, g=52/3-1/80, v = v0, c = c0, mu = mu0)


# Conditions initiales 
# SIRage
initial.state_age <- c(SE = N*(1-v0)*m0/(m0+c0), SA = N*(1-v0)*c0/(m0+c0), IE = 1, IA = 1, RE = N*v0*m0/(m0+c0), RA = N*v0*c0/(m0+c0))
# SIRvimm
# Équilibre sans virus
VE0_Vimm <- v0*m0*N/(m0+c0+mu0)
VA0_Vimm <- c0*VE0_Vimm/(m0+mu0)
SE0_Vimm <- ((1-v0)*m0*N + mu0*VE0_Vimm)/(m0+c0)
SA0_Vimm <- (c0*SE0_Vimm + mu0*VA0_Vimm)/m0
# Condition initiale
initial.state_Vimm <- c(VE = VE0_Vimm, VA = VA0_Vimm, SE = SE0_Vimm-1, SA = SA0_Vimm-1, IE = 1, IA = 1, RE = 0, RA = 0)

# On simule les solutions de nos équations
times <- seq(0, 200, by = 0.1)
out_age <- ode(y = initial.state_age, times = times, func = SIRage, parms = parameters_age)
out_Vimm <- ode(y = initial.state_Vimm, times = times, func = SIRvimm, parms = parameters_Vimm)

# On récupère les variables
# SIRage
SEage = out_age[,"SE"] 
SAage = out_age[,"SA"] 
IEage = out_age[,"IE"] 
IAage = out_age[,"IA"] 
REage = out_age[,"RE"] 
RAage = out_age[,"RA"]
tage = out_age[,"time"] 
# SIRvimm
VEvimm = out_Vimm[,"VE"]
VAvimm = out_Vimm[,"VA"]
SEvimm = out_Vimm[,"SE"] 
SAvimm = out_Vimm[,"SA"] 
IEvimm = out_Vimm[,"IE"] 
IAvimm = out_Vimm[,"IA"] 
REvimm = out_Vimm[,"RE"] 
RAvimm = out_Vimm[,"RA"]
tvimm = out_Vimm[,"time"] 

# On trace les solutions
par(mfrow=c(3,5), mar=c(2,2,1,1))
plot(out_age, mfrow = NULL, mfcol = NULL, mar = NULL)
plot(out_Vimm, mfrow = NULL, mfcol = NULL)


# Solution de SIRage et SIRvimm sur un même graphique
t <- out_age[,1]
out_df <- data.frame(t, VEvimm, VAvimm, SEage, SEvimm, SAage, SAvimm, IEage, IEvimm, IAage, IAvimm, REage, REvimm, RAage, RAvimm)

out_names <- names(out_df)
par(mfrow = c(2,3),mar=c(4,4,1,1))
t_end <- length(t)
for (i in seq(4, (length(out_names)-1), by = 2)) {
  i1 <- i
  i2 <- i+1
  p1 <- out_df[,i1]
  p2 <- out_df[,i2]
  M <- max(p1,p2)
  m <- max(min(p1,p2),0)
  brnY <- c(m,M)
  name <- out_names[i]
  plot(t, p1, 'l', col = 'black', ylim = brnY, ylab = name)
  lines(t, p2, 'l', col = 'red')
  
}
legend("bottomright", legend = c('age','v imm'), col = c('black', 'red'), lty = c(1, 1))

# Zoom sur IA
par(mfrow=c(1,1), mar=c(4,4,1,1))
p1 <- out_df[,"IAage"]
p2 <- out_df[, "IAvimm"]
M <- 10000
m <- max(min(p1,p2),0)
brnY <- c(m,M)
plot(t, p1, 'l', col = 'black', ylim = brnY, ylab = "IA")
lines(t, p2, 'l', col = 'red')
legend("topright", legend = c('age','v imm'), col = c('black', 'red'), lty = c(1, 1))

