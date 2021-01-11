library(deSolve) # importation de la librarie diferential equation Solve 
library(phaseR)

rm(list = ls()) # Ré-initialisation de toutes les variables
#mIA <- matrix(data = NA, nrow = 11, ncol = 2)
# Parametres du modèles
N = 1000000 # popultation 
c0 = 1/18 # gama passage à l'âge adulte 18ans 
m0 = 1/80 # taux de mortalité/natalité, essperance de vie de 80ans
v0 = 0.7 # couverture vaccinale de 80%
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
SIRVage <-function(t, state, parameters) { 
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

# même modèle, pour vérifier d'éventuelles fautes de frappe
SIRVage2 <- function(t, x, parameters)
{
  m <- parameters[1] 
  v <- parameters[2]
  beta <- parameters[3]
  g <- parameters[4]
  mu <- parameters[5]
  c <- parameters[6]
  
  SE <- x[1]
  SA <- x[2]
  IE <- x[3]
  IA <- x[4]
  RE <- x[5]
  RA <- x[6]
  VE <- x[7]
  VA <- x[8]
  
  dSEdt <- (1-v)*m*N + (mu-beta*(IA+IE))*VE - m*SE - beta*SE*(IA+IE) - c*SE # equation de SE
  dSAdt <- (mu-beta*(IA+IE))*VA - m*SA - beta*SA*(IA+IE) + c*SE # equation de SA
  dIEdt <- beta*SE*(IA+IE) - m*IE - g*IE - c*IE # equation de IE
  dIAdt <- beta*SA*(IA+IE) - m*IA - g*IA + c*IE # equation de IA
  dREdt <- g*IE - m*RE - c*RE # equation de RE
  dRAdt <- g*IA - m*RA + c*RE # equation de RA
  dVEdt <- v*m*N - m*VE - (mu-beta*(IA+IE))*VE - c*VE # equation de VE
  dVAdt <- c*VE - m*VA - (mu-beta*(IA+IE))*VA # equation de VA
  
  list(c(dSEdt, dSAdt, dIEdt, dIAdt, dREdt, dRAdt, dVEdt, dVAdt))
}


# Paramètres d'entrée du modèle 
# SIRage
parameters_age <- c(beta = 6.5*(52/3)/N, b = m0, m =m0, g=52/3-1/80, v = v0, c = c0)
# SIRVage
parameters_Vimm <- c(beta = 6.5*(52/3)/N, b = m0, m =m0, g=52/3-1/80, v = v0, c = c0, mu = mu0)
# SIRVage2
parameters_SIRVage2 <- c(m0, v0, 6.5*(52/3)/N, 52/3, mu0, c0)


# Conditions initiales 
# SIRage
initial.state_age <- c(SE = N*(1-v0)*m0/(m0+c0), SA = N*(1-v0)*c0/(m0+c0), IE = 1, IA = 1, RE = N*v0*m0/(m0+c0), RA = N*v0*c0/(m0+c0))
# SIRVage
# Équilibre sans virus
VE0_Vimm <- v0*m0*N/(m0+c0+mu0)
VA0_Vimm <- c0*VE0_Vimm/(m0+mu0)
SE0_Vimm <- ((1-v0)*m0*N + mu0*VE0_Vimm)/(m0+c0)
SA0_Vimm <- (c0*SE0_Vimm + mu0*VA0_Vimm)/m0
# Condition initiale
initial.state_Vimm <- c(VE = VE0_Vimm, VA = VA0_Vimm, SE = SE0_Vimm-1, SA = SA0_Vimm-1, IE = 1, IA = 1, RE = 0, RA = 0)
# SIRVage2
VE0 <- v0*m0*N/(m0+c0+mu0)
VA0 <- c0*VE0/(m0+mu0)
SE0 <- ((1-v0)*m0*N + mu0*VE0)/(m0+c0)
SA0 <- (c0*SE0 + mu0*VA0)/m0
initial.state_SIRVage2 <- c(SE0-1, SA0-1, 1, 1, 0, 0, VE0, VA0)

# On simule les solutions de nos équations
times <- seq(0, 300, by = 0.1)
out_age <- ode(y = initial.state_age, times = times, func = SIRage, parms = parameters_age, method = "ode45")
out_Vimm <- ode(y = initial.state_Vimm, times = times, func = SIRVage, parms = parameters_Vimm, method = "lsodes")
out_SIRVage2 <- ode(y = initial.state_SIRVage2, times = times, func = SIRVage2, parms = parameters_SIRVage2)

# On récupère les variables
# SIRage
SEage = out_age[,"SE"] 
SAage = out_age[,"SA"] 
IEage = out_age[,"IE"] 
IAage = out_age[,"IA"] 
REage = out_age[,"RE"] 
RAage = out_age[,"RA"]
tage = out_age[,"time"] 
# SIRVage
VEvimm = out_Vimm[,"VE"]
VAvimm = out_Vimm[,"VA"]
SEvimm = out_Vimm[,"SE"] 
SAvimm = out_Vimm[,"SA"] 
IEvimm = out_Vimm[,"IE"] 
IAvimm = out_Vimm[,"IA"] 
REvimm = out_Vimm[,"RE"] 
RAvimm = out_Vimm[,"RA"]
tvimm = out_Vimm[,"time"] 
# SIRVage2
VE_SIRVage2 = out_SIRVage2[,"7"]
VA_SIRVage2 = out_SIRVage2[,"8"]
SE_SIRVage2 = out_SIRVage2[,"1"] 
SA_SIRVage2 = out_SIRVage2[,"2"] 
IE_SIRVage2 = out_SIRVage2[,"3"] 
IA_SIRVage2 = out_SIRVage2[,"4"] 
RE_SIRVage2 = out_SIRVage2[,"5"] 
RA_SIRVage2 = out_SIRVage2[,"6"]
t_SIRVage2 = out_SIRVage2[,"time"] 

# On trace les solutions
par(mfrow=c(3,5), mar=c(2,2,1,1))
plot(out_age, mfrow = NULL, mfcol = NULL, mar = NULL)
plot(out_Vimm)
plot(out_SIRVage2)

# Solution de SIRage et SIRVage sur un même graphique
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
# on peut observer une plus haure valeure de IA pour le modèle avec immunité vaccinale limitée

# IA en fonction de la couverture vaccinale SIRVage
Vvec=seq(from = 0, to = 1, by = 0.05)
IAvec = matrix(0,21,1)
IEvec = matrix(0,21,1)
Ivec = matrix(0,21,1)

# mu0 <- 1/30 # test avec différentes durrées d'immunité vaccinale

methods = c("lsoda", "lsode", "lsodes","lsodar","vode", "daspk", "euler", "rk4", "ode23", "ode45", "radau", "bdf", "bdf_d", "adams", "impAdams", "impAdams_d")
used_meth = c()

for (k in 1:length(Vvec)){ 
  vi=Vvec[k]
  
  parameters_Vimm_i <- c(beta = 6.5*(52/3)/N, b = m0, m =m0, g=52/3-1/80, v = vi, c = c0, mu = mu0)
  
  # Équilibre sans virus
  VE0_Vimm_i <- vi*m0*N/(m0+c0+mu0)
  VA0_Vimm_i <- c0*VE0_Vimm_i/(m0+mu0)
  SE0_Vimm_i <- ((1-vi)*m0*N + mu0*VE0_Vimm_i)/(m0+c0)
  SA0_Vimm_i <- (c0*SE0_Vimm_i + mu0*VA0_Vimm_i)/m0
  # Condition initiale
  initial.state_Vimm_i <- c(VE = VE0_Vimm_i, VA = VA0_Vimm_i, SE = SE0_Vimm_i-1, SA = SA0_Vimm_i-1, IE = 1, IA = 1, RE = 0, RA = 0)
  times <- seq(0, 500, by = 0.1)
  
  # beaucoups de problèmes avec les différentes méthodes de résolution pour la fonction ode
  # on teste donc les différentes méthodes jusqu'à avoir une valeur de IA supperieure à 0 
  for (meth in methods) {
    out <- ode(y = initial.state_Vimm_i, times = times, func = SIRVage, parms = parameters_Vimm_i, method = meth)
    #SE = out[,"SE"] 
    #SA = out[,"SA"] 
    IE = out[,"IE"] 
    IA = out[,"IA"]  
    #RE = out[,"RE"] 
    #RA = out[,"RA"] 
    #t = out[,"time"] 

    IAvec[k]= max(IA[length(IA)],0)
    IEvec[k]= max(IE[length(IE)],0)
    Ivec[k] = IAvec[k] + IEvec[k]
    if (IAvec[k] > 10 && !is.nan(IAvec[k])) {
      used_meth <- c(used_meth, meth) # on récupère les méthodes ayant fonctionné 
      break
    }
  }
}

# On trace pour finir comment le rapport IA/N dépend de v:
par(mfrow=c(1,1))

# fraction de la pop dans I en fonction de la couverture vaccinale 
M <- max(Ivec, IAvec, IEvec)/N
m <- min(Ivec, IAvec, IEvec)/N
plot(Vvec,Ivec/N, type = "l", ylim = c(m,M)) 
lines(Vvec, IAvec/N, col = 2)
lines(Vvec, IEvec/N, col = 3)
legend("topright", legend = c('I/N','IA/N', 'IE/N'), col = c(1,2,3), lty = c(1, 1, 1))

# fraction de la pop non vaccinée dans I en fonction de la couverture vaccinale 
plot(Vvec,(Ivec/N)/(1-Vvec), type = "l")
lines(Vvec, (IAvec/N)/(1-Vvec), col = 2)
lines(Vvec, (IEvec/N)/(1-Vvec), col = 3)
legend("topleft", legend = c('(I/N)/(1-Vvec)','(IA/N)/(1-Vvec)', '(IE/N)/(1-Vvec)'), col = c(1,2,3), lty = c(1, 1, 1))


# la fraction de la pop dans IA augmente avec la couvereture vaccinale -> il ne faudrait donc pas vacciner ??? C'est louche ! 
# il y a peut-être un pb dans la conception du modèle 
# dans SIRV sans les classes d'âge on a bien une décroissance de I avec l'augmentation de la couverture vaccinale
# est-ce que lorqu'on vaccine à la naissance on reste immunisé pendant l'enfance puis BIM! quand on est adulte on est baisé ???
# plot IE/N et I/N
# On a bien I=f(v) décroissante comme dans SIRV
# comparer à SIRage 
# Test avec différentes valeurs de mu
# avec mu = 0 on devrais avoir les mêmes résultats que SIRage, ce n'est pas DU TOUT le cas ! -> revoir le modèle 




# IA en fonction de la couverture vaccinale SIRVage2
IAvec2 = matrix(0,21,1)

methods = c("lsoda", "lsode", "lsodes","lsodar","vode", "daspk", "euler", "rk4", "ode23", "ode45", "radau", "bdf", "bdf_d", "adams", "impAdams", "impAdams_d")
used_meth2 = c()

for (k in 1:length(Vvec)){ 
  vi=Vvec[k]
  parameters_SIRVage2_i <- c(m0, vi, 6.5*(52/3)/N, 52/3, mu0, c0)

  # SIRVage2
  VE0_i <- vi*m0*N/(m0+c0+mu0)
  VA0_i <- c0*VE0_i/(m0+mu0)
  SE0_i <- ((1-vi)*m0*N + mu0*VE0_i)/(m0+c0)
  SA0_i <- (c0*SE0_i + mu0*VA0_i)/m0
  initial.state_SIRVage2_i <- c(SE0_i-1, SA0_i-1, 1, 1, 0, 0, VE0_i, VA0_i)
  times <- seq(0, 500, by = 0.1)
  
  # beaucoups de problèmes avec les différentes méthodes de résolution pour la fonction ode
  # on teste donc les différentes méthodes jusqu'à avoir une valeur de IA supperieure à 0 
  for (meth in methods) {
    out <- ode(y = initial.state_SIRVage2_i, times = times, func = SIRVage2, parms = parameters_SIRVage2_i, method = meth)
    #SE <- out[,"1"] 
    #SA <- out[,"2"] 
    #IE <- out[,"3"] 
    IA <- out[,"4"]  
    #RE <- out[,"5"] 
    #RA <- out[,"6"]
    #VE <- out[,"7"]
    #VA <- out[,"8"]
    #t = out[,"time"] 
    
    IAvec2[k]= max(IA[length(IA)],0)
    if (IAvec2[k] > 10 && !is.nan(IAvec2[k])) {
      used_meth2 <- c(used_meth2, meth) # on récupère les méthodes ayant fonctionné 
      break
    }
  }
}

# On trace pour finir comment le rapport IA/N dépend de v:
par(mfrow=c(1,1))
plot(Vvec,IAvec2/N, type = "l") # fraction de la pop dans IA en fonction de la couverture vaccinale 
plot(Vvec,(IAvec2/N)/(1-Vvec), type = "l") # fraction de la pop non vaccinée dans IA en fonction de la couverture vaccinale 
used_meth2
# on obtient les même résultats que pour SIRVage, le problème ne viens donc surement pas d'une faute de frappe dans le modèle