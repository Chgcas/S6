rm(list=ls())

n<-10 # Nombre de variables utilisées pour simuler une fois A_{t,n}
m<-100 # Nombre de simulations de A_{t,n}
t_max<-100 # Nombre de pas de temps
estim <-1000 # Nombre de variables utilisées pour calculer E(abs(X)^3)

# On estime E(abs(X)^3) pour X suivant une loi Stud(4) et X suivant une loi N(0,1)

estim_moment3 <- function(estim){
  phi1 <- rnorm(estim)
  phi2 <- rnorm(estim) / (sqrt(rchisq(n,4)/4))
  return(data.frame(moment3norm=mean(abs(phi1)^3), moment3stud4=mean(abs(phi2)^3)))
}

moment3norm=estim_moment3(estim)$moment3norm
moment3stud4=estim_moment3(estim)$moment3stud4

# ------------- Mouvement brownien -------------

brown <- function(t_max) {
  tps <- 1:t_max
  tps <- tps[-1] - tps[-t_max] # On construit la suite tn - tn-1
  w <- rnorm(t_max - 1, 0, sqrt(tps))
  return(c(0, cumsum(w)))
}

bt <- brown(t_max)

# ------------- Approximation avec des lois N(0,1) -------------

brown_norm <- function(n, t, moment3stud4, moment3norm){ # t peut être un réel ou un vecteur
  nt <- round(n*t) # Fonction partie entière
  echantillon_normale <- rnorm(nt[length(nt)])
  x <- echantillon_normale / sqrt(n)
  cumx <- cumsum(x)
  return(data.frame(echantillon = cumx[nt], borne_sup = 3*moment3norm / sqrt(nt)))
  # On retourne une trajectoire utilisant des lois normales et une borne de Berry-Esseen 3*rho/(sigma^3 * sqrt(nt))
}

b_norm <- brown_norm(n, c(1:t_max), moment3stud4, moment3norm)

# ------------- Approximation avec des lois Student(3) -------------

brown_stud3 <- function(n, t, moment3stud4, moment3norm){
  nt <- round(n*t) # Fonction partie entière
  x <- (1 / sqrt(n)) * rnorm(nt[length(nt)]) / (sqrt(rchisq(nt[length(nt)], 3) / 3))
  # x contient nt[length(nt)] simulations de la loi Student(3)
  cumx <- cumsum(x)
  sigma <- sqrt(3) # sqrt(k / (k-2))
  return(data.frame(echantillon=cumx[nt] / sigma, borne_sup = 1))
}

b_stud3 <- brown_stud3(n, c(1:t_max), moment3stud4, moment3norm)

# ------------- Approximation avec des lois Student(4) -------------

brown_stud4 <- function(n, t, moment3stud4, moment3norm){
  nt <- round(n*t) # Fonction partie entière
  echantillon_student <- rnorm(nt[length(nt)]) / (sqrt(rchisq(nt[length(nt)], 4) / 4))
  x <- (1 / sqrt(n)) * echantillon_student
  # x contient nt[length(nt)] simulations de la loi Student(4)
  cumx <- cumsum(x)
  sigma <- sqrt(2) # sqrt(k / (k-2))
  return(data.frame(echantillon = cumx[nt] / sigma, borne_sup = 3*moment3stud4 / (sigma^3*sqrt(nt))))
}

b_stud4 <- brown_stud4(n, c(1:t_max), moment3stud4, moment3norm)

# ------------- Approximation avec des lois Rademacher -------------

brown_rade <- function(n, t, moment3stud4, moment3norm){
  nt <- round(n*t) # Fonction partie entière
  x <- (2*rbinom(nt[length(nt)], 1, 1/2) - 1) / sqrt(n)
  # x contient nt[length(nt)] simulations de la loi Rademacher
  cumx <- cumsum(x)
  return(data.frame(echantillon = cumx[nt], borne_sup = 1/sqrt(nt)))
}

b_rade <- brown_rade(n, c(1:t_max), moment3stud4, moment3norm)

# ------------- Approximation avec des lois uniformes -------------

brown_unif <- function(n, t, moment3stud4, moment3norm){
  nt <- round(n*t) # Fonction partie entière
  x <- runif(nt[length(nt)], -1, 1) / sqrt(n)
  # x contient nt[length(nt)] simulations de la loi uniforme
  cumx <- cumsum(x)
  sigma <- 1/sqrt(3)
  rho <- 1/4
  return(data.frame(echantillon = cumx[nt] / sigma, borne_sup = (3*rho) / ((sigma^3) * sqrt(nt))))
}

b_unif <- brown_unif(n, c(1:t_max), moment3stud4, moment3norm)

# ------------- Calcul des bornes du graphique des trajectoires -------------

hist_min <- min(min(bt), min(b_norm$echantillon), min(b_stud3$echantillon), min(b_stud4$echantillon), min(b_rade$echantillon), min(b_unif$echantillon))
hist_max <- max(max(bt), max(b_norm$echantillon), max(b_stud3$echantillon), max(b_stud4$echantillon), max(b_rade$echantillon), max(b_unif$echantillon))

# ------------- Graphique des trajectoires -------------

hist_brown <- function(t_max, bt, b_norm, b_stud3, b_stud4, b_rade, b_unif, hist_min, hist_max){
  palette <- c("grey", "tomato3", "orange", "green", "blue", "purple")
  # création d'une palette de couleurs
  par(mfrow = c(1, 1))
  plot(1:t_max, bt, type = 'l', col = palette[1], xlab = 't', ylab = "X_t", 
       ylim = c(hist_min,hist_max), main = "Approximations du mouvement brownien")
  # mouvement brownien
  lines(1:t_max, rep(0,t_max), lty = 3)
  # droite horizontale y=0
  lines(1:t_max, b_norm, type = 'l', col = palette[2])
  # à l'aide de la N(0,1)
  lines(1:t_max, b_stud3, type = 'l', col = palette[3])
  # à l'aide de la Student(3)
  lines(1:t_max, b_stud4, type = 'l', col = palette[4])
  # à l'aide de la Student(4)
  lines(1:t_max, b_rade, type = 'l', col = palette[5])
  # à l'aide de la Rademacher
  lines(1:t_max, b_unif, type = 'l', col = palette[6])
  # à l'aide de la loi uniforme
  
  legend("bottomleft", c("Brownien", "N(0,1)", "Student(3)", "Student(4)", 
                         "Rademacher", "Uniforme"), col = palette, bg = "grey97", 
                          lwd = 1.5, box.lty = 0)
}

  

hist_brown(t_max, bt, b_norm$echantillon, b_stud3$echantillon, b_stud4$echantillon, b_rade$echantillon, b_unif$echantillon, hist_min, hist_max)

affiche_M_brownien <- function(n, M, t_max){
  plot(1:t_max, bt, type = 'l', col = "grey", xlab = 't', ylab = "Xt", ylim = c(-3*sqrt(t_max),3*sqrt(t_max)) , main = "Mvt Brownien") 
  
  for(i in 1:(M-1)){
    bt<-brown(t_max)
    lines(1:t_max, bt, type='l', col="grey")
  }
  
  lines(1:t_max , sqrt(1:t_max), col="black", lwd=2)
  lines(1:t_max , -sqrt(1:t_max), col="black", lwd=2)
  legend("topleft", c("Mvt Brownien", "t->sqrt(t)"), col = c("grey", "black"), bg = "grey97", lwd = 1.5, box.lty = 0)
}

affiche_M_approx <- function(n, M, t_max, loi, moment3stud4, moment3norm, nom_loi, couleur){
  echantillon<-loi(n,c(1:t_max), moment3stud4, moment3norm)
  plot(1:t_max, array(echantillon$echantillon), type = 'l', col = couleur, xlab = 't', ylab = "Xt", ylim = c(-3*sqrt(t_max),3*sqrt(t_max)), main = nom_loi) 
  
  for(i in 1:(M-1)){
    echantillon<-loi(n,c(1:t_max), moment3stud4, moment3norm)
    lines(1:t_max, echantillon$echantillon, type='l', col = couleur)
  }
  
  lines(1:t_max , sqrt(1:t_max), col="black", lwd=2)
  lines(1:t_max , -sqrt(1:t_max), col="black", lwd=2)
  legend("topleft", c(nom_loi, "t->sqrt(t)"), col = c(couleur, "black"), bg = "grey97", lwd = 1.5, box.lty = 0)
}


par(mfrow = c(2, 3))

affiche_M_brownien(n,100,t_max)
affiche_M_approx(n,100,t_max, brown_norm, moment3stud4, moment3norm, "Loi Normale", "tomato3")
affiche_M_approx(n,100,t_max, brown_stud3, moment3stud4, moment3norm, "Student 3", "orange")
affiche_M_approx(n,100,t_max, brown_stud4, moment3stud4, moment3norm, "Student 4", "green")
affiche_M_approx(n,100,t_max, brown_rade, moment3stud4, moment3norm, "Rademacher", "blue")
affiche_M_approx(n,100,t_max, brown_unif, moment3stud4, moment3norm, "Loi Uniforme", "purple")




# ------------- Fonction de densité et de répartition de la variable A_{t,n}, à t fixé -------------

# On choisit un certain réel t qui est fixé

# La fonction estim_fct_empirique trace 2 graphiques : 
# 1 - la fonction de densité empirique de A_{t,n} ;
# 2 - la f.d.r. empirique de A_{t,n} et la f.d.r. d'une N(0,t) encadrée par la borne de Berry-Esseen
# Elle renvoie l'écart maximal trouvé entre la f.d.r. empirique et la f.d.r. de la loi N(0,t)

estim_fct_empirique <- function(brown_simu, m, n, t, moment3stud4, moment3norm){
  vble <- numeric(length=m)
  for (i in 1:m){
    vble[i] <- brown_simu(n,t, moment3stud4, moment3norm)$echantillon
    # vble[i] contient une simulation de A_{t,n}
  }
  par(mfrow = c(1, 2))
  hist(vble, breaks=40, freq=FALSE, xlim=c(-10,10), ylim=c(0,0.2), main = "Densité empirique de A_{10,n}", xlab = 'x', ylab = "Densité")
  # Cet histogramme resprésente la densité empirique de A_{t,n}
  box()
  abs <- seq(-15,15,0.01)
  lines(abs, dnorm(abs,0,sqrt(t)))
  # On trace sur l'histogramme la fonction de densité d'une loi N(0,t)
  fdr_emp <- numeric(length=length(abs))
  ecart <- numeric(length=length(abs))
  j<-1
  for (i in abs){
    fdr_emp[j] <- sum(vble<=i)/m
    ecart[j] <- abs(fdr_emp[j] - pnorm(i,0,sqrt(t)))
    j<-j+1
  }
  borne_sup_emp <- max(ecart)
  borne_berry_essen <- brown_simu(n,t, moment3stud4, moment3norm)$borne_sup
  plot(abs,fdr_emp, type="l",col="red", main = "Fdr empirique de A_{10,n}", xlab = 'x', ylab = 'F')
  lines(abs,pnorm(abs,0,sqrt(t)), col="green")
  lines(abs, pnorm(abs,0,sqrt(t)) + brown_simu(n,t, moment3stud4, moment3norm)$borne_sup, col='blue')
  lines(abs, pnorm(abs,0,sqrt(t)) - brown_simu(n,t, moment3stud4, moment3norm)$borne_sup, col='blue')
  legend("topleft", c("Fdr_emp", "Fdr N(0,t)", "Borne B-E"), col = c("red", "green", "blue"), bg = "grey97", 
         lwd = 1.5, box.lty = 0)
  return(data.frame(sup_empirique = borne_sup_emp, sup_berry_essen = borne_berry_essen))
}


estim_fct_empirique(brown_stud4,100,100,10, moment3stud4, moment3norm)
estim_fct_empirique(brown_stud4,100,1000,10, moment3stud4, moment3norm)
estim_fct_empirique(brown_stud4,100,10000,10, moment3stud4, moment3norm)
estim_fct_empirique(brown_stud4,100,100,10, moment3stud4, moment3norm)
estim_fct_empirique(brown_stud4,1000,100,10, moment3stud4, moment3norm)
estim_fct_empirique(brown_stud4,10000,100,10, moment3stud4, moment3norm)



# ------------- Fonction de répartition empirique de la variable A_{t,n}, à t et n fixés -------------

# On modifie la fonction estim_fct_empirique 

estim_fct_empirique2 <- function(brown_simu, m, n, t, moment3stud4, moment3norm,titre){
  vble <- numeric(length=m)
  for (i in 1:m){
    vble[i] <- brown_simu(n,t, moment3stud4, moment3norm)$echantillon
    # vble[i] contient une simulation de A_{t,n}
  }
  #box()
  abs <- seq(-15,15,0.01)
  fdr_emp <- numeric(length=length(abs))
  ecart <- numeric(length=length(abs))
  j<-1
  for (i in abs){
    fdr_emp[j] <- sum(vble<=i)/m
    ecart[j] <- abs(fdr_emp[j] - pnorm(i,0,sqrt(t)))
    j<-j+1
  }
  borne_sup_emp <- max(ecart)
  plot(abs,fdr_emp, type="l",col="red", main = titre, xlab = 'x', ylab = 'F')
  lines(abs,pnorm(abs,0,sqrt(t)), col="green")
  lines(abs, pnorm(abs,0,sqrt(t)) + brown_simu(n,t, moment3stud4, moment3norm)$borne_sup, col='blue')
  lines(abs, pnorm(abs,0,sqrt(t)) - brown_simu(n,t, moment3stud4, moment3norm)$borne_sup, col='blue')
  #legend("topleft", c("Fdr_emp", "Fdr N(0,t)"), col = c("red", "green"), bg = "grey97", 
  #lwd = 1.5, box.lty = 0)
  return(borne_sup_emp)
}




affiche_tableau_fct_empirique <- function(m1, m2, m3, n, t, moment3stud4, moment3norm){
  ecart_emp <- matrix(rep(0, 12), nrow = 4)
  borne_berry_esseen <- rep(0, 4)
  
  par(mfrow = c(2,3))
  
  # student(4)
  borne_berry_esseen[1] <- brown_stud4(n,t, moment3stud4, moment3norm)$borne_sup
  ecart_emp[1,1] <- estim_fct_empirique2(brown_stud4, m1, n, t, moment3stud4, moment3norm, paste0("Student(4) ; m=",toString(m1)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[1,1], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[1], digits = 8)))), bg = "grey97", box.lty = 0)
  ecart_emp[1,2] <- estim_fct_empirique2(brown_stud4, m2, n, t, moment3stud4, moment3norm, paste0("Student(4) ; m=",toString(m2)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[1,2], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[1], digits = 8)))), bg = "grey97", box.lty = 0)
  ecart_emp[1,3] <- estim_fct_empirique2(brown_stud4, m3, n, t, moment3stud4, moment3norm, paste0("Student(4) ; m=",toString(m3)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[1,3], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[1], digits = 8)))), bg = "grey97", box.lty = 0)
  
  # normale standard
  borne_berry_esseen[2] <- brown_norm(n,t, moment3stud4, moment3norm)$borne_sup
  ecart_emp[2,1] <- estim_fct_empirique2(brown_norm, m1, n, t, moment3stud4, moment3norm, paste0("N(0,1) ; m=",toString(m1)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[2,1], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[2], digits = 8)))), bg = "grey97", box.lty = 0)
  ecart_emp[2,2] <- estim_fct_empirique2(brown_norm, m2, n, t, moment3stud4, moment3norm, paste0("N(0,1) ; m=",toString(m2)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[2,2], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[2], digits = 8)))), bg = "grey97", box.lty = 0)
  ecart_emp[2,3] <- estim_fct_empirique2(brown_norm, m3, n, t, moment3stud4, moment3norm, paste0("N(0,1) ; m=",toString(m3)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[2,3], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[2], digits = 8)))), bg = "grey97", box.lty = 0)
  
  par(mfrow = c(2,3))
  
  # uniforme([-1,1])
  borne_berry_esseen[4] <- brown_unif(n,t, moment3stud4, moment3norm)$borne_sup
  ecart_emp[4,1] <- estim_fct_empirique2(brown_unif, m1, n, t, moment3stud4, moment3norm, paste0("Uniforme ; m=",toString(m1)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[4,1], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[4], digits = 8)))), bg = "grey97", box.lty = 0)
  ecart_emp[4,2] <- estim_fct_empirique2(brown_unif, m2, n, t, moment3stud4, moment3norm, paste0("Uniforme ; m=",toString(m2)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[4,2], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[4], digits = 8)))), bg = "grey97", box.lty = 0)
  ecart_emp[4,3] <- estim_fct_empirique2(brown_unif, m3, n, t, moment3stud4, moment3norm, paste0("Uniforme ; m=",toString(m3)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[4,3], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[4], digits = 8)))), bg = "grey97", box.lty = 0)
  
  # Rademacher
  borne_berry_esseen[3] <- brown_rade(n,t, moment3stud4, moment3norm)$borne_sup
  ecart_emp[3,1] <- estim_fct_empirique2(brown_rade, m1, n, t, moment3stud4, moment3norm, paste0("Rademacher ; m=",toString(m1)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[3,1], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[3], digits = 8)))), bg = "grey97", box.lty = 0)
  ecart_emp[3,2] <- estim_fct_empirique2(brown_rade, m2, n, t, moment3stud4, moment3norm, paste0("Rademacher ; m=",toString(m2)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[3,2], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[3], digits = 8)))), bg = "grey97", box.lty = 0)
  ecart_emp[3,3] <- estim_fct_empirique2(brown_rade, m3, n, t, moment3stud4, moment3norm, paste0("Rademacher ; m=",toString(m3)))
  legend("topleft", c(paste0("sup_emp = ", toString(round(ecart_emp[3,3], digits = 8))), paste0("borne_theo = ", toString(round(borne_berry_esseen[3], digits = 8)))), bg = "grey97", box.lty = 0)
}



affiche_tableau_fct_empirique(m1 = 100, m2 = 1000, m3 = 10000, n = 100, t = 10, moment3stud4, moment3norm)





# On crée un processus qui ressemble à A_n pour des U([-1,1]) mais avec des pas de temps aléatoires, 
# qui sont iid de loi exponentielle(n).

pas_de_temps <- function(n,t){
  tau_k <- rexp(2*n*t, n) 
  # En général, on utilise environ n*t variables pour tracer A_n sur l'intervalle [0,t].
  # On tire donc 2*n*t variables exp(n) pour être sûr d'en avoir suffisamment.
  T_n <- cumsum(tau_k)
  while(T_n[length(T_n)] < t){
    tau_k <- cumsum(rexp(n*t, n))
    T_n <- append(T_n, T_n[length(T_n)] + tau_k)
    # Cette boucle n'est pas utilisé en général, 
    # elle sert seulement si 2*n*t variables exp(n) ne suffisent pas à atteindre t.
  }
  nb_pas <- length(which(T_n <= t))
  return(data.frame(temps_saut = T_n[1:nb_pas], nb_pas = nb_pas))
}

compound_process <- function(temps_saut, nb_pas, n, t){
  U_k <- runif(nb_pas,-1,1) * (sqrt(3/n))
  A_hat <- numeric(length=length(t))
  for (i in t){
    comptage <- sum(temps_saut <= i)
    A_hat[i] <- sum(U_k[1:comptage])
  }
  return(A_hat)
}

hist_unif_pas_alea <- function(n, t_max, b_unif){
  pas_tps_unif <- pas_de_temps(n, t_max)
  b_unif_pas_alea <- compound_process(pas_tps_unif$temps_saut, pas_tps_unif$nb_pas, n, 1:t_max)
  
  hist_min2 <- min(min(b_unif$echantillon), min(b_unif_pas_alea))
  hist_max2 <- max(max(b_unif$echantillon), max(b_unif_pas_alea))
  
  par(mfrow = c(1, 1))
  plot(1:t_max, b_unif_pas_alea, type = 'l', col = "red", main = "Approximations en utilisant des lois uniformes", xlab = 't', ylab = 'X_t', ylim = c(hist_min2,hist_max2))
  lines(1:t_max, b_unif$echantillon, type = 'l', col = "purple")
  lines(1:t_max, rep(0,t_max), lty = 3)
  # droite horizontale y=0
  legend("topleft", c("Pas de temps fixe", "Pas de temps aléatoire"), 
         col = c("purple", "red"), bg = "grey97", lwd = 1.5, box.lty = 0)
}

hist_unif_pas_alea(n, t_max, b_unif)
