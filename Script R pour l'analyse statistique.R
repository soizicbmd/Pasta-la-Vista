#ANALYSE DE LA SIGNIFICANCE DE NOS EXPÉRIENCES

temoin<-rnorm(8,mean=35, sd=1) 

# cq(DIL10) : 26.84

hist(temoin)

MuT<-mean(temoin)
sdT <-sd(temoin)

# Calcul de la probabilité d'avoir un Cq de 26,01 sans avoir fait le protocol d'enrichissement
dnorm(26.84, mean=MuT, sd=sdT)
#2.188393e-08

#Calcul de la probabilité d'avoir un Cq inférieur ou égale à 26,01 sans notre protocol
pnorm(26,84, mean=MuT, sd=sdT)
#1.221007e-10





# ETUDE DU NOMBRE DE COPIES DE READ


#Visualisation des données
donnees <- c(1, 7, 1, 8, 3, 2, 1, 1, 1, 6, 1, 4, 5, 3, 2, 3, 1, 1, 1, 1, 1, 5, 4, 2, 1, 2, 2, 1, 2, 4)
hist(donnees, breaks = seq(min(donnees) - 0.5, max(donnees) + 0.5, 1), col = "lightblue", main = "Depht of each homologous sequence")
lambda_estime <- mean(donnees)
points(0:max(donnees), dpois(0:max(donnees), lambda = lambda_estime) * length(donnees), type = "l", col = "red")
legend("topright", legend = "Poisson(lambda_estimed)", col = "red", pch = c(16, 19), lty = c(0, 0), bty = "n")

#La loi de Poisson semble bien correspondre à notre jeu de données

# Test de la loi de Poisson avec la méthode des chi-carré
observed <- donnees
expected <- dpois(0:29, lambda = mean(observed)) * sum(observed)
chisq.test(observed, p = expected / sum(expected))

#data:  observed
# X-squared = 3.2882e+19, df = 29, p-value < 2.2e-16

## Test de la loi de Poisson avec la méthode des QQ plot
qqplot(observed, y = qpois(ppoints(length(observed)), lambda = mean(observed)), main = "QQ plot test")

# Les points s'alignent le long d'une droite

# On teste la probabilité d'avoir 1 copie de gène 
dpois(1, lambda_estime, log = FALSE)
# 0.1970971


# On teste la probabilité d'avoir 2 copies de gène 
dpois(2, lambda_estime, log = FALSE)
#0.2529413

# On teste la probabilité d'avoir 3 copies de gène 
dpois(3, lambda_estime, log = FALSE)
#0.2164053

# On teste la probabilité d'avoir 4 copies de gène 
dpois(4, lambda_estime, log = FALSE)
#0.1388601

# lambda_estime = 2.566667
# Il semblerait donc que l'on aurait en moyenne 2,5 individus.

#21 profondeur de read sur le ribosome de PG

# On calcule la probabilité de vraisemblance

X <- c(1 : 15)
Y <- c(log10(dpois(21, lambda_estime, log = FALSE)),
log10(dpois(21, 2*lambda_estime, log = FALSE)),
log10(dpois(21, 3*lambda_estime, log = FALSE)),
log10(dpois(21, 4*lambda_estime, log = FALSE)),
log10(dpois(21, 5*lambda_estime, log = FALSE)),
log10(dpois(21, 6*lambda_estime, log = FALSE)),
log10(dpois(21, 7*lambda_estime, log = FALSE)),
log10(dpois(21, 8*lambda_estime, log = FALSE)),
log10(dpois(21, 9*lambda_estime, log = FALSE)),
log10(dpois(21, 10*lambda_estime, log = FALSE)),
log10(dpois(21, 11*lambda_estime, log = FALSE)),
log10(dpois(21, 12*lambda_estime, log = FALSE)),
log10(dpois(21, 13*lambda_estime, log = FALSE)),
log10(dpois(21, 14*lambda_estime, log = FALSE)),
log10(dpois(21, 15*lambda_estime, log = FALSE)))

maxY <- max(Y)


# Tracer les points
plot(X, Y, type = "o", col = "blue", pch = 16, main = "Probability of resemblance", xlab = "Number of copies", ylab = "Log(Probability of getting a depth of 21 reads knowing lambda_estimated)")

# Relier les points avec des lignes
lines(X, Y, type = "o", col = "blue")

#Maximum de vraissemblance
points(8, maxY, col = "red", pch = 19)
legend("bottomright", legend = "Maximum of resemblance", col = "red", pch = c(16, 19), lty = c(0, 0), bty = "n")


#Le maximum de ressemblane se trouve 
