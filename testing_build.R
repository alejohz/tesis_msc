# Importing packages
library("tidyr")
library("ICSNP")
library("robustbase") #Qn and Sn
library("ggplot2") #Plotting
library("robets") #modified robets
#source("robets/R/robets.R")

set.seed(seed = 10042022) #seed es la fecha
e <- rnorm(500, mean = 0, sd = 1) # error

alpha <- 0.5
beta <- 0.003
theta <- 0.9
m_1 <- 0
for (i in 2:length(e)) {
  m_1[i] <- alpha + beta * i + theta * m_1[i - 1] + e[i]
}
x_y <- m_1
# Stack overflow  for time series simulation (Study simulated data by CREVTIS)

contaminate <- function(data_cont, i){
  set.seed(seed = 13032022)
  e <- rnorm(1, mean = 0, sd = 1)
  i <- 0.1
  alfa <- round(length(data_cont) * i)
  for (j in seq(1, alfa, 1)) {
      w <- round(runif(1, 1, length(data_cont)))
      if (w %% 2 == 0) data_cont[w] <- data_cont[w] + e * 2
      else data_cont[w] <- data_cont[w] - e * 2
  }
  return(data_cont)
  }
data_cont_0.1 <- contaminate(m_1, 0.1)
data_cont_0.4 <- contaminate(m_1, 0.4)

# Psi function comparison
# Qn as scale estimator
# Data wirh 0 % Contamination
a <- c("Huber", "Bisquare", "Hampel", "Welsh1", "Welsh2")
l <- c()
m <- c()
for (i in a) {
  m_ahenao <- robets(x_y, model = "AAN", scale.estimator = "Qn", psifun = i)
  l[[paste0(i)]] <- m_ahenao$robaicc
  m[[paste0(i)]] <- m_ahenao
}
# l <- cbind(a, l)
# l <- rbind(c("psifun", "lik", "aic", "bic", "aicc"), l)
print("SIMULATED DATA")
print(which.min(l))
# print(l)

# Psi function comparison
# Qn as scale estimator
# Data wirh 10 % Contamination
a <- c("Huber", "Bisquare", "Hampel", "Welsh1", "Welsh2")
l <- c()
for (i in a) {
  m_ahenao <- robets(data_cont_0.1, model = "AAN", scale.estimator = "Qn",
   psifun = i)
  l[[paste0(i)]] <- m_ahenao$robaicc
}
# l <- cbind(a, l)
# l <- rbind(c("psifun", "lik", "aic", "bic", "aicc"), l)
print("CONTAMINATED DATA WITH 0.1")
# print(l)
print(which.min(l))

a <- c("Huber", "Bisquare", "Hampel", "Welsh1", "Welsh2")
l <- c()
for (i in a) {
  m_ahenao <- robets(data_cont_0.4, model = "AAN", scale.estimator = "Qn",
   psifun = i)
  l[[paste0(i)]] <- m_ahenao$robaicc
}
# l <- cbind(a, l)
# l <- rbind(c("psifun", "lik", "aic", "bic", "aicc"), l)
print("CONTAMINATED DATA WITH 0.4")
# print(l)
print(which.min(l))

values <- c(27, 27, 7, 24, 39, 40, 24, 45, 36, 37, 31, 47, 16, 24, 6, 21,
35, 36, 21, 40, 32, 33, 27, 42, 14, 21, 5, 19, 31, 32, 19, 36,
29, 29, 24, 42, 15, 24, 21)
a <- c("Huber", "Bisquare", "Hampel", "Welsh1", "Welsh2")
l <- c()
for (i in a) {
  m_ahenao <- robets(values, model = "AAN", scale.estimator = "Qn",
   psifun = i)
  l[[paste0(i)]] <- m_ahenao$robaicc
}
# l <- cbind(a, l)
# l <- rbind(c("psifun", "lik", "aic", "bic", "aicc"), l)
print("Crevits simulated data (Outlier near end)")
print(which.min(l))

# print(l)
# m_a <- robets(X_Y, model = "AAN", scale.estimator = "Qn", psifun = "Huber")
# print(m_a)

# m_a <- robets(X_Y, model = "AAN", scale.estimator = "Qn", psifun = "Bisquare")
# print(m_a)
# m_cont <- robets(data_cont, model = "AAN", scale.estimator = "Qn",
#   psifun = "Bisquare")
# print(m_cont)