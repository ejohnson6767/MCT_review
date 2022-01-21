
# functions for simulating and determining initial conditions -------------

# population map of species 1, the gleaner
N1Dynam <- function(N1, R, pars){
  N1 * ((pars[["w"]] * pars[["mu1"]] * R)/(pars[["K1"]] + R) - pars[["d"]])
}

# population map of species 2, the opportunist
N2Dynam <- function(N2, R, pars){
  N2 * ((pars[["w"]] * pars[["mu2"]] * R)/(pars[["K2"]] + R) - pars[["d"]])
}

# map for the resource
RDynam <- function(N1, N2, R, pars){
  pars[["d"]]*(pars[["S"]] - R) - ((pars[["mu1"]] * R * N1)/(pars[["K1"]] + R)) - ((pars[["mu2"]] * R * N2)/(pars[["K2"]] + R))
}

# function for a single time-step
step <- function(N1, N2, R, pars){
  # local growth
  N1 <- N1 + N1Dynam(N1, R, pars)*pars[["dt"]];
  N2 <- N2 + N2Dynam(N2, R, pars)*pars[["dt"]];
  R <- R + RDynam(N1, N2, R, pars)*pars[["dt"]];
  
  # dispersal dynamics
  disp1 <- (1-pars[["q1"]])*N1*pars[["dt"]]
  disp2 <- (1-pars[["q2"]])*N2*pars[["dt"]]
  N1 <- N1 - disp1 + sum(disp1)/pars[["num_patches"]]
  N2 <- N2 - disp2 + sum(disp2)/pars[["num_patches"]]
  
  list(N1 = N1, N2 = N2, R = R)
}

# functions for simulating multiple time-steps
sim <- function(N1, N2, R, pars, num_steps){
  for(i in 1:num_steps)
  {
    temp <- step(N1, N2, R, pars)
    N1 <- temp$N1
    N2 <- temp$N2
    R <- temp$R
  }
  return(list(N1 = N1, N2 = N2, R = R))
}

# functions for simulating multiple time-steps and recording density/concentration data
sim_record <- function(N1_init, N2_init, R_init, pars, num_steps){
  
  N1_mat <- matrix(NA, nrow = pars[["num_patches"]], ncol = num_steps+1)
  N2_mat <- matrix(NA, nrow = pars[["num_patches"]], ncol = num_steps+1)
  R_mat <- matrix(NA, nrow = pars[["num_patches"]], ncol = num_steps+1)
  N1_mat[,1] <- N1_init
  N2_mat[,1] <- N2_init
  R_mat[,1] <- R_init
  
  for(i in 1:num_steps)
  {
    temp <- step(N1_mat[,i], N2_mat[,i], R_mat[,i], pars)
    N1_mat[,i+1] <- temp[["N1"]]
    N2_mat[,i+1] <- temp[["N2"]]
    R_mat[,i+1] <- temp[["R"]]
  }
  return(list(N1 = N1_mat, N2 = N2_mat, R = R_mat))
}

# equilibrium consumer density when there is only one resident
NEq <- function(w, d, mu, K, S0)
{
  (w * (d * (K + S0) - mu * S0 * w))/ (d - mu * w)
}

# equilibrium resource concentration when there is only one resident
REq <- function(w, d, mu, K, S0)
{
  (d * K)/(mu * w - d)
}

# proof of concept: coexistence -------------------------------------------

# specify parameters
num_patches <- 1000
S0 <- 3 # mean resource supply point
sigma <- 0.5
S <- rnorm(num_patches, S0, sigma)
S <- ifelse(S < 0,0,S)

pars <- list(
  num_patches = num_patches,
  dt = 0.01,
  d = 0.59,
  w = 1,
  S0 = S0,
  mu1 = 1,
  K1 = 0.4,
  q1 = 0.1,
  mu2 = 1000,
  K2 = 1000,
  q2 = 0.1,
  sigma = sigma,
  S = S
  )

# plot opporutnist-gleaner trade-off
xmax <- 1
x <- seq(0,xmax,0.01)
b1 <- pars[["w"]]*pars[["mu1"]]*x / (pars[["K1"]] + x)
b2 <- pars[["w"]]*pars[["mu2"]]*x / (pars[["K2"]] + x)
plot(x, b1, type = "l", col = "orange", lwd = 3, xlim = c(0,xmax))
lines(x, b2, type = "l", col = "green", lwd = 3)
abline(h=pars[["d"]], col = "Blue")

N1Eq <- NEq(pars[["w"]], pars[["d"]], pars[["mu1"]], pars[["K1"]], pars[["S0"]])
N2Eq <- NEq(pars[["w"]], pars[["d"]], pars[["mu2"]], pars[["K2"]], pars[["S0"]])
REq_inv2 <- REq(pars[["w"]], pars[["d"]], pars[["mu1"]], pars[["K1"]], pars[["S0"]])
REq_inv1 <- REq(pars[["w"]], pars[["d"]], pars[["mu2"]], pars[["K2"]], pars[["S0"]])


### INVASION ANALYSIS

# invasion analysis parameters
tmax <- 10000 
thin <- 100 # thin data archive for faster plotting

## species 1 is the invader
# intial conditions 
N1_init <- rep(0.01, pars[["num_patches"]])
N2_init <- rep(N2Eq, pars[["num_patches"]])
R_init <- rep(REq_inv1, pars[["num_patches"]])
# simulate
res <- sim_record(N1_init, N2_init, R_init, pars, tmax)
# plot invasion time series of species 1
plot(colMeans(res$N1[,seq(1,tmax+1, by = thin)]), type = "l", col = "orange", lwd = 3, main = "invasion time series of species 1", xlab = "time", ylab = "density", xaxt='n')

## species 2 is the invader
# intial conditions
N1_init <- rep(N1Eq, pars[["num_patches"]])
N2_init <- rep(0.01, pars[["num_patches"]])
R_init <- rep(REq_inv2, pars[["num_patches"]])
# simulate
res <- sim_record(N1_init, N2_init, R_init, pars, tmax)
plot(colMeans(res$N2[,seq(1,tmax+1, by = thin)]), type = "l", col = "green", lwd = 3, main = "invasion time series of species 2", xlab = "time", ylab = "density", xaxt='n')


### BOTH SPECIES AT TYPICAL ABUNDANCES

### intial conditions when both species are abundance
N1_init <- rep(N1Eq/2, pars[["num_patches"]])
N2_init <- rep(N1Eq/2, pars[["num_patches"]])
R_init <- rep(REq_inv1, pars[["num_patches"]])

tmax <- 100000
res <- sim_record(N1_init, N2_init, R_init, pars, tmax)
plot(colMeans(res$N1[,seq(1,tmax+1, by = thin)]), type = "l", col = "orange", lwd = 3, ylim = c(0.5,1.5), xlab = "Time", ylab = "Density", xaxt='n')
lines(colMeans(res$N2[,seq(1,tmax+1, by = thin)]), type = "l", col = "green", lwd = 3)
lines(colMeans(res$R[,seq(1,tmax+1, by = thin)]), type = "l", col = "purple", lwd = 3)
legend("topright", 
       legend = c("Consumer 1: gleaner", "Consumer 2: opportunist", "Resource"), 
       col = c("orange","green","purple"),
       lty= 1,
       lwd= 3)


# resource variation as a function of local retention --------------------

# specify parameters
num_patches <- 1000
S0 <- 3 # mean resource supply point
sigma <- 0.1
S <- rnorm(num_patches, S0, sigma)
S <- ifelse(S < 0,0,S)
q_vec <- seq(0,1, by = 0.05)
R_var <- numeric(length(q_vec))
tmax <- 1000

# simulations
for(i in seq_along(q_vec))
{
  pars <- list(
    num_patches = num_patches,
    dt = 0.01,
    d = 0.59,
    w = 1,
    S0 = S0,
    mu1 = 1,
    K1 = 0.4,
    q1 = q_vec[i],
    mu2 = 1000,
    K2 = 1000,
    q2 = 0.1,
    sigma = sigma,
    S = S
  )
  res <- sim(1, 0, 1, pars, tmax)
  R_var[i] <- var(res[["R"]])
}

# approximation
q_range <- seq(0,1,by = 0.001)
R_var_estim <- pars[["sigma"]]^2 *(((pars[["d"]] * pars[["K1"]] * (1 - q_range) * pars[["mu1"]] * pars[["w"]])/(
  pars[["d"]]^2 * (1 + (-1 + pars[["d"]]) * q_range) * (pars[["K1"]] + pars[["S0"]]) - 
    pars[["d"]] * pars[["mu1"]] * (pars[["d"]] * pars[["K1"]] * q_range + 2 * (1 + (-1 + pars[["d"]]) * q_range) * pars[["S0"]]) * pars[["w"]] + (1 + (-1 + pars[["d"]]) * q_range) * pars[["mu1"]]^2 * pars[["S0"]] * pars[["w"]]^2)))^2

# plot results
ymax <- max(c(R_var,R_var_estim))
plot(x = q_vec, y = R_var, pch = 16, ylim = c(0,ymax), xlab = "Local retention", ylab = "Resource variation")
lines(x = q_range, y = R_var_estim, col = "Red", pch = 16, lwd = 2)
legend("bottomleft", 
       legend = c("Simulation", "Approximation"), 
       col = c("Black","Red","purple"),
       lty=c(NA, 1),
       pch=c(16, NA),
       lwd= 3)
### interestingly, for reasons we are unsure of, the analytical approximation performs better when there is high local retention


# resource variation as a function of opportunism --------------------

# specify parameters
num_patches <- 1000
S0 <- 3 # mean resource supply point
sigma <- 0.1
S <- rnorm(num_patches, S0, sigma)
S <- ifelse(S < 0,0,S)

# fixed parameters
d = 0.59
w = 1
q1 = 0.9
R0 <- 0.6
tmax <- 1000

# varying parameters
K1_vec <- seq(0.1,1, by = 0.05)
mu1_vec <- d*(K1_vec + R0)/(R0*w)

R_var <- numeric(length(K1_vec))
R_mean <- numeric(length(K1_vec))

# simulations
for(i in seq_along(K1_vec))
{
  pars <- list(
    num_patches = num_patches,
    dt = 0.1,
    d = d,
    w = w,
    S0 = S0,
    mu1 = mu1_vec[i],
    K1 = K1_vec[i],
    q1 = q1,
    mu2 = 1000,
    K2 = 1000,
    q2 = 0.1,
    sigma = sigma,
    S = S
  )
  res <- sim(1, 0, 1, pars, tmax)
  R_var[i] <- var(res[["R"]])
  R_mean[i] <- mean(res[["R"]])
}

# approximation
K1_range <- seq(min(K1_vec),max(K1_vec),by = 0.001)

R_var_estim <- ((-1 + pars[["q1"]])^2 * R0^2 * (K1_range + 
                    R0)^2 * pars[["sigma"]]^2)/(pars[["d"]] * K1_range * pars[["q1"]] * (R0 - pars[["S0"]]) + (-1 + pars[["q1"]]) * (R0^2 + K1_range * pars[["S0"]]))^2

## plot results

# the mean of R should be consistently close to R0
plot(x = K1_vec, y = R_mean, pch = 16, xlab = expression("Gleaner" %<->% "Opportunist"), ylab = "Resource variation")

# plot variance
ymax <- max(c(R_var,R_var_estim))
plot(x = K1_vec, y = R_var, pch = 16, ylim = c(0,ymax), xlab = expression("Gleaner" %<->% "Opportunist"), ylab = "Resource variation")
lines(x = K1_range, y = R_var_estim, col = "Red", pch = 16, lwd = 2)
legend("bottomleft", 
       legend = c("Simulation", "Approximation"), 
       col = c("Black","Red","purple"),
       lty=c(NA, 1),
       pch=c(16, NA),
       lwd= 3)

