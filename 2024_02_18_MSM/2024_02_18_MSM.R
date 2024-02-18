# ----- Method of Simulated Moments -----
# M. Möhler, February 2024
# original version: M. Möhler, P. Bertazzoni, V. Kazakova (Feb. 2021)

# Model with endogeneity: 
#  M. Carrasco - A Regularization Approach to the Many Instruments Problem, 
#  Journal of Econometrics (2012)
# Code snippet for GMM:
#  P. Chauss - Computing Generalized Method of Moments and Generalized Empirical 
#  Likelihood with R, Journal of Statistical Software (2010)


library(mvtnorm)  # creation of the synthetic population
library(gmm)      # GMM estimator
library(ggplot2)  # plotting


# ----- GMM Estimator -------------------------------------------------------- #
# y:   dependent variable
# w:   independent variable
# x:   variable to be used for instruments

EstGMM <- function(y, w, x){
  
  h <- cbind(x, x^2, x^3)   # set of instruments for moment conditions
  g <- y ~ w                # model object as input to the GMM function
  
  # GMM call with instrumental variables
  res <- gmm::gmm(g, x=h, type="iterative", crit=1e-10, itermax=1000) 
  
  coef(res)[2]
}


# ----- MSM Estimator -------------------------------------------------------- #
# y:    dependent variable
# w:    independent variable
# R:    number of draws to use in the simulator
# anti: use antithetic sampling (yes/no)

ESTMSM <- function(y, w, R, anti = FALSE){
  
  m_y <- mean(y) # empirical moment
  
  # objective function:
  objective <- function(parameter, R) {
    
    if(anti) {
      Rand1 <- rnorm(R/2) 
      Rand2 <- rnorm(R/2)
      Rand1 <- c(Rand1, -Rand1)
      Rand2 <- c(Rand1, -Rand1)
    } else {
      Rand1 <- rnorm(R)
      Rand2 <- rnorm(R)
    }
    
    w_sim <- exp(-rnorm(R)^2) + Rand1
    y_sim <- parameter * w_sim + Rand2
    
    m_y_sim <- mean(y_sim)
    
    score <- (m_y_sim - m_y) %*% (m_y_sim - m_y)
    score
  }
  
  # optimization step
  res <- optimize(objective, interval = c(-1, 2), R = R) 
  
  res$minimum
}


# ----- 1. Generating synthetic population -----

Sig <- matrix(c(1, 0.5, 0.5, 1),2,2) # covariance matrix
N   <- 4000                          # population size
X   <- rnorm(N)                      # x is standard normally distributed
Err <- rmvnorm(N, sigma=Sig)         # error term (includes both u and epsilon)

delta <- 0.1                         # parameter to be estimated

W   <- exp(-X^2) + Err[,1]
Y   <- delta * W + Err[,2]

Pop <- data.frame(cbind(Y, W, X))


# ----- 2. Simulation -----

set.seed(20240217)

n <- 400       # sample size
SimRep <- 1000 # number of Simulation runs
types <- c("OLS", "GMM", "MSM1", "MSM2", "MSM3")

ResMSM1 <- ResMSM2 <- ResMSM3 <- ResGMM  <- ResOLS  <- numeric(SimRep)

for(i in 1:SimRep){
  
  # Drawing sample (SRS)
  s <- sample(1:N, size=n, replace=FALSE)
  samp <- Pop[s, ]
  
  y <- samp[,1]
  w <- samp[,2]
  x <- samp[,3]
  
  # Estimating:
  ResMSM1[i] <- ESTMSM(y, w, R = 100,   anti = FALSE)
  ResMSM2[i] <- ESTMSM(y, w, R = 10000, anti = FALSE)
  ResMSM3[i] <- ESTMSM(y, w, R = 100,   anti = TRUE)
  ResGMM[i]  <- EstGMM(y, w, x)
  ResOLS[i]  <- coef(lm(y ~ w))[2]
  
  cat("iteration: ", i, "\n")
}

SimResults <- data.frame(run      = rep(1:SimRep, length(types)),
                         estimate = c(ResOLS, ResGMM, ResMSM1, ResMSM2, ResMSM3),
                         type     = factor(rep(types, each = SimRep), levels = types)) 

SimMeans <- data.frame(type     = factor(types, levels = types),
                       estimate = round(c(mean(ResOLS), 
                                          mean(ResGMM), 
                                          mean(ResMSM1),
                                          mean(ResMSM2), 
                                          mean(ResMSM3)), 4),
                       sd       = round(c(sd(ResOLS), 
                                          sd(ResGMM), 
                                          sd(ResMSM1), 
                                          sd(ResMSM2), 
                                          sd(ResMSM3)), 4),
                       mse      = round(c(mean((ResOLS  - delta)^2),
                                          mean((ResGMM  - delta)^2),
                                          mean((ResMSM1 - delta)^2),
                                          mean((ResMSM2 - delta)^2),
                                          mean((ResMSM3 - delta)^2)), 4))


# ----- 3. Visualising Results -----

# colors following Paul Tol ('bright' palette)
# https://cran.r-project.org/web/packages/khroma/vignettes/tol.html
cols <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE")

# graph without MSM3, MSM4
set.seed(100)
ggplot(SimResults[SimResults$type != "MSM3", ], 
       aes(x = type, y = estimate, fill = type)) +
  geom_jitter(aes(color = type), alpha = 0.3, height = 0, show.legend = FALSE) +
  geom_boxplot(width = 0.2, show.legend = FALSE, outlier.color = NA) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  geom_point(data = SimMeans[SimMeans$type != "MSM3", ], 
             cex = 1.5, pch = 21, color = "grey20", show.legend = FALSE) +
  geom_hline(yintercept = delta, lty = "dashed") +
  theme_classic() +
  xlab(NULL)

# graph with MSM3
set.seed(100)
ggplot(SimResults, 
       aes(x = type, y = estimate, fill = type)) +
  geom_jitter(aes(color = type), alpha = 0.3, height = 0, show.legend = FALSE) +
  geom_boxplot(width = 0.2, show.legend = FALSE, outlier.color = NA) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  geom_point(data = SimMeans, 
             cex = 1.5, pch = 21, color = "grey20", show.legend = FALSE) +
  geom_hline(yintercept = delta, lty = "dashed") +
  theme_classic() +
  xlab(NULL)


# ----- 4. Comparing computation times -----

runs <- c(seq(1e+3, 9e+3, 1e+3), seq(1e+4, 9e+4, 1e+4), seq(1e+5, 9e+5, 1e+5), 1e+6)

set.seed(2563)

s <- sample(1:N, size=n, replace=FALSE)
samp <- Pop[s, ]
y <- samp[,1]
w <- samp[,2]
x <- samp[,3]

Var_Est1 <- Var_Est2 <- Time_Est1 <- Time_Est2 <- numeric(length(runs))

for(j in 1:length(runs)){
  
  r <- runs[j] # varying the strength of the simulator
  Est1 <- Est2 <- numeric(100)
  
  # Normal MSM Estimator:
  Time_Est1[j] <- system.time({
    for(i in 1:100){
      Est1[i] <- ESTMSM(y, w, r)
    }
    Var_Est1[j] <- var(Est1)
  })[3]
  
  # MSM Estimator with antithetic variables: 
  Time_Est2[j] <- system.time({
    for(i in 1:100){
      Est2[i] <- ESTMSM(y, w, r)
    }
    Var_Est2[j] <- var(Est2)
  })[3]
  
  print(paste("Done:", j, "To Do:", length(runs) - j))
}

CompResults <- data.frame(R        = rep(runs, 2), 
                          variance = c(Var_Est1, Var_Est2),
                          time     = c(Time_Est1, Time_Est2),
                          type     = rep(c("RV_iid", "RV_antit"), each = length(runs)))
CompResults$variance <- CompResults$variance / CompResults$variance[1]
CompResults$time <- CompResults$time / CompResults$time[1]
  
# displaying the variance-computation-trade-off
ggplot(CompResults[CompResults$type == "RV_iid", ], aes(x = R)) +
  geom_smooth(aes(y = variance), lty = "dashed") +
  geom_smooth(aes(y = time / 700), color = "red") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(sec.axis = sec_axis(trans = ~. * 700, name="time")) +
  theme_classic()

