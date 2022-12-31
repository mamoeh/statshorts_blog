# ----- Ridge regression in the logit model: an application to birth weight data -----
# M. Moehler, Dec. 2022
# Based on joint work with V. Kazakova in March 2022
# Data: G.A.F. Seber & C.J. Wild, Nonlinear regression, Wiley, 1989


library(glmnet)    # for ridge regression
library(ggplot2)   # for plotting
library(gganimate) # for animated path plot


# (1) Data preparation -----

seber <- read.csv("Seber.csv", sep = ";") # adjust to local setting
# BW : Birth weight
# BPD: Biparietal diameter
# AC : Abdominal circumference

seber$VLBW   <- ifelse(seber$BW < 1500, 1, 0) # creating binary target variable
seber$noVLBW <- 1 - seber$VLBW
seber$VLBW_f <- factor(seber$VLBW, levels = c(1, 0), labels = c("yes", "no"))
seber$rank   <- rank(seber$BW)


# birth weight plot
ggplot(seber, aes(rank, BW)) +
  geom_point(aes(color = VLBW_f, pch = VLBW_f), show.legend = FALSE) +
  geom_hline(yintercept = 1500, lty = "dashed") +
  theme_minimal()

# BPD vs. AC plot
ggplot(seber, aes(BPD, AC)) +
  geom_point(aes(color = VLBW_f, pch = VLBW_f)) +
  scale_color_discrete(name = "VLBW") +
  scale_shape_discrete(name = "VLBW") +
  theme_minimal()


# the two predictors are highly correlated
cor(seber$AC, seber$BPD, method = "pearson")
cor(seber$AC, seber$BPD, method = "spearman")

seber[, 2:4] <- scale(seber[, 2:4]) # scaling is required for ridge


# (2) Fitting models on original data -----

# IRLS
modGlm <- glm(cbind(VLBW, noVLBW) ~ BPD + AC, data = seber, family = "binomial")
coef(modGlm)

# Ridge
mmX <- model.matrix(BW ~ BPD + AC, data = seber) # model matrix (including intercept)
cv <- cv.glmnet(mmX, seber$VLBW, alpha = 0, family = "binomial") # Cross-validation
plot(cv)
(lam <- cv$lambda.min) # optimal lambda from CV
modRdg <- glmnet(mmX, seber$VLBW, alpha = 0, family = "binomial")
coef(modRdg, s = lam)


## Animated path plot

lams <- seq(0.01, 10, 0.01) # lambda values to show
resLam <- data.frame(lam   = rep(lams, each = 2),
                     coeff = rep(c("BPD", "AC"), length(lams)),
                     estim = numeric(length = length(lams)*2))
resLam$coeff <- factor(resLam$coeff, levels = unique(resLam$coeff))
# extracting solution paths for given lambda range
for(i in seq(lams)) { resLam[resLam$lam == lams[i], ]$estim <- coef(modRdg, s = lams[i])[3:4, ] }
resLam$frame <- log10(resLam$lam)

aniRdg <- ggplot(resLam, aes(lam, estim, color = coeff, group = coeff)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0) +
  xlab("lambda") +
  ylab("estimate") +
  theme_minimal() +
  scale_color_manual(values = c("red", "darkgreen"), name = "coefficient") +
  transition_reveal(frame)

aniRdg # view animation


# (3) Fitting models on bootstrap data -----

runs <- 100 # number of bootstrap samples
resGLM <- resRDG <- data.frame(Interc = numeric(length = runs), 
                               BPD    = numeric(length = runs),
                               AC     = numeric(length = runs),
                               pred   = numeric(length = runs)) # empty results data.frame
i    <- 1 # running index
nErr <- 0 # count convergence errors
set.seed(100)
while(i <= runs){
  
  sam <- sample(1:50, size = 50, replace = TRUE) # draw bootstrap sample
  samSeb <- seber[sam, ]                         # subset to chosen sample
  
  ## fitting models
  
  # glm
  err <- tryCatch(
    modGLM <- glm(cbind(VLBW, noVLBW) ~ BPD + AC, data = samSeb, family = "binomial"),
    warning = function(w) TRUE)
  if(class(err)[1] == "logical"){
    nErr <- nErr + 1
    next} # upon failure to converge: re-try with new sample
  
  # ridge
  mmX <- model.matrix(BW ~ BPD + AC, data = samSeb)
  cv <- cv.glmnet(mmX, samSeb$VLBW, nfolds = 10, alpha = 0, family = "binomial")
  lam <- cv$lambda.min
  modRDG <- glmnet(mmX, samSeb$VLBW, alpha = 0, lambda = lam, family = "binomial")
  
  resGLM[i, 1:3] <- modGLM$coefficients
  resRDG[i, 1:3] <- coef(modRDG)[-2]
  
  ## measuring predictive accuracy
  
  predGLM <- predict.glm(modGLM, newdata = seber, type = "response")
  predRDG <- predict.glmnet(modRDG, newx = model.matrix(BW ~ BPD + AC, data = seber))
  pred <- data.frame(GLM = predGLM, RDG = as.numeric(predRDG), truth = seber$VLBW)
  
  pred$GLM <- ifelse(pred$GLM > 0.5, 1, 0)
  pred$RDG <- ifelse(pred$RDG > 0, 1, 0)
  
  resGLM[i, 4] <- sum(pred$GLM == pred$truth) / 50
  resRDG[i, 4] <- sum(pred$RDG == pred$truth) / 50
  
  print(paste("Done:", i, "To Do:", runs - i))
  i <- i + 1
}

# collating results
res <- data.frame(estimate = c(resGLM$Interc, resGLM$BPD, resGLM$AC,
                               resRDG$Interc, resRDG$BPD, resRDG$AC),
                  coeff    = rep(rep(c("Interc.", "BDP", "AC"), each = runs), 2),
                  method   = rep(c("IRLS", "Ridge"), each = runs*3))
res$coeff  <- factor(res$coeff, levels = unique(res$coeff))
res$method <- as.factor(res$method)

ggplot(res[res$coeff != "Interc.", ], aes(coeff, estimate)) +
  geom_boxplot(aes(color = method)) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "darkorange")) +
  xlab("coefficient")

# assessing accuracy gain
resPred <- data.frame(accuracy = c(resGLM$pred, resRDG$pred),
                      method   = rep(c("IRLS", "Ridge"), each = runs))

ggplot(resPred, aes(accuracy, method, color = method)) +
  geom_boxplot(show.legend = FALSE) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "darkorange"))

