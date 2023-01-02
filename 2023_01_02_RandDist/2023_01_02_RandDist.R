# ----- The average distance between points randomly placed in diagonally adjacent unit squares -----
# M. Moehler, Jan. 2023


library(ggplot2)   # for plotting
library(gganimate) # for animated plots


# ----- (1) Monte Carlo approximation of mean distance -----

n <- 1e+6 # length of Monte Carlo sequence

# drawing the required no. of random points
set.seed(400)
pts <- data.frame(x1 = runif(n, 0, 1), y1 = runif(n, 0, 1),
                  x2 = runif(n, 1, 2), y2 = runif(n, 1, 2))

pts$dist <- sqrt((pts$x1 - pts$x2)^2 + (pts$y1 - pts$y2)^2) # euclidean distances
mean(pts$dist) # Monte-Carlo approximation of mean distance


# ----- (2) Plots / Animations -----

# for the animation, a pre-computed sequence is used
load("ptsMC.R")
n <- nrow(ptsMC)

MC <- data.frame(x = c(ptsMC$x1, ptsMC$x2), y = c(ptsMC$y1, ptsMC$y2),
                 n = rep(1:n, 2))
anibckgr <- data.frame(x = rep(c(0, 1, 1, 2, 2, 0, 0), n),
                       y = rep(c(0, 0, 2, 2, 1, 1, 0), n),
                       n = rep(1:n, each = 7))

aniMC <- ggplot(MC, aes(x, y, group = n)) +
  geom_path(data = anibckgr) +
  geom_line() +
  theme_void() +
  transition_time(n/10) +
  shadow_mark(alpha = 0.3)

aniMC <- animate(aniMC, start_pause = 5, end_pause = 5, duration = 7)
aniMC

bckgr <- data.frame(n = 1:n, interc = 1.47356)

aniPath <- ggplot(ptsMC, aes(n, avgdist)) +
  geom_line(color = "blue") +
  geom_point() +
  geom_hline(data = bckgr, aes(yintercept = interc)) +
  theme_minimal() +
  ylab("mean dist.") +
  ylim(c(1.45, 1.525)) +
  transition_reveal(n)
aniPath <- animate(aniPath, start_pause = 5, end_pause = 5, duration = 7)
aniPath

