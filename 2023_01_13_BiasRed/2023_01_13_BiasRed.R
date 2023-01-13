# ----- Bias-reduced estimation from geodata masked for privacy -----
# M. Möhler, Jan. 2023


library(spatstat)   # for spatial data operations
library(ggplot2)    # for plotting
library(ggforce)    # - " -
library(gganimate)  # for animated example of random perturbation
library(sdcSpatial) # source for the 'dwellings' data set


# ----- Helper functions -----

# function to calculate random new coordinates for an observation 
# based on uniform displacement distance and random angle
#   x: vector of x-coordinates
#   y: vector of y-coordinates
#   rmax: maximum displacement distance
#   rmin: minimum displacement distance (0 by default)
relocate <- function(x, y, rmax, rmin = 0) {
  
  # draw random angle and distance
  ang <- runif(1, 0, 360)
  dis <- runif(1, rmin, rmax)

  # compute new coordinates
  xn <- x + cos((ang * pi) / 180) * dis
  yn <- y + sin((ang * pi) / 180) * dis
  
  return(c(xn, yn))
}


# function to calculate for a set of locations each one's distance to
# the closest of a set of points of interest
#   loc: n by 2 matrix of location coordinates
#   poi: m by 2 matrix of POI coordinates
get_mindist <- function(loc, poi) {
  
  # transform locations and POI to point patterns
  wind <- owin(c(min(c(loc[, 1], poi[, 1])), max(c(loc[, 1], poi[, 1]))),
               c(min(c(loc[, 2], poi[, 2])), max(c(loc[, 2], poi[, 2]))))
  pppLoc <- ppp(loc[, 1], loc[, 2], window = wind)
  pppPoi <- ppp(poi[, 1], poi[, 2], window = wind)
  
  # find and return distance to nearest POI for each location
  nnD <- nncross(pppLoc, pppPoi, what = "dist")
  return(nnD)
}

# function to get approximate location based on simulated points
#   X:     a (perturbed) coordinate pair
#   ppsim: simulated points (ppp)
#   dmax:  maximum displacement distance of previous perturbation
approximate_location <- function(X, ppsim, dmax) {
  
  X <- as.numeric(X) # a coordinate pair
  # compute distance to all simulated points
  sdist <- sqrt((ppsim$x - X[1])^2 + (ppsim$y - X[2])^2)
  
  # If no points are simulated within the possible range, the coordinates are
  # kept, otherwise we calculate a centroid
  npt <- sum(sdist <= dmax)
  xc <- ifelse(npt > 0, mean(ppsim$x[sdist <= dmax]), X[1])
  yc <- ifelse(npt > 0, mean(ppsim$y[sdist <= dmax]), X[2])
  
  return(cbind(xc, yc))
}


# ----- Animation -----

# create an original location pattern
set.seed(100)
pp <- as.data.frame(rMatClust(kappa = 10, scale = 0.1, mu = 15))
le <- nrow(pp) # no. of locations

# create points of interest
poi <- data.frame(x = c(0.2, 0.6, 0.7),
                  y = c(0.7, 0.2, 0.6))

# original location pattern with POI
ggplot(pp, aes(x, y)) +
  geom_point() +
  geom_point(data = poi, color = "blue", pch = 3, cex = 2) +
  theme_void()

# preparing frame-specific data
pts <- data.frame(x = numeric(le * (le + 1)),
                  y = numeric(le * (le + 1)),
                  type = "original", 
                  step = rep(1:(le + 1), each = le),
                  dist = 99)
pois <- data.frame(x = rep(poi$x, le + 1),
                   y = rep(poi$y, le + 1),
                   step = rep(1:(le + 1), each = nrow(poi)))
mdis <- data.frame(step = 1:(le + 1),
                   mdis = 99)

rmax <- 0.1 # set maximum displacement distance

# initialize with original locations
pts[pts$step == 1, c("x", "y")] <- as.matrix(pp)
pts[pts$step == 1, "dist"] <- get_mindist(pp[, 1:2], poi[, 1:2])
mdis[mdis$step == 1, "mdis"] <- mean(pts$dist[pts$step == 1])

for(i in 1:nrow(pp)) {
  
  # initialize step with locations from previous step
  pts[pts$step == i + 1, 1:3] <- pts[pts$step == i, 1:3]
  
  # change only one location per step (for animation)
  xy <- pts[pts$step == i + 1, c("x", "y")][i, ]
  
  # calculate new location
  xy_new <- as.numeric(relocate(xy[1], xy[2], rmax = rmax))
  
  pts[pts$step == i + 1, c("x", "y")][i, ] <- xy_new
  pts[pts$step == i + 1, "type"][i] <- "new"
  
  # calculate new distances
  pts[pts$step == i + 1, "dist"] <- get_mindist(cbind(pts$x[pts$step == i + 1],
                                                      pts$y[pts$step == i + 1]),
                                                poi[, 1:2])
  mdis[mdis$step == i + 1, "mdis"] <- mean(pts[pts$step == i + 1, "dist"])
  
  print(paste("Done:", i, "To Do:", nrow(pp) - i))
}
pts$type <- factor(pts$type, levels = unique(pts$type))

aniLoc <- ggplot(pts, aes(x, y, color = type)) +
  xlim(c(0 - rmax, 1 + rmax)) +
  ylim(c(0 - rmax, 1 + rmax)) +
  geom_point() +
  theme_void() +
  scale_color_manual(values = c("black", "red")) +
  geom_point(data = pois, color = "blue", pch = 3, cex = 2) +
  transition_states(step)

aniLoc <- animate(aniLoc, start_pause = 5, end_pause = 5, duration = 38)
aniLoc

aniMdis <- ggplot(mdis, aes(step, mdis)) +
  geom_line(color = "blue") +
  geom_point() + 
  xlab("# locations changed") +
  ylab("avg. dist. to nearest POI") +
  theme_minimal() +
  transition_reveal(step)

aniMdis <- animate(aniMdis, start_pause = 5, end_pause = 5, duration = 38)
aniMdis


# ----- Simulation -----

runs <- 1000 # no. of Monte Carlo runs
dmax <- 500  # maximum displacement distance
n <- 100     # sample size
s <- 1e+6    # no. of points to simulate for location approximation

data("dwellings") # data set to be used
dwellings <- dwellings[, c("x", "y")] # only coordinates are used

# prepare boundaries
xrng <- c(min(dwellings$x), max(dwellings$x))
yrng <- c(min(dwellings$y), max(dwellings$y))

# a minimal jitter (1m) is applied to make each location unique
set.seed(100)
suppressWarnings(
dwellings <- as.data.frame(rjitter(ppp(dwellings$x, dwellings$y, owin(xrng, yrng)), 
                                   radius = 1, retry = TRUE)))

# prepare quadrants
mx <- mean(dwellings$x)
my <- mean(dwellings$y)
dwellings$quadrant <- 1
dwellings$quadrant[dwellings$y < my] <- 4
dwellings$quadrant[dwellings$x < mx] <- 2
dwellings$quadrant[dwellings$x < mx & dwellings$y < my] <- 3

# prepare results data.frames
resEst <- data.frame(run   = rep(1:runs, each = 2*3),
                     model = rep(rep(c("true", "naive", "calibrated"), each = 2), runs),
                     coeff = rep(c("Intercept", "distance"), runs*3),
                     estim = 0)
resMSE <- data.frame(run   = rep(1:runs, each = 2),
                     model = rep(c("naive", "calibrated"), runs),
                     mse   = 0)

set.seed(200)
for(i in 1:runs) {
  
  # select at random 2 POIs per quadrant
  poiQ1 <- sample(which(dwellings$quadrant == 1), 2)
  poiQ2 <- sample(which(dwellings$quadrant == 2), 2)
  poiQ3 <- sample(which(dwellings$quadrant == 3), 2)
  poiQ4 <- sample(which(dwellings$quadrant == 4), 2)
  poiQ <- c(poiQ1, poiQ2, poiQ3, poiQ4)
  
  POI <- dwellings[poiQ, ]    # points of interest
  dwell <- dwellings[-poiQ, ] # remaining respondents
  dwell <- dwell[sample(1:nrow(dwell), n), ] # subsample
  
  # create density map
  den <- density.ppp(ppp(dwell$x, dwell$y, window = owin(xrng, yrng)), 
                     sigma = 150, positive = TRUE)
  
  # displace locations (geomasking)
  ang <- runif(n, 0, 360)
  r <- runif(n, 0, dmax)
  dwell$xn <- dwell$x + cos((ang * pi) / 180) * r
  dwell$yn <- dwell$y + sin((ang * pi) / 180) * r
  
  # simulate locations based on density map
  psim <- rpoint(s, f = den)
  
  # approximate locations
  aloc <- t(apply(cbind(dwell$xn, dwell$yn), 1, approximate_location,
                  ppsim = psim, dmax = dmax))
  dwell$xc <- aloc[, 1]
  dwell$yc <- aloc[, 2]
  
  # calculate distances to POI
  dwell$dist  <- get_mindist(dwell[, c("x", "y")],   POI[, c("x", "y")])
  dwell$distn <- get_mindist(dwell[, c("xn", "yn")], POI[, c("x", "y")])
  dwell$distc <- get_mindist(dwell[, c("xc", "yc")], POI[, c("x", "y")])
  
  # make target variable
  dwell$targ <- 50 + dwell$dist + rnorm(n, 0, 100)
  
  # estimating models
  modOrig <- lm(targ ~ dist, data = dwell)
  modNaiv <- lm(targ ~ distn, data = dwell)
  modCali <- lm(targ ~ distc, data = dwell)
  
  # MSE
  mseNaiv <- mean((dwell$dist - dwell$distn)^2)
  mseCali <- mean((dwell$dist - dwell$distc)^2)
  
  # storing results
  resEst[resEst$run == i & resEst$model == "true", "estim"] <- modOrig$coefficients
  resEst[resEst$run == i & resEst$model == "naive", "estim"] <- modNaiv$coefficients
  resEst[resEst$run == i & resEst$model == "calibrated", "estim"] <- modCali$coefficients
  
  resMSE[resMSE$run == i & resMSE$model == "naive", "mse"] <- mseNaiv
  resMSE[resMSE$run == i & resMSE$model == "calibrated", "mse"] <- mseCali
  
  print(paste("Done:", i, "To Do:", runs- i))
}


# formatting results
resEst$model <- factor(resEst$model, levels = unique(resEst$model)[c(2, 3, 1)])
resMSE$model <- factor(resMSE$model, levels = unique(resMSE$model))
resMSE$mse <- as.numeric(resMSE$mse)

# colors following Paul Tol ('high contrast' palette)
# https://cran.r-project.org/web/packages/khroma/vignettes/tol.html
cols <- c("#BB5566", "#DDAA33", "#004488")

ggplot(resEst[resEst$coeff == "distance", ], aes(x = model, y = estim, group = model)) +
  geom_violin(aes(fill = model), trim = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.1) +
  geom_hline(yintercept = 1, lty = "dashed") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  ylab("estimate") +
  xlab(NULL)

ggplot(resMSE, aes(model, mse, fill = model)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = cols[1:2]) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme_classic() +
  ylab("MSE") +
  xlab(NULL)


# ----- Explanatory graphic -----

set.seed(300)
denFull <- density.ppp(ppp(dwellings$x, dwellings$y, window = owin(xrng, yrng)),
                       sigma = 50, positive = TRUE)

# picking a small no. of simulated points to explain the procedure
psimFull <- as.data.frame(rpoint(3e+5, f = denFull))

# picking a part of the map to explaint the procedure on
xprng <- c(156500, 157500)
yprng <- c(462500, 463500)

psimFull <- psimFull[psimFull$x >= xprng[1] & psimFull$x <= xprng[2] &
                       psimFull$y >= yprng[1] & psimFull$y <= yprng[2], ]

obs  <- 717 # picking a point to explain the procedure with
maxr <- 200
omod <- c(cos((375 * pi) / 180) * 150, sin((375 * pi) / 180) * 150)

# original and displaced location
locs <- data.frame(x = c(psimFull$x[obs],
                         psimFull$x[obs] + omod[1]),
                   y = c(psimFull$y[obs],
                         psimFull$y[obs] + omod[2]))

psimFull <- psimFull[-obs, ]
psimFull$dobs <- sqrt((psimFull$x - locs$x[2])^2 + (psimFull$y - locs$y[2])^2)
psimFull$poss <- psimFull$dobs < maxr

# adding estimated location
locs <- rbind(locs, c(mean(psimFull$x[psimFull$poss]), mean(psimFull$y[psimFull$poss])))
arrs <- data.frame(x = locs$x[1:2], y = locs$y[1:2], xend = locs$x[2:3], yend = locs$y[2:3])

ggplot(dwellings, aes(x, y)) +
  geom_point(data = psimFull, aes(color = poss), alpha = .5, show.legend = FALSE) +
  scale_color_manual(values = c("gray80", "gray40")) +
  geom_density_2d(show.legend = FALSE, color = "gray30") +
  geom_segment(data = arrs, aes(xend = xend, yend = yend), lty = "dashed",
               color = c("red", "blue")) +
  geom_point(data = locs, color = c("black", "red", "blue"), pch = 3, cex = 5) +
  geom_circle(data = locs[2, ], aes(x0 = x, y0 = y, r = maxr)) +
  xlim(xprng - c(0, 250)) +
  ylim(yprng + c(250, 0)) +
  theme_void()

