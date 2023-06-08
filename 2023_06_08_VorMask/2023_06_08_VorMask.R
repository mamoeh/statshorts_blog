# ----- On the reversibility of Voronoi geomasking -----
# M. Möhler, June 2023


library(deldir)     # Voronoi tessellation
library(spatstat)   # NN-computation etc.
library(sdcSpatial) # example data
library(ggplot2)    # plotting
library(dplyr)      # convenience


# ----- Toy example -----

# create random point pattern
set.seed(20230528)
pts <- spatstat.random::rpoispp(win = owin(c(0, 1), c(0, 1)), lambda = 20) %>%
  as.data.frame() %>%
  mutate(id = 1:n())

pol <- deldir(pts) # compute Delauny triangulation 

# Voronoi masking is implemented in the package 'sdcSpatial', but the latter is
# still in an early version. Should the function later be recalled, 
# here's a substitute:
#   mask_voronoi <- function(x) {
#     x_nn <- x[nnwhich(x), ]
#     cbind((x[, 1] + x_nn[, 1]) / 2, (x[, 2] + x_nn[, 2]) / 2)
#   }

pts_m <- pts
pts_m[, c("x", "y")] <- sdcSpatial::mask_voronoi(pts[, c("x", "y")])


# visualize toy example

loc <- ggplot() + 
  geom_segment(data = pol$dirsgs, 
               aes(x1, y1, xend = x2, yend = y2),
               color = "#619CFF", lwd = .8) +
  geom_segment(data = pol$delsgs, 
               aes(x1, y1, xend = x2, yend = y2),
               lty = "dotted", color = "grey50", lwd = .8) +
  theme_void()

# before masking
loc + geom_point(data = pts, aes(x, y), pch = 15, cex = 3)

# after masking
loc + 
  geom_point(data = pts, aes(x, y), pch = 15, cex = 3, alpha = 0.3) +
  geom_point(data = pts_m, aes(x, y), pch = 15, cex = 3) +
  geom_segment(data = cbind(pts, rename(pts_m[, c("x", "y")], x_new = x, y_new = y)), 
               aes(x, y, xend = x_new, yend = y_new),
               color = "red", lwd = .8, 
               arrow = arrow(type = "closed", length = unit(0.15, "cm")))


# ----- Reversal experiment 1 -----

## (a) prepare example data

data("dwellings")
dwellings$ID <- 1:nrow(dwellings)

# transform to (0,0) origin for fewer-digit coordinates
dwellings$x <- dwellings$x - min(dwellings$x)
dwellings$y <- dwellings$y - min(dwellings$y)
# apply small jitter do avoid accidental duplicates after rounding
set.seed(12345)
dwellings$x <- round(dwellings$x + runif(nrow(dwellings), -0.1, 0.1), 3)
dwellings$y <- round(dwellings$y + runif(nrow(dwellings), -0.1, 0.1), 3)


## (b) anonymize

dwellings_m <- dwellings %>% rename(ID_orig = ID) # de-identify
dwellings_m[, c("x", "y")] <- sdcSpatial::mask_voronoi(dwellings[, c("x", "y")])

# measure displacement distances
disp_dist <- sqrt((dwellings$x - dwellings_m$x)^2 + (dwellings$y - dwellings_m$y)^2)
summary(disp_dist)
# measure orig. nearest neighbor distances
dist_nn <- nndist(dwellings[, c("x", "y")])
summary(dist_nn)
# verify: displacement distances are exactly 1/2 of nearest neighbor distances
summary(disp_dist) / summary(dist_nn)

# shuffle around masked locations
dwellings_m <- dwellings_m[sample(1:nrow(dwellings_m), nrow(dwellings_m), replace = FALSE), ]
rownames(dwellings_m) <- 1:nrow(dwellings_m) # new IDs of masked points

# how many singleton coordinates vs. twin coordinates?
dwellings_m$duplicate <- duplicated(dwellings_m[, c("x", "y")]) |
  duplicated(dwellings_m[, c("x", "y")], fromLast = TRUE)

table(dwellings_m$duplicate)                 # absolute
table(dwellings_m$duplicate)/nrow(dwellings) # relative

# verify: no. of dupl. locations after VM == no. of mutual nearest neighbors before
dwellings$nn <- nnwhich(dwellings[, c("x", "y")])
dwellings$duplicate <- dwellings[dwellings$nn, ]$nn == dwellings$ID
sum(dwellings$duplicate)


## (c) de-anonymize

dwellings_m$ID_found <- NA

# transform to ppp
win <- owin(c(min(dwellings$x), max(dwellings$x)), 
            c(min(dwellings$y), max(dwellings$y)))
ppp_dwell <- ppp(dwellings$x, dwellings$y, 
                 marks = dwellings$duplicate, 
                 window = win)
ppp_dwell_m <- ppp(dwellings_m$x, dwellings_m$y, 
                   marks = dwellings_m$duplicate, 
                   window = win)

## (c1) duplicates

# find 2 nearest neighbors from shortened candidate sets
dupl <- nncross(ppp_dwell_m[ppp_dwell_m$marks],
                ppp_dwell[ppp_dwell$marks],
                k = 1:2)

# For duplicates we just randomly assign the two nearest neighbors
# and are done

dupl$true_loc <- NA

# sort pairs to have lower-ID nearest neighbor first, even if found second
sc <- dupl$which.1 < dupl$which.2
sorted <- cbind(sc * dupl$which.1 + (!sc) * dupl$which.2, 
                sc * dupl$which.2 + (!sc) * dupl$which.1)
dupl[, c("which.1", "which.2")] <- sorted

# 1st member of each duplicate gets 1 of 2 locations at random ...
member <- which(!duplicated(dupl[, c("which.1", "which.2")]))

set.seed(12345)
rand_loc <- sample(1:2, size = nrow(dupl) / 2, replace = TRUE)

dupl[member, ]$true_loc <- ifelse(rand_loc == 1, 
                                  dupl[member, ]$which.1, 
                                  dupl[member, ]$which.2)

# ... 2nd member gets the other
dupl[-member, ]$true_loc <- ifelse(rand_loc == 1,
                                   dupl[-member, ]$which.2,
                                   dupl[-member, ]$which.1)

# search IDs of found locations in candidate set and assign as found 
dwellings_m[dwellings_m$duplicate, "ID_found"] <-
  dwellings[dwellings$duplicate, "ID"][dupl$true_loc]


## (c2) singletons

sing <- nncross(ppp_dwell_m[!ppp_dwell_m$marks], 
                ppp_dwell[!ppp_dwell$marks],
                k = 1:2)

sing$true_loc <- NA

# find true locations
sing$equi <- round(sing$dist.1, 4) == round(sing$dist.2, 4)
sing$dis1 <- sing$dist.1 < sing$dist.2
sing$dis1[sing$equi] <- NA
sing$true_loc <- sing$dis1 * sing$which.1 + !sing$dis1 * sing$which.2

# fill in the rest by successive elimination
while(any(is.na(sing$true_loc))) {
  
  # update already assigned locations
  tl <- sing$true_loc[!is.na(sing$true_loc)]
  
  # update to be assigned locations 
  sing$which.1[!is.na(sing$true_loc) | sing$which.1 %in% tl] <- 0
  sing$which.2[!is.na(sing$true_loc) | sing$which.2 %in% tl] <- 0
  
  # find one or more fits
  found <- which(xor(sing$which.1 == 0, sing$which.2 == 0))
  
  # assign
  sing[found, ]$true_loc <- sing[found, ]$which.1 + sing[found, ]$which.2
}

# search IDs of found locations in candidate set and assign as found 
dwellings_m[!dwellings_m$duplicate, "ID_found"] <- 
  dwellings[!dwellings$duplicate, "ID"][sing$true_loc]

## (c3) inspect success

dwellings_m$identified <- dwellings_m$ID_found == dwellings_m$ID_orig
dwellings_m$identified <- factor(dwellings_m$identified, levels = c(TRUE, FALSE))

sum(!dwellings$duplicate) + sum(dwellings$duplicate) / 2 # expected identifications
table(dwellings_m$identified)                            # actual identifications
table(dwellings_m$identified) / nrow(dwellings_m)        # share among all

# by duplicate
(succ <- table(dwellings_m$duplicate[dwellings_m$identified == TRUE])) # abs.
succ / table(dwellings_m$duplicate)                                    # rel.

# plot success
ggplot(filter(dwellings_m, x > 5e+3 & x < 6e+3 & y > 5e+3, y < 6e+3), 
       aes(x, y, color = identified)) +
  geom_point(cex = 0.01) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  theme_void()


# ----- Reversal experiment 2 -----

ssize <- 5000 # set sample size

# draw sample
set.seed(20230530)
dw_samp <- dwellings[sample(1:nrow(dwellings), size = ssize, replace = FALSE), 1:4]
dw_samp$ID <- 1:nrow(dw_samp)

# anonymize
dw_samp_m <- dw_samp
dw_samp_m[, c("x", "y")] <- sdcSpatial::mask_voronoi(dw_samp_m[, c("x", "y")])

# displacement distances (curiosity only)
disp_dist2 <- sqrt((dw_samp$x - dw_samp_m$x)^2 + (dw_samp$y - dw_samp_m$y)^2)
summary(disp_dist2)

# shuffle
dw_samp_m <- dw_samp_m[sample(1:nrow(dw_samp_m), nrow(dw_samp_m), replace = FALSE), ]
dw_samp_m$ID <- 1:nrow(dw_samp)

# number of expected identifications
(nloc <- nrow(unique(dw_samp_m[, c("x", "y")])))
nloc - (ssize - nloc) # ... of which singletons
2 * (ssize - nloc)          # ... of which duplicates


# ----- fun: reverse VM ------------------------------------- #
# X: masked coordinates (n by 2 matrix or data.frame)
# Y: coordinates of candidates (N by 2 matrix or data.frame)

reverse_voronoi <- function(X, Y) {
  
  lX <- nrow(X)
  lY <- nrow(Y)
  
  if(lX > lY) {stop("X must be shorter or equal in length to Y")}
  
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  Y$id <- 1:lY
  
  # compute duplicates
  X$dup <- duplicated(X) | duplicated(X, fromLast = TRUE)
  
  if(lX < lY) {
    
    ## Reconstruction step
    
    # upper bound for search distance
    nn1 <- nndist(X[, 1:2])
    nn2 <- nndist(X[, 1:2], k = 2)
    X$srch <- ifelse(X$dup, nn2, nn1)
    
    print("Reconstructing sample ...")
    
    # Find for each coordinate pair in the masked set two points in the candidate
    # set that are exactly equidistant from it. They are the true original location
    # and the original location's nearest neighbor from the sample.
    j <- rep(0, lX)
    for (i in 1:lX) {
      
      d <- round(sqrt((X[i, 1] - Y[, 1])^2 + (X[i, 2] - Y[, 2])^2), 6)
      
      jd <- which((duplicated(d) | duplicated(d, fromLast = TRUE)) & d < X$srch[i])
      
      # In the extremely unlikely case that we find two pairs of exactly equi-
      # distant locations among candidates, we choose the nearer ones.
      if(length(jd) > 2) {
        print(paste("Found", length(jd) / 2, "pairs for location", i))
        print(paste("Distances:", paste(d[jd], collapse = " , ")))
        print("Choosing closest pair.")
        jd <- jd[order(d[jd])]
      }
      
      if(!jd[1] %in% j) {j[j == 0][1] <- jd[1]}
      if(!jd[2] %in% j) {j[j == 0][1] <- jd[2]}
    }
    
    Y <- Y[j, ]
  }
  
  ## De-anomymizing step
  
  print("De-anonymizing ...")
  
  # compute (to-be) duplicates for candidate set
  nnY <- nnwhich(Y[, 1:2])
  Y$dup <- nnY[nnY] == 1:nrow(Y)
  
  w <- owin(c(min(Y[, 1]), max(Y[, 1])),
            c(min(Y[, 2]), max(Y[, 2])))
  
  suppressWarnings(
    p_X <- ppp(X[, 1], X[, 2], window = w, marks = X$dup))
    p_Y <- ppp(Y[, 1], Y[, 2], window = w, marks = Y$dup) 
  
  # de-anonymize (ca. half of) duplicates
    
  dupl <- nncross(p_X[p_X$marks], p_Y[p_Y$marks], k = 1:2)
  
  dupl$fnd <- NA
  
  sc <- dupl$which.1 < dupl$which.2
  srtd <- cbind(sc * dupl$which.1 + (!sc) * dupl$which.2,
                sc * dupl$which.2 + (!sc) * dupl$which.1)
  dupl[, c("which.1", "which.2")] <- srtd
  
  mbr <- which(!duplicated(dupl[, c("which.1", "which.2")]))
  rdm <- sample(c(TRUE, FALSE), size = nrow(dupl) / 2, replace = TRUE)
  
  dupl$fnd[mbr]  <- ifelse(rdm, dupl$which.1[mbr],  dupl$which.2[mbr])
  dupl$fnd[-mbr] <- ifelse(rdm, dupl$which.2[-mbr], dupl$which.1[-mbr])
  
  # de-anonymize singletons
  
  sing <- nncross(p_X[!p_X$marks], p_Y[!p_Y$marks], k = 1:2)
  
  sing$fnd <- NA
  
  eq <- round(sing$dist.1, 4) == round(sing$dist.2, 4)
  d1 <- sing$dist.1 < sing$dist.2
  d1[eq] <- NA
  sing$fnd <- d1 * sing$which.1 + !d1 * sing$which.2
  
  while(any(is.na(sing$fnd))) {
    
    fnd <- sing$fnd[!is.na(sing$fnd)]
    
    sing$which.1[!is.na(sing$fnd) | sing$which.1 %in% fnd] <- 0
    sing$which.2[!is.na(sing$fnd) | sing$which.2 %in% fnd] <- 0
    
    ft <- which(xor(sing$which.1 == 0, sing$which.2 == 0))
    sing$fnd[ft] <- sing$which.1[ft] + sing$which.2[ft]
  }
  
  # return de-identified locations
  
  X$fnd <- NA
  
  X$fnd[X$dup]  <- Y$id[Y$dup][dupl$fnd]
  X$fnd[!X$dup] <- Y$id[!Y$dup][sing$fnd]
  
  X[, c("x", "y", "fnd")]
}


# reversal with sample reconstruction
set.seed(1234)
dw_samp_fnd <- reverse_voronoi(dw_samp_m[, c("x", "y")], dwellings[, c("x", "y")])


# check success
dw_samp_fnd$identified <- dw_samp_fnd$fnd == rownames(dw_samp_fnd)

dwellings$sampled <- rownames(dwellings) %in% rownames(dw_samp)
dwellings$identified <- FALSE
dwellings$identified[dwellings$sampled] <- 
  dw_samp_fnd$identified[order(rownames(dw_samp_fnd))]

sum(dwellings$identified) / ssize # rate of de-anonymization

#plot success
dwellings$status <- ifelse(!dwellings$sampled, "not sampled", "sampled, not identified")
dwellings$status <- ifelse(dwellings$identified, "sampled, identified", dwellings$status)
dwellings$status <- as.factor(dwellings$status)

ggplot(filter(dwellings, x > 5e+3 & x < 7e+3 & y > 5e+3, y < 7e+3), 
       aes(x, y, color = status)) +
  geom_point(cex = 0.01) +
  scale_color_manual(values = c("lightgrey", "#F8766D", "#00BFC4")) +
  theme_void()

