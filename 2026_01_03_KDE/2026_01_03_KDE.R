##### Support recovery for KDE #################################################
# this version:     2026-01-03
# original version: 2025-03-28
############################################################################## #

library(ks)
library(glmnet)
library(plotmo)
library(ggplot2)


# ----- Helper function -----

make_K_mat <- function(x, eval.points, H, normalize = FALSE, thresh = NA) {
  
  # decide 1D or 2D case
  dimH <- dim(as.matrix(H))[1]
  
  # N x N matrix K_h from Hut (2020) (2.5)
  if(dimH == 1) {
    
    # 1D
    N <- length(x)
    G <- length(eval.points)
    K_mat <- matrix(data = 0, nrow = G, ncol = N)
    
    for(i in 1:N) {
      K_mat[, i] <- ks::kde(x[i], h = H, eval.points = eval.points)$estimate
    }
    
  } else {
    
    # 2D
    N <- nrow(x)
    G <- nrow(eval.points)
    K_mat <- matrix(data = 0, nrow = G, ncol = N)
    
    for(i in 1:N) {
      
      K_mat[, i] <- ks::kde(x[i, ], H = H, eval.points = eval.points)$estimate
    }
  }
  
  # C_h from Hut (2020) (2.6)
  if(normalize) {
    K_mat <- K_mat / rowSums(K_mat)
  }
  # option to set all entries below some threshold to 0
  if(!is.na(thresh)) {
    K_mat[K_mat < thresh] <- 0
  }
  
  K_mat
}


# ----- Minimal example: Hut method (2D) -----

locs2d <- cbind(c(-0.55, -0.30, -0.35, 0.00, 0.05, 0.60, 0.70),
                c(-0.50, -0.40,  0.00, 0.00, 0.45, -0.15, -0.25))

locs2d_w <- c(5, 30, 10, 20, 20, 5, 30)

pop2d <- rbind(locs2d, c(0.00, -0.60), c(-0.45, 0.45), c(0.10, 0.35),
               c(-0.65, -0.35), c(-0.25, -0.25), c(0.30, 0.15))

grid_size <- 0.05
ev_pts2d <- as.matrix(expand.grid(round(seq(-1, 1, grid_size), 2),
                                  round(seq(-1, 1, grid_size), 2)))

# diagonal bandwidth matrix
H <- ks::Hpi.diag(locs2d)
# diagonal symmetric bandwidth matrix
H <- diag(rep(mean(diag(H)), 2))
# non-diagonal bandwidth matrix
#H <- ks::Hpi(locs2d)

y2d   <- ks::kde(x = locs2d, H = H, eval.points = ev_pts2d)
y2d_w <- ks::kde(x = locs2d, H = H, eval.points = ev_pts2d, w = (locs2d_w / sum(locs2d_w)) * length(locs2d_w))

kde2d   <- data.frame(x = ev_pts2d[, 1], y = ev_pts2d[, 2], dens.est. = y2d$estimate)
kde2d_w <- data.frame(x = ev_pts2d[, 1], y = ev_pts2d[, 2], dens.est. = y2d_w$estimate)

ggplot(kde2d, aes(x, y)) +
  geom_raster(aes(fill = dens.est.), alpha = 0.5) +
  scale_fill_viridis_c() +
  geom_point(data = data.frame(x = locs2d[, 1], y = locs2d[, 2]),
             pch = 21, fill = "red", color = "white", size = 2.5) +
  geom_point(data = data.frame(x = pop2d[8:13, 1], y = pop2d[8:13, 2]),
             pch = 21, fill = "black", color = "white", size = 2.5) +
  xlab("E") + ylab("N") +
  theme_bw() +
  theme(legend.position = "bottom", legend.key.height = unit(0.3, "cm")) +
  ggtitle("true")

ggplot(kde2d_w, aes(x, y)) +
  geom_raster(aes(fill = dens.est.)) +
  scale_fill_viridis_c() +
  geom_point(data = data.frame(x = locs2d[, 1], y = locs2d[, 2], w = locs2d_w),
             color = "red", aes(size = w)) +
  theme_bw()


# ----- (1) Density weighted, membership known -----

y2d_w$estimate <- y2d_w$estimate * sum(locs2d_w) # in the case of KDE: to KDS

places <- apply(locs2d, 1, function(x) {which(ev_pts2d[, 1] == x[1] & ev_pts2d[, 2] == x[2])})
m_vec <- y2d_w$estimate[places]

K_h <- make_K_mat(x = as.data.frame(locs2d), eval.points = locs2d, H = H)
res <- m_vec %*% solve(K_h)

rbind(res, locs2d_w)


# ----- (2) Density unweighted, membership sensitive -----

places <- apply(pop2d, 1, function(x) {which(ev_pts2d[, 1] == x[1] & ev_pts2d[, 2] == x[2])})
m_vec <- y2d$estimate[places]

C_h <- make_K_mat(x = as.data.frame(pop2d), eval.points = pop2d, H = H, normalize = TRUE)

res <- m_vec %*% solve(C_h)
pick <- which(zapsmall(res) > 0)
cbind(pop2d[pick, ], locs2d)

res_pop2d <- data.frame(x = pop2d[, 1], y = pop2d[, 2],
                        included = c(rep(TRUE, 7), rep(FALSE, 6)),
                        reconstr = round(as.numeric(res), 2))
res_pop2d$inferred <- FALSE
res_pop2d$inferred[pick] <- TRUE

ggplot(kde2d, aes(x, y)) +
  geom_raster(aes(fill = dens.est.), alpha = 0.5) +
  scale_fill_viridis_c() +
  geom_text(data = res_pop2d, aes(label = reconstr, color = inferred), show.legend = FALSE,
            size = 4, fontface = "bold") +
  scale_color_manual(values = c("black", "red")) +
  xlab("E") + ylab("N") +
  theme_bw() +
  theme(legend.position = "bottom", legend.key.height = unit(0.3, "cm")) +
  ggtitle("inferred")


# ----- (3) via Lasso shrinkage -----

# apply lasso
res_las <- glmnet(x = C_h, y = m_vec, intercept = FALSE, alpha = 1, 
                  lower.limits = 0, upper.limits = 1)
plot(res_las, xvar = "lambda", label = TRUE) # coefficient profiles
plotmo::plot_glmnet(res_las)

cand_sizes <- apply(res_las$beta, 2, function(x) sum(x != 0))
res_las$beta[, min(which(cand_sizes == nrow(locs2d)))] # estimate for true active set size
res_las$beta[, ncol(res_las$beta)]                     # estimate with undershrinkage

# prepare for plotting
las_df <- as.data.frame(expand.grid(loc = 1:nrow(pop2d), 
                                    step = 1:length(res_las$lambda), 
                                    lambda = 0, w_i = 0, nloc = 0, shrinkage = 0, 
                                    samp = FALSE))
for(i in 1:length(res_las$lambda)) {
  las_df$lambda[las_df$step == i] <- res_las$lambda[i]
  las_df$w_i[las_df$step == i] <- res_las$beta[, i]
  las_df$nloc[las_df$step == i] <- sum(abs(res_las$beta[, i]) > 0)
}
las_df$samp[las_df$loc %in% 1:7] <- TRUE # true locations
las_df$shrinkage <- log10(las_df$lambda) # plot by standardized shrinkage factor
las_df$shrinkage <- (las_df$shrinkage - min(las_df$shrinkage)) / (max(las_df$shrinkage) - min(las_df$shrinkage))

# make colors and labels for shrinkage profiles
loc_cols <- heat.colors(14)[1:7]
pop_cols <- grey.colors(14, start = 0, end = 0.3)[1:7]

loc_pnts <- data.frame(shrinkage = min(las_df$shrinkage[las_df$shrinkage >= 0.375]),
                       w_i = las_df$w_i[las_df$shrinkage == min(las_df$shrinkage[las_df$shrinkage >= 0.375])],
                       loc = 1:13)
loc_labs <- data.frame(shrinkage = 0.35, 
                       w_i = las_df$w_i[las_df$shrinkage == min(las_df$shrinkage[las_df$shrinkage >= 0.375])],
                       loc = 1:13)
# tweak labels for plotting
loc_labs$w_i[4]  <- loc_labs$w_i[4] + 0.005
loc_labs$w_i[10] <- loc_labs$w_i[10] - 0.005
loc_labs$shrinkage[c(8, 9, 11, 13)] <- loc_labs$shrinkage[c(8, 9, 11, 13)] - c(0, 0.03, 0, 0.05)
loc_labs$w_i[c(8, 9, 11, 13)] <- loc_labs$w_i[c(8, 9, 11, 13)] + c(0.015, 0.015, - 0.015, -0.015)

# plot shrinkage profiles
ggplot(las_df[las_df$shrinkage >= 0.375, ], aes(shrinkage, w_i, group = loc, color = as.factor(loc))) +
  geom_line(show.legend = FALSE) +
  geom_text(data = loc_labs, aes(label = loc), show.legend = FALSE, size = 3) +
  geom_point(data = loc_pnts, show.legend = FALSE) +
  #geom_vline(xintercept = min(las_df$shrinkage[las_df$nloc == 7]), lty = "dashed", color = "blue") +
  scale_color_manual(values = c(loc_cols, pop_cols)) +
  scale_x_reverse(limits = c(1.0, 0.3)) +
  theme_minimal() +
  xlab("penalty strength") +
  ggtitle("variable selection (Lasso)")

# select either the true-size active set or the undershrunk version
pick_las <- which(res_las$beta[, ncol(res_las$beta)] > 0)
las_sel <- data.frame(x = pop2d[pick_las, 1], y = pop2d[pick_las, 2])

# plot selection
ggplot(kde2d, aes(x, y)) +
  geom_raster(aes(fill = dens.est.), alpha = 0.5) +
  scale_fill_viridis_c() +
  geom_point(data = data.frame(x = locs2d[, 1], y = locs2d[, 2]),
             pch = 21, fill = "red", color = "white", size = 2.5) +
  geom_point(data = data.frame(x = pop2d[8:13, 1], y = pop2d[8:13, 2]),
             pch = 21, fill = "black", color = "white", size = 2.5) +
  #geom_point(data = las_sel, pch = 1, size = 3, color = "blue") +
  geom_point(data = las_sel, size = 2.5, color = "blue") +
  geom_text(data = data.frame(x = locs2d[, 1], y = locs2d[, 2], lbl = 1:7),
            aes(label = lbl), color = "red", size = 3.5, 
            nudge_x = -0.10, nudge_y = -0.04) +
  geom_text(data = data.frame(x = pop2d[8:13, 1], y = pop2d[8:13, 2], lbl = 8:13),
            aes(label = lbl), color = "black", size = 3.5, 
            nudge_x = -0.10, nudge_y = -0.04) +
  xlab("E") + ylab("N") +
  theme_bw() +
  theme(legend.position = "bottom", legend.key.height = unit(0.3, "cm")) +
  ggtitle("inferred (Lasso)")

## assess fit of solution

y2d_t <- ks::kde(pop2d[pick_las, ], H = H, eval.points = ev_pts2d)

kde2d$dens.cand. <- y2d_t$estimate
kde2d$dens.diff. <- abs(kde2d$dens.cand. - kde2d$dens.est.)
(mise_curr <- mean((y2d_t$estimate - y2d$estimate)^2)) # non-optimal MISE

ggplot(kde2d, aes(x, y)) +
  geom_raster(aes(fill = dens.diff.), alpha = 0.5) +
  geom_point(data = data.frame(x = pop2d[pick_las, 1], y = pop2d[pick_las, 2]),
             pch = 21, fill = "black", color = "white", size = 2.5) +
  scale_fill_viridis_c(option = "plasma") +
  xlab("E") + ylab("N") +
  theme_bw() +
  theme(legend.position = "bottom", legend.key.height = unit(0.3, "cm"))

# We select too much with the Lasso
# --> do stepwise improvement (backward elimination)

n_el <- length(pick_las) - nrow(locs2d) # how many to eliminate?
bw_sel <- as.data.frame(expand.grid(i         = c(pick_las),
                                    iteration = 1:n_el,
                                    mise      = NA))

pick_curr <- pick_las # initialize with solution that's too large
for(j in 1:n_el) {
  
  mise_i <- vector("numeric", length(pick_curr))
  for(i in 1:length(pick_curr)) {
    
    y2d_i <- ks::kde(pop2d[pick_curr[-i], ], H = H, eval.points = ev_pts2d)
    mise_i[i] <- mean((y2d_i$estimate - y2d$estimate)^2)
    bw_sel$mise[bw_sel$iteration == j & bw_sel$i == pick_curr[i]] <- mise_i[i]
  }
  
  # deselect location that leads to best mise improvement
  pick_curr <- pick_curr[-which.min(mise_i)] 
  
  print(paste("iteration:", j, "current MISE:", min(mise_i)))
}

pick_curr # refined solution vector

bw_sel$step <- as.factor(paste("iteration", bw_sel$iteration))
bw_sel$true_loc <- FALSE
bw_sel$true_loc[bw_sel$i %in% 1:7] <- TRUE
# data only for the deleted points (when and why)
bw_delete <- bw_sel[(bw_sel$iteration == 1 & bw_sel$i == 10) | (bw_sel$iteration == 2 & bw_sel$i == 12), ]

# plot updating process
ggplot(bw_sel, aes(as.factor(i), color = true_loc)) +
  geom_hline(data = bw_delete, aes(yintercept = mise), lty = "dashed", color = "blue", alpha = 0.6) +
  geom_segment(aes(xend = as.factor(i), y = 0, yend = mise), show.legend = FALSE) +
  geom_point(aes(y = mise), show.legend = FALSE) +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(~step, nrow = 2) +
  xlab("i") + ylab("MISE") +
  theme_minimal()

# plot refined solution
bw_las_sel <- data.frame(x = pop2d[pick_curr, 1], y = pop2d[pick_curr, 2])

ggplot(kde2d, aes(x, y)) +
  geom_raster(aes(fill = dens.est.), alpha = 0.5) +
  scale_fill_viridis_c() +
  geom_point(data = data.frame(x = locs2d[, 1], y = locs2d[, 2]),
             pch = 21, fill = "red", color = "white", size = 2.5) +
  geom_point(data = data.frame(x = pop2d[8:13, 1], y = pop2d[8:13, 2]),
             pch = 21, fill = "black", color = "white", size = 2.5) +
  geom_point(data = bw_las_sel, size = 2.5, color = "blue") +
  geom_text(data = data.frame(x = locs2d[, 1], y = locs2d[, 2], lbl = 1:7),
            aes(label = lbl), color = "red", size = 3.5, 
            nudge_x = -0.10, nudge_y = -0.04) +
  geom_text(data = data.frame(x = pop2d[8:13, 1], y = pop2d[8:13, 2], lbl = 8:13),
            aes(label = lbl), color = "black", size = 3.5, 
            nudge_x = -0.10, nudge_y = -0.04) +
  xlab("E") + ylab("N") +
  theme_bw() +
  theme(legend.position = "bottom", legend.key.height = unit(0.3, "cm")) +
  ggtitle("inferred (Lasso + Backward elim.)")

