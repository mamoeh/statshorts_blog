# ----- CART synthesis of small-scale georeferences -----
# M. Möhler, Jan. 2023

library(raster)   # for creating geographic grid cells
library(sf)       # for plotting geographic data
library(ggplot2)  # for general plotting
library(dplyr)    # convenience tools
library(spatstat) # for simulating point locations
library(synthpop) # contains the CART synthesizer


# ----- Loading AMELIA data -----

## loading variables (customize path to local setting)

pathAMELIA <- "2023_01_20_SynthGeo/AMELIA_HH_level_v0.2.3"
varsAMELIA <- paste(pathAMELIA, list.files(path = pathAMELIA), sep = "/")
for(i in varsAMELIA) { load(i) }

vnames <- c("CIT", "DIS", "DOU", "EDI", "EHHS", "HH050", "HHS", "HID",
            "HS010", "HS020", "HS030", "HS040", "HS050", "HS060", "HS070", "HS080",
            "HS100", "HS110", "HY010", "HY020", "HY025", "HY030", "HY040", "HY050",
            "HY060", "HY070", "HY080", "HY090", "HY100", "HY110", "HY120", "HY130",
            "HY140", "INC", "PROV", "REG", "SOC", "UEP")

# bind to single data.frame
for(i in seq(vnames)) {
  
  if(i == 1) {
    AMELIAv023 <- as.data.frame(get(vnames[i]))
  } else {
    AMELIAv023 <- cbind(AMELIAv023, get(vnames[i]))
  }
}
names(AMELIAv023) <- vnames
rm(list = vnames)
AMELIAv023 <- AMELIAv023[order(AMELIAv023$CIT), ] # order by city

# formatting
AMfactors <- c("PROV", "REG", "DIS", "UEP", "DOU")
AMELIAv023[, AMfactors] <- lapply(AMELIAv023[, AMfactors], as.factor)


## Loading the map (customize path to local setting)

load("2023_01_20_SynthGeo/AMAP_v0.2.1.RData")
AMAP <- st_as_sf(AMAP)

# add count data to the map
CITc <- as.data.frame(table(AMELIAv023$CIT))
names(CITc) <- c("CIT", "HH_n")
AMAP <- merge(AMAP, CITc, by= "CIT")
AMAP$HH_dens <- AMAP$HH_n / (AMAP$AREA / 1e+6)


# ----- Creating orig. georeferences -----

# prepare an empty data frame for positioning households
AMELIAxy <- data.frame(CIT = AMELIAv023$CIT,
                       x = 0, y = 0)

# simulate uniformly distributed point locations
set.seed(123)
for(i in seq(AMAP$CIT)) {
  
  city <- AMAP[AMAP$CIT == AMAP$CIT[i], ]
  # create simulation window from AMELIA map
  w <- as.owin(st_geometry(city))
  # simulate point data & assign to households
  rp <- as.data.frame(rpoint(city$HH_n, f = city$HH_dens, win = w))
  AMELIAxy[AMELIAxy$CIT == AMAP$CIT[i], c("x", "y")] <- rp
  
  print(paste("Cities done:", i, "To do:", nrow(AMAP) - i))
}
AMELIAv023 <- cbind(AMELIAv023, AMELIAxy[, c("x", "y")]) # attach coordinates
rm(AMELIAxy, CITc)

# coordinates to grid cell
cellwidth <- 1000
AMELIAv023$xc <- ceiling(AMELIAv023$x / cellwidth)
AMELIAv023$yc <- ceiling(AMELIAv023$y / cellwidth)
AMELIAv023$cell <- paste0(AMELIAv023$xc, ",", AMELIAv023$yc)


# ----- Synthesizing -----

# custom predictor matrix
pm <- matrix(ncol = ncol(AMELIAv023), nrow = ncol(AMELIAv023), data = 0)
rownames(pm) <- colnames(pm) <- colnames(AMELIAv023)
pm["cell", ] <- 1
pm["cell", c("cell", "x", "y", "HID", "PROV", "REG", "DIS", "CIT", "xc", "yc")] <- 0

# custom method set (use only CART)
me <- rep("", ncol(AMELIAv023))
me[which(colnames(AMELIAv023) == "cell")] <- "cart"

AMp2 <- AMELIAv023[AMELIAv023$PROV == 2, ] # choose only 2nd AMELIA-province
dstr <- unique(AMp2$DIS)
synSeed <- c(1111, 2222, 3333)

# synthesize per district
AMp2s <- AMp2
for(i in seq(dstr)) {
  
  # subset data to district
  AMd <- AMp2[AMp2$DIS == dstr[i], ]
  AMd$cell <- as.factor(AMd$cell)
  mf <- length(unique(AMd$cell))
  
  # synthesize grid cell via CART for select district
  AMdSyn <- syn(AMd, 
                visit.sequence = c("cell"),
                method = me,
                predictor.matrix = pm,
                seed = synSeed[i],
                maxfaclevels = mf,
                cart.minbucket = 10)
  
  # store synthetic cell
  AMp2s[AMp2s$DIS == dstr[i], "cell"] <- as.character(AMdSyn$syn$cell) 
  
  print(paste("Districts done:", i, "To Do:", length(dstr) - i))
}


# ----- Assessing synthesis -----

# building grid
crng_x <- min(AMp2$xc):max(AMp2$xc)
crng_y <- min(AMp2$yc):max(AMp2$yc)

grid <- raster(ncols = length(crng_x), nrows = length(crng_y),
               xmn = (crng_x[1] - 1) * cellwidth,
               ymn = (crng_y[1] - 1) * cellwidth,
               xmx = crng_x[length(crng_x)] * cellwidth,
               ymx = crng_y[length(crng_y)] * cellwidth)
grid <- st_as_sf(rasterToPolygons(grid))[, -1]

cells <- expand.grid(xc = crng_x, yc = rev(crng_y))
grid$cell <- paste0(cells$xc, ",", cells$yc)

# aggregating by grid cell
AMp2_sum <- AMp2 %>% group_by(cell) %>%
  summarise(HHS_sum = sum(HHS), HY020_sum = sum(HY020))
AMp2s_sum <- AMp2s %>% group_by(cell) %>%
  summarise(HHS_sum = sum(HHS), HY020_sum = sum(HY020))

# prepare for visual comparison
gridorig <- merge(grid, AMp2_sum, by = "cell", all.x = FALSE)
gridsynt <- merge(gridorig[, -(2:3)], AMp2s_sum, by = "cell", all.x = TRUE)
gridorig$type <- "original"
gridsynt$type <- "synthetic"
grid <- rbind(gridorig, gridsynt)
grid$type <- as.factor(grid$type)

# (1) heatmaps

ggplot(grid) +
  geom_sf(aes(fill = HHS_sum), color = NA) +
  geom_sf(data = AMAP[AMAP$PROV == 2, ], fill = NA, lwd = 0.1) +
  facet_wrap(~ type) +
  theme_void() +
  theme(legend.position = "top", legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.15, "cm")) +
  scale_fill_viridis_c(name = "head count", option = "F")

ggplot(grid) +
  geom_sf(aes(fill = HY020_sum), color = NA) +
  geom_sf(data = AMAP[AMAP$PROV == 2, ], fill = NA, lwd = 0.1) +
  facet_wrap(~ type) +
  theme_void() +
  theme(legend.position = "top", legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.15, "cm")) +
  scale_fill_viridis_c(name = "disp. inc.", trans = "log10",option = "F")

# (2) scatterplots

names(gridsynt)[2:3] <- c("HHS_syn_sum", "HY020_syn_sum")
grid2 <- cbind(gridorig[, c("cell", "HHS_sum", "HY020_sum")], 
               st_drop_geometry(gridsynt[, c("HHS_syn_sum", "HY020_syn_sum")]))

ggplot(grid2, aes(HHS_sum, HHS_syn_sum)) +
  geom_abline(slope = 1, lty = "dashed", color = "red") +
  geom_point(cex = 0.2) +
  theme_minimal() +
  ggtitle("head count") +
  xlab("original") +
  ylab("synthetic")

ggplot(grid2, aes(HY020_sum, HY020_syn_sum)) +
  geom_abline(slope = 1, lty = "dashed", color = "red") +
  geom_point(cex = 0.2) +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle("disp. income") +
  xlab("original") +
  ylab("synthetic")

