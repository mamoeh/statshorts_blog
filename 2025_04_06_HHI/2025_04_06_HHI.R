# ----- Herfindahl-Hirschman-Index for Country-of-birth data at municipality level -----
# M. Möhler, April 2025
# based on:
# P. Prenzel: 'Kann man das überhaupt messen?' - Der Bedarf an detaillierten 
#              räumlichen Bevölkerungsdaten in der Wirtschaftsgeographie, 2024

library(dplyr)        # data wrangling
library(readxl)       # for reading in xlsx files
library(sf)           # for working with shapefiles
library(ggplot2)      # plotting
library(RColorBrewer) # for diverging color scale


# ----- (1) Prepare Data -----

# data freely available from:
# https://ergebnisse.zensus2022.de/datenbank/online/url/9a631b03
cob_path <- "statshorts_blog/2025_04_06_HHI/1000A-1016_en_flat.csv"

# shapefile of Germany for plotting from:
# https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-250-000-stand-01-01-vg250-01-01.html
shp1_path <- "statshorts_blog/2025_04_06_HHI/vg250_01-01.utm32s.shape.ebenen/vg250_ebenen_0101/VG250_GEM.shp"
shp2_path <- "statshorts_blog/2025_04_06_HHI/vg250_01-01.utm32s.shape.ebenen/vg250_ebenen_0101/VG250_LAN.shp"

# additional classifications of municipalities by BBSR
# https://www.bbsr.bund.de/BBSR/DE/forschung/raumbeobachtung/Raumabgrenzungen/downloads/download-referenzen.html
bbsr_path <- "statshorts_blog/2025_04_06_HHI/raumgliederungen-referenzen-2022.xlsx"

# classification of municipalities to election districts & district-level election results
# https://www.bundeswahlleiterin.de/dam/jcr/aa868597-0e60-476c-bd2b-279c1e9a142a/btw25_wkr_gemeinden_20241130_utf8.csv
# https://www.bundeswahlleiterin.de/dam/jcr/f49a47a1-735b-4e9b-b4e1-4c73cad2292e/btw25_kerg2.csv
elec1_path <- "statshorts_blog/2025_04_06_HHI/20200415_btw21_wkr_gemeinden_utf8.csv"
elec2_path <- "statshorts_blog/2025_04_06_HHI/btw25_kerg2.csv"


# read in data
cob_data <- read.csv2(cob_path) %>%
  filter(X1_variable_code == "GEOGM4" & 
           value_variable_code == "PRS018" & 
           value != "-") %>%
  select(X1_variable_attribute_code, 
         X1_variable_attribute_label,
         X2_variable_attribute_code,
         X2_variable_attribute_label,
         value)

cob_data$value <- as.numeric(cob_data$value)

# discard category 'unknown / stateless' (see Prenzel, 2024, p.79)
cob_data <- cob_data %>% filter(cob_data$X2_variable_attribute_code != "LAND999")


## country level as benchmark

cob_full <- cob_data %>% group_by(X2_variable_attribute_code) %>%
  summarise(value = sum(value))

# HHI for Germany
(hhi_de <- sum((cob_full$value / sum(cob_full$value))^2))


## NUTS-1 level as benchmark

cob_data$NUTS1 <- substr(cob_data$X1_variable_attribute_code, 1, 2)
cob_nuts <- cob_data %>% group_by(NUTS1, X2_variable_attribute_code) %>%
  summarise(value = sum(value))
# make list of unique NUTS-1 IDs
hhi_nuts <- data.frame(NUTS1 = sort(unique(cob_nuts$NUTS1)), hhi_nuts = NA)
for(i in 1:nrow(hhi_nuts)) {
  # select unique NUTS-1 region
  cob_sub <- cob_nuts[cob_nuts$NUTS1 == hhi_nuts$NUTS1[i], ]
  # compute Herfindahl-Hirschman index for country of birth
  hhi_nuts$hhi_nuts[i] <- sum((cob_sub$value / sum(cob_sub$value))^2)
}


## Municipal level

# make list of unique region IDs
frac_data <- data.frame(ars = unique(cob_data$X1_variable_attribute_code), hhi = NA)
for(i in 1:nrow(frac_data)) {
  # select unique region
  cob_sub <- cob_data[cob_data$X1_variable_attribute_code == frac_data$ars[i], ]
  # compute Herfindahl-Hirschman index for country of birth
  frac_data$hhi[i] <- sum((cob_sub$value / sum(cob_sub$value))^2)
}

boxplot(frac_data$hhi) # there is a clear outlier
frac_data$hhi[which.min(frac_data$hhi)] <- NA # discard for interpretability

# quantiles of HHI, see Prenzel (2024), p.81
qtls_hhi <- round(fivenum(frac_data$hhi), 2)
frac_data$hhi_qtl <- cut(frac_data$hhi, breaks = qtls_hhi)
# HHI over benchmark
frac_data$hhi_rel <- frac_data$hhi - hhi_de


# ----- (2) Prepare plotting -----

# read in shapefile and join with index data
de_gem <- read_sf(shp1_path) %>%
  filter(GF == 4) # filter for land area only
de_gem <- merge(x = de_gem, y = frac_data, 
                by.x = "ARS_0", by.y = "ars", 
                all.x = TRUE, all.y = FALSE)

# add NUTS-1-level benchmarks
de_gem <- merge(x = de_gem, y = hhi_nuts,
                by.x = "SN_L", by.y = "NUTS1",
                all.x = TRUE)
de_gem$hhi_rel_nuts <- de_gem$hhi - de_gem$hhi_nuts


# read in shapefile for NUTS-1 borders (for plotting purposes only)
de_lan <- read_sf(shp2_path) %>% filter(GF == 4)


# ----- (3) Plot maps -----

# HHI at municipal level (continuous)
ggplot(de_gem) +
  geom_sf(aes(fill = hhi), color = NA) +
  geom_sf(data = de_lan, color = "black", fill = NA) +
  scale_fill_viridis_c(name = "HHI", na.value = "grey50") +
  theme_void() +
  theme(legend.key.width = unit(.25, "cm"))
#ggsave("HHI_cont.png", width = 800, height = 1000, units = "px")

# HHI at municipal level (quantiles)
ggplot(de_gem) +
  geom_sf(aes(fill = hhi_qtl), color = NA) +
  geom_sf(data = de_lan, color = "black", fill = NA) +
  scale_fill_viridis_d(name = "HHI (Quantile)", na.value = "grey50") +
  theme_void() +
  theme(legend.key.width = unit(.25, "cm"))
#ggsave("HHI_quant.png", width = 1000, height = 1000, units = "px")

# HHI at municipal level (difference from country average)
div_cols <- rev(RColorBrewer::brewer.pal(11, name = "PiYG"))[c(1, 6, 11)]
ggplot(de_gem) +
  geom_sf(aes(fill = hhi_rel), color = NA) +
  geom_sf(data = de_lan, color = "black", fill = NA) +
  scale_fill_gradient2(low = div_cols[1], mid = div_cols[2], high = div_cols[3],
                       midpoint = 0, name = "Diff.") +
  theme_void() +
  theme(legend.key.width = unit(.25, "cm"))
#ggsave("HHI_diff1.png", width = 800, height = 1000, units = "px")

# HHI at municipal level (difference from NUTS-1 average)
ggplot(de_gem) +
  geom_sf(aes(fill = hhi_rel_nuts), color = NA) +
  geom_sf(data = de_lan, color = "black", fill = NA) +
  scale_fill_gradient2(low = div_cols[1], mid = div_cols[2], high = div_cols[3],
                       midpoint = 0, name = "Diff.") +
  theme_void() +
  theme(legend.key.width = unit(.25, "cm"))
#ggsave("HHI_diff2.png", width = 800, height = 1000, units = "px")


# ----- (4) Further analyses ----

## Special aggregations (BBSR)

# read in BBSR municipality classifications
bbsr_raumg <- readxl::read_xlsx(bbsr_path, sheet = "Gemeindereferenz (inkl. Kreise)") %>%
  select(GEM2022_RS, GTU2022, GTU_NAME, RLG2022, RLG_NAME, GWS2022, GWS_NAME)

# correcting leading zeros
bbsr_raumg$GEM2022_RS <- ifelse(nchar(bbsr_raumg$GEM2022_RS) == 11,
                                paste0("0", bbsr_raumg$GEM2022_RS), 
                                bbsr_raumg$GEM2022_RS)

# merge BBSR categories to COB-data
cob_data <- merge(x = cob_data, y = bbsr_raumg, 
                  by.x = "X1_variable_attribute_code", by.y = "GEM2022_RS",
                  all.x = TRUE, all.y = FALSE)

# aggregate by BBSR categories
cob_gtu <- cob_data %>% group_by(GTU2022, X2_variable_attribute_code) %>%
  summarise(value = sum(value), GTU_NAME = unique(GTU_NAME))
cob_rlg <- cob_data %>% group_by(RLG2022, X2_variable_attribute_code) %>%
  summarise(value = sum(value), RLG_NAME = unique(RLG_NAME))
cob_gws <- cob_data %>% group_by(GWS2022, X2_variable_attribute_code) %>%
  summarise(value = sum(value), GWS_NAME = unique(GWS_NAME))

# make list of unique categories
hhi_gtu <- data.frame(GTU2022 = unique(cob_gtu$GTU2022),
                      GTU_NAME = unique(cob_gtu$GTU_NAME), hhi_gtu = NA)
hhi_rlg <- data.frame(RLG2022 = unique(cob_rlg$RLG2022),
                      RLG_NAME = unique(cob_rlg$RLG_NAME), hhi_rlg = NA)
hhi_gws <- data.frame(GWS2022 = unique(cob_gws$GWS2022),
                      GWS_NAME = unique(cob_gws$GWS_NAME), hhi_gws = NA)

# compute HHIs
for(i in 1:nrow(hhi_gtu)) {
  cob_sub <- cob_gtu[cob_gtu$GTU2022 == hhi_gtu$GTU2022[i], ]
  hhi_gtu$hhi_gtu[i] <- sum((cob_sub$value / sum(cob_sub$value))^2)
}
for(i in 1:nrow(hhi_rlg)) {
  cob_sub <- cob_rlg[cob_rlg$RLG2022 == hhi_rlg$RLG2022[i], ]
  hhi_rlg$hhi_rlg[i] <- sum((cob_sub$value / sum(cob_sub$value))^2)
}
for(i in 1:nrow(hhi_gws)) {
  cob_sub <- cob_gws[cob_gws$GWS2022 == hhi_gws$GWS2022[i], ]
  hhi_gws$hhi_gws[i] <- sum((cob_sub$value / sum(cob_sub$value))^2)
}

# benchmark against national HHI
hhi_gtu$hhi_gtu_diff <- hhi_gtu$hhi_gtu - hhi_de
hhi_rlg$hhi_rlg_diff <- hhi_rlg$hhi_rlg - hhi_de
hhi_gws$hhi_gws_diff <- hhi_gws$hhi_gws - hhi_de

hhi_gtu
hhi_rlg
hhi_gws


## Election results

# read in electoral district classification
elec_class <- read.csv2(elec1_path, skip = 7, colClasses = "character")
# correct ARS parts
elec_class$RGS_GemVerband <- ifelse(elec_class$RGS_GemVerband == "0000",
                                    paste0("0", elec_class$RGS_Gemeinde),
                                    elec_class$RGS_GemVerband)
# create ARS
elec_class$ARS <- paste0(elec_class$RGS_Land, elec_class$RGS_RegBez, elec_class$RGS_Kreis,
                         elec_class$RGS_GemVerband, elec_class$RGS_Gemeinde)

# merge electoral districts to COB-data
cob_data <- merge(x = cob_data, 
                  y = elec_class[!duplicated(elec_class$ARS), c("Wahlkreis.Nr", "ARS")],
                  by.x = "X1_variable_attribute_code", by.y = "ARS",
                  all.x = TRUE, all.y = FALSE)

# aggregate by electoral district
cob_dist <- cob_data %>% group_by(Wahlkreis.Nr, X2_variable_attribute_code) %>%
  summarise(value = sum(value))
cob_dist <- cob_dist[!is.na(cob_dist$Wahlkreis.Nr), ]

# make list of unique electoral district IDs
hhi_dist <- data.frame(Wahlkreis.Nr = unique(cob_dist$Wahlkreis.Nr), hhi_dist = NA)
# compute HHIs
for(i in 1:nrow(hhi_dist)) {
  cob_sub <- cob_dist[cob_dist$Wahlkreis.Nr == hhi_dist$Wahlkreis.Nr[i], ]
  hhi_dist$hhi_dist[i] <- sum((cob_sub$value / sum(cob_sub$value))^2)
}

# read in election results
cClass <- c(rep("character", 9), rep("numeric", 8), rep("character", 2))
elec_res <- read.csv2(elec2_path, skip = 9, colClasses = cClass) %>% 
  filter(Gebietsart == "Wahlkreis" & Stimme == 2 & Gruppenart == "Partei") %>%
  select(Gebietsnummer, Gruppenart, Gruppenname, Anzahl)

# municipalities with several electoral districts need to have theirs combined
# in the results data (since we only have HHI at municipality level)
double_dist <- unique(elec_class$ARS[duplicated(elec_class$ARS)])

for(i in seq(double_dist)) {
  # re-name electoral districts to reflect combination
  ed_nrs <- elec_class$Wahlkreis.Nr[elec_class$ARS == double_dist[i]]
  elec_res$Gebietsnummer[elec_res$Gebietsnummer %in% ed_nrs] <- ed_nrs[1]
}

# re-compute percentages
elec_res <- elec_res %>% group_by(Gebietsnummer, Gruppenname) %>%
  summarise(Anzahl = sum(Anzahl, na.rm = TRUE))
elec_sums <- elec_res %>% group_by(Gebietsnummer) %>%
  summarise(Total = sum(Anzahl))
elec_res <- merge(elec_res, elec_sums, by = "Gebietsnummer")
elec_res$Prozent <- (elec_res$Anzahl / elec_res$Total) * 100

# merge election results and electoral districts HHIs
elec_res <- merge(x = elec_res, y = hhi_dist,
                  by.x = "Gebietsnummer", by.y = "Wahlkreis.Nr",
                  all.x = TRUE, all.y = FALSE)

# subset to relevant parties
elec_res <- elec_res %>% filter(Gruppenname %in% c("AfD", "CDU", "FDP", "SPD", 
                                                   "GRÜNE", "BSW"))
party_cols <- c("#009ee0", "#5F316E", "#151518", "#F9E03A", "#409A3C", "#E3000F")

ggplot(elec_res, aes(hhi_dist, Prozent, color = Gruppenname)) +
  geom_point(show.legend = FALSE, size = .4) +
  scale_color_manual(values = party_cols) +
  facet_wrap(~Gruppenname, nrow = 2) +
  xlab("HHI") +
  theme_bw()
#ggsave("HHI_election.png", width = 1500, height = 1200, units = "px")

