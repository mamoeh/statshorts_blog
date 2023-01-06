# ----- Best practices Datenvisualisierung -----
# M. Möhler, Jan. 2023
# Hinweis: Die verwendeten Daten sind manuell aus den Beispieldiagrammen im jew.
#          Original abgelesen; sie sind ausschließlich zur Übung vorgesehen und 
#          sollten nicht inhaltlich interpretiert werden.


library(readxl)       # zum Einlesen von xlsx-Dateien
library(ggplot2)      # zum Erstellen von Grafiken
library(directlabels) # für intelligente Labelling-Methoden
library(khroma)       # für barrierefreie Farbskalen


# ----- JWB 2022, S.41 -----

# einlesen
ladep <- read_xlsx("2023_01_06_DataViz/JWB_eyeball_data.xlsx", sheet = "2022_41")

# formatieren
ladep$Jahr <- as.Date(ladep$Jahr)

# Farben festlegen
sb5Blau <- rgb(0, 115/255, 190/255)
sb5Text <- "grey30"

# Grafik
ggplot(ladep, aes(Jahr, Ladep)) +
  geom_bar(stat = "identity", fill = sb5Blau, width = 150) +
  geom_hline(yintercept = (1:5)*1e+4, color = "white") +
  geom_text(aes(label = format(Ladep, big.mark = ".", decimal.mark = ",")), 
            nudge_y = 3e+3, cex = 2.5) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme_minimal() +
  theme(line = element_blank(),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.title    = element_text(color = sb5Text),
        plot.subtitle = element_text(color = sb5Text),
        plot.caption  = element_text(color = sb5Text)) +
  guides(y = "none") + 
  ylab(NULL) +
  xlab(NULL) +
  labs(title    = "Schaubild 5: Anzahl der Elektro-Ladepunkte in Deutschland",
       subtitle = "(Stichtag jew. 31.12.)",
       caption  = "Quelle: Bundesnetzagentur")


# ----- JWB 2022, S.70 -----

# einlesen
arbmrkt <- read_xlsx("2023_01_06_DataViz/JWB_eyeball_data.xlsx", sheet = "2022_70")

# formatieren
arbmrkt$Monat <- as.Date(arbmrkt$Monat)

# indexieren
basis <- as.Date("2020-03-01") # wähle Basismonat (= 100)
arbmrkt$Arbeitsl_ind   <- arbmrkt$Arbeitsl   / arbmrkt$Arbeitsl[arbmrkt$Monat == basis]   * 100
arbmrkt$Unterbesch_ind <- arbmrkt$Unterbesch / arbmrkt$Unterbesch[arbmrkt$Monat == basis] * 100
arbmrkt$Beschaeft_ind  <- arbmrkt$Beschaeft  / arbmrkt$Beschaeft[arbmrkt$Monat == basis]  * 100

# zum Plotten aufbereiten
arb <- data.frame(Monat = rep(arbmrkt$Monat, 3),
                  Wert  = c(arbmrkt$Arbeitsl, arbmrkt$Unterbesch, arbmrkt$Beschaeft),
                  Index = c(arbmrkt$Arbeitsl_ind, arbmrkt$Unterbesch_ind, arbmrkt$Beschaeft_ind),
                  Typ   = rep(c("Arbeitslosigkeit", "Unterbeschäftigung", "SV Beschäftigung"),
                              each = nrow(arbmrkt)))
arb$Typ <- factor(arb$Typ, levels = unique(arb$Typ),
                  labels = c("Arbeitslosigkeit", "Unterbeschäftigung (ohne Kurzarbeit)",
                             "SV-pflicht. Beschäftigung"))

# abs. Veränderungen berechnen
arb$lag <- c(NA, arbmrkt$Arbeitsl[1:(nrow(arbmrkt) - 1)],
             NA, arbmrkt$Unterbesch[1:(nrow(arbmrkt) - 1)],
             NA, arbmrkt$Beschaeft[1:(nrow(arbmrkt) - 1)])
arb$Aend <- arb$Wert - arb$lag
  
# Farben festlegen
sb9Text <- "grey30"
sb9cols <- c(rgb(25/255, 130/255, 195/255),
             rgb(5/255, 70/255, 120/255),
             rgb(235/255, 115/255, 80/255))

# label-Punkte für die Datumsachse
brks  <- arbmrkt$Monat[seq(1, 36, 3)] 
brks2 <- arbmrkt$Monat[seq(2, 36, 3)]

#Grafik - Version 1
ggplot(arb[arb$Monat != min(arb$Monat), ], aes(Monat, Aend/1e+5, fill = Typ)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sb9cols) +
  scale_x_date(breaks = brks2, date_labels = "'%y-%m", date_minor_breaks = "1 month") +
  scale_y_continuous(n.breaks = 6) +
  theme_minimal() +
  labs(title    = "Schaubild 9: Saisonbereinigte Entwicklung des Arbeitsmarkts",
       subtitle = "(Veränderung zum Vormonat)",
       caption  = "Quelle: Statistik der Bundesagentur für Arbeit") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title    = element_text(color = sb9Text),
        plot.subtitle = element_text(color = sb9Text),
        plot.caption  = element_text(color = sb9Text),
        axis.text.x = element_text(angle = 90)) +
  ylab("100 Tsd. Personen") +
  xlab(NULL) +
  guides(fill = guide_legend(nrow = 2))

# Grafik - Version 2
ggplot(arb, aes(Monat, Index, color = Typ)) +
  geom_line() +
  scale_color_manual(values = sb9cols) +
  scale_x_date(breaks = brks, date_labels = "'%y-%m", date_minor_breaks = "1 month") +
  theme_minimal() +
  labs(title    = "Schaubild 9: Saisonbereinigte Entwicklung des Arbeitsmarkts",
       subtitle = "(März 2020 = 100)",
       caption  = "Quelle: Statistik der Bundesagentur für Arbeit") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title    = element_text(color = sb9Text),
        plot.subtitle = element_text(color = sb9Text),
        plot.caption  = element_text(color = sb9Text),
        axis.text.x = element_text(angle = 90)) +
  ylab(NULL) +
  xlab(NULL) +
  guides(color = guide_legend(nrow = 2))


# ----- JWB 2020, S.10 -----

# einlesen
trbhg <- read_xlsx("2023_01_06_DataViz/JWB_eyeball_data.xlsx", sheet = "2020_10")

# formatieren
trbhg$Jahr <- as.Date(as.character(trbhg$Jahr), format = "%Y")
trbhg$Gebiet <- as.factor(trbhg$Gebiet)
trbhg$Gebiet_lang <- factor(trbhg$Gebiet, levels = unique(trbhg$Gebiet),
                            labels = c("USA", "EU", "Russland", "Sonstige", 
                                       "Asien & Ozeanien (ohne China, Indien)",
                                       "China", "Afrika", "Indien"))

# Farben festlegen
sb1Text <- "grey30"

# Grafik
ggplot(trbhg, aes(Jahr, Emissionen_MioT/1000, color = Gebiet_lang, group = Gebiet)) +
  geom_line() +
  labs(title    = "Schaubild 1: Globale Treibhausgasemissionen",
       subtitle = "(CO2-Äquivalent)",
       caption  = "Quelle: Internationale Energieagentur") +
  ylab("Milliarden Tonnen") +
  xlab(NULL) +
  khroma::scale_color_muted() +
  scale_y_continuous(n.breaks = 10, limits = c(0, 10)) +
  scale_x_date(breaks = as.Date(as.character(seq(1970, 2015, 5)), format = "%Y"), 
               date_labels = "%Y") +
  geom_dl(aes(label = Gebiet), method = "top.points") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title    = element_text(color = sb1Text),
        plot.subtitle = element_text(color = sb1Text),
        plot.caption  = element_text(color = sb1Text))

