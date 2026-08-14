#####           Co-phylogenetic analysis          #####
#####          eMPRess Cost distributions:        #####
#######################################################

### Libraries
library(ape); library(vegan)
library(paco); library(tidyverse)
library(ggridges); library(viridis)
library(readxl); library(here)

### Old World - New World:
##########################
OW.ridge.data <- readxl::read_xlsx(here::here('3_Cophylo_analysis/5_Event_based_cophylo/4_Cost_Distribution_plots/eMPRess_Ridgeline_Data.xlsx'), sheet = 1)

# basic example
ggplot(OW.ridge.data, aes(x = Events, y = EventType, fill = EventType)) +
  geom_density_ridges(jittered_points = TRUE,
                      position = position_points_jitter(width = 0, height = 0.1),
                      point_shape = '*', point_size = 4, point_alpha = 0.7, alpha = 0.8,) +
  scale_fill_manual(values = c("#05598C", "#6E868C", "#B9B88C", "#FEFEB2"))+
  theme_ridges() + 
  theme(legend.position = "none")

cols <- scico::scico(n=7, palette = 'nuuk')
"#05598C" "#386981" "#6E868C" "#A0A597" "#B9B88C" "#D1D083" "#FEFEB2"


### New World - Viannia:
########################
NW.ridge.data <- readxl::read_xlsx(here::here('3_Cophylo_analysis/5_Event_based_cophylo/4_Cost_Distribution_plots/eMPRess_Ridgeline_Data.xlsx'), sheet = 3)

# basic example
ggplot(NW.ridge.data, aes(x = Events, y = EventType, fill = EventType)) +
  geom_density_ridges(jittered_points = TRUE,
                      position = position_points_jitter(width = 0, height = 0.1),
                      point_shape = '*', point_size = 4, point_alpha = 0.7, alpha = 0.8,) +
  scale_fill_manual(values = c("#2F5D8D", "#83A0BE", "#8DB28C", "#327431"))+
  theme_ridges() + 
  theme(legend.position = "none")

cols <- scico::scico(n=7, palette = 'cork')
"#2C194C" "#2F5D8D" "#83A0BE" "#E5ECEC" "#8DB28C" "#327431" "#0E2802"


### L. braziliensis:
####################
Lb.ridge.data <- readxl::read_xlsx(here::here('3_Cophylo_analysis/5_Event_based_cophylo/4_Cost_Distribution_plots/eMPRess_Ridgeline_Data.xlsx'), sheet = 4)

# basic example
ggplot(Lb.ridge.data, aes(x = Events, y = EventType, fill = EventType)) +
  geom_density_ridges(jittered_points = TRUE,
                      position = position_points_jitter(width = 0, height = 0.1),
                      point_shape = '*', point_size = 4, point_alpha = 0.7, alpha = 0.8,) +
  scale_fill_manual(values = c("#28386A", "#2C5788", "#7396B6", "#A3B9CF"))+
  theme_ridges() + 
  theme(legend.position = "none")

cols <- scico::scico(n=14, palette = 'cork')
"#2C194C" "#28386A" "#2C5788" "#4B76A0" "#7396B6" "#A3B9CF" "#D6E0E8" 


### L. guyanensis:
##################
Lg.ridge.data <- readxl::read_xlsx(here::here('3_Cophylo_analysis/5_Event_based_cophylo/4_Cost_Distribution_plots/eMPRess_Ridgeline_Data.xlsx'), sheet = 5)

# basic example
ggplot(Lg.ridge.data, aes(x = Events, y = EventType, fill = EventType)) +
  geom_density_ridges(jittered_points = TRUE,
                      position = position_points_jitter(width = 0, height = 0.1),
                      point_shape = '*', point_size = 4, point_alpha = 0.7, alpha = 0.8,) +
  scale_fill_manual(values = c("#154B10", "#2B6F2A", "#7DA87E", "#AEC9AD"))+
  theme_ridges() + 
  theme(legend.position = "none")

cols <- scico::scico(n=14, palette = 'cork')
"#DBE7DC" "#AEC9AD" "#7DA87E" "#528B51" "#2B6F2A" "#154B10" "#0E2802"





