#############################
#####   Creating Maps   #####
#############################

### Libraries:
##############
library(readxl);library(data.table); library(stringr)
library(tidyverse); library(sf); library(ggplot2)
library(ggspatial);library(geodist); library(here)

### Import data:
################
SampCoords.data <- read_xlsx(here::here('1_LRV_analysis/7_Maps/Sample_Coordinates.xlsx'), sheet = 1)

## Remove NA's:
SampCoords.data <- SampCoords.data[!is.na(SampCoords.data$Latitude),]

### Set CRS:
############
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"

### Convert data frames to sf objects:
########################################
SampCoords.sf <- st_as_sf(x = SampCoords.data, coords = c("Longitude", "Latitude"), crs = projcrs)

### Read in base maps of South America / Brazil: -- UNZIP FOLDER FIRST !!!
###################################################
ARG.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/ARG_adm0.shp'))
BOL.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/BOL_adm0.shp'))
BOL.admin1 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/BOL_adm1.shp'))
BRA.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/BRA_adm0.shp'))
BRA.admin1 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/BRA_adm1.shp'))
BRA.admin2 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/BRA_adm2.shp'))
CHL.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/CHL_adm0.shp'))
COL.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/COL_adm0.shp'))
ECU.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/ECU_adm0.shp'))
GUF.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/GUF_adm0.shp'))
GUY.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/GUY_adm0.shp'))
PER.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/PER_adm0.shp'))
PER.admin1 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/PER_adm1.shp'))
PRY.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/PRY_adm0.shp'))
SUR.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/SUR_adm0.shp'))
URY.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/URY_adm0.shp'))
VEN.admin0 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/VEN_adm0.shp'))
VEN.admin1 <- st_read(here::here('1_LRV_analysis/7_Maps/SHAPEFILES/VEN_adm1.shp'))

### Plot Map of Northern South America:
#######################################
SA.plot <- ggplot()+
  geom_sf(data= SA_admin) +
  geom_sf(data = BRA.admin0,  lwd= 1) +
  geom_sf(data= BRA.admin1, fill = 'gray80') + 
  geom_sf(data = VEN.admin0,  lwd= 1) + 
  geom_sf(data = VEN.admin1, fill = 'gray80') +
  xlim(-82, -34) + ylim(-55,13) +
  #ggspatial::annotation_scale(data = SA_admin, location = 'br') +
  theme_classic()

### Plotting sample presence/absence of LRV1 on map -- Supp. Fig. 1:
### All samples:
SampCoords.All <- SA.plot + geom_point(data = SampCoords.data, aes(x=Longitude, y=Latitude, shape = study, color = LRV1, size = 0.75, alpha = 0.8)) +
  scale_color_manual(values = c("#C4A681" , "#4E6E1B")) + scale_shape_manual(values = c(18, 16, 17))

### L. braziliensis samples:
SampCoords.Lb <- SA.plot + geom_point(data = subset(SampCoords.data, SampCoords.data[,4] == 'Leishmania (Viannia) braziliensis'), aes(x=Longitude, y=Latitude, shape = study, color = LRV1, size = 0.75, alpha = 0.8)) +
  scale_color_manual(values = c("#C4A681" , "#4E6E1B")) + scale_shape_manual(values = c(18, 16, 17)) + xlim(-82, -34) + ylim(-32,13) +
  ggtitle('L. braziliensis') + theme(legend.position="none")

### L. guyanensis samples:
SampCoords.Lg <- SA.plot + geom_point(data = subset(SampCoords.data, SampCoords.data[,4] == 'Leishmania (Viannia) guyanensis'), aes(x=Longitude, y=Latitude, shape = study, color = LRV1, size = 0.75, alpha = 0.8)) +
  scale_color_manual(values = c("#C4A681" , "#4E6E1B")) + scale_shape_manual(values = c(16, 17)) + xlim(-82, -34) + ylim(-32,13) +
  ggtitle('L. guyanensis') + theme(legend.position="none")

### L. naiffi samples:
SampCoords.Ln <- SA.plot + geom_point(data = subset(SampCoords.data, SampCoords.data[,4] == 'Leishmania (Viannia) naiffi'), aes(x=Longitude, y=Latitude, shape = study, color = LRV1, size = 0.75, alpha = 0.8)) +
  scale_color_manual(values = c("#C4A681" , "#4E6E1B")) + scale_shape_manual(values = c(16)) + xlim(-82, -34) + ylim(-32,13) +
  ggtitle('L. naiffi') + theme(legend.position="none")


SampCoords.All
SampCoords.Lb
SampCoords.Lg
SampCoords.Ln


### Clade distribution map - Main Fig. 1C: 
###########################

## Create pie charts:
#####################
pie.data <- read_xlsx(here::here('1_LRV_analysis/7_Maps/Sample_Coordinates.xlsx'), sheet = 2)
pie.data$Clade <- as.factor(pie.data$Clade)
levels(pie.data$Clade)
pie.data$Proportion <- as.numeric(pie.data$Proportion)
pie.data$Location <- as.factor(pie.data$Location)
pie.data$study <- as.factor(pie.data$study) 

ggplot(data=pie.data, aes(x=" ", y=Proportion, group=Clade,  fill=Clade)) + 
  scale_fill_manual(values = c("#97999c", "#97999c","#97999c","#97999c","#97999c","#97999c", "#f05922", "#333333", "#ff00ff", "#4d94cc", "#ffa30b", "#c6b099", "#0070ba", "#d71819", "#d9e7bf", "#223b50",  
                               "#395560", "#4c6459", "#005d62",  "#007b9a", "#b25f30", "#82683a", "#676926", "#80813e", "#7d9a42", "#00ae51",  "#f0a676", "#f4be3e",  "#f8bccb" )) +
  #scale_color_manual(values = rep('white',29)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) + 
  facet_grid(.~ Location) +
  theme_void() 

SA.plot <- ggplot()+
  geom_sf(data= SA_admin) +
  geom_sf(data = BRA.admin0,  lwd= 1) +
  geom_sf(data= BRA.admin1, fill = 'gray80') + 
  geom_sf(data = VEN.admin0,  lwd= 1) + 
  geom_sf(data = VEN.admin1, fill = 'gray80') +
  geom_sf(data = BOL.admin0,  lwd= 1) + 
  geom_sf(data = BOL.admin1) +
  geom_sf(data = PER.admin0,  lwd= 1) + 
  geom_sf(data = PER.admin1) +
  xlim(-82, -34) + ylim(-55,13) +
  #ggspatial::annotation_scale(data = SA_admin, location = 'br') +
  theme_classic()

Clade.plot <- SA.plot + xlim(-82, -34) + ylim(-32,13) +
  geom_point(data=pie.data, aes(x=Longitude, y=Latitude), size = 0.3)