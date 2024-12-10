#Remote Sensing Landsat imagery - Spermonde Archipelago summer 2022
#Author: Elsa B. Girard


# check and clean memory
ls()
rm(list=ls())
gc()
ls()

#load packages
library(raster)
library(stringr)
library(ncdf4)
library(viridis)
library(landsat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(readxl)
library(naniar)
library(reshape2)
library(lubridate)
library(padr)
library(geosphere)
library(ggmap)
library(terra)
library(patchwork)

"%not%" <- Negate("%in%") #to create a "not in" sign, the opposite of %in%


#----------------load geotiff files and add extent (test on one file)---------------------------

#set local path
setwd("~downloads/Dataset_OLCI_Sentinel3/cropped_files_sentinel_3")
getwd()

#load geotiff
data <- stack("OLL2WFR_20220901T020947Z_20220901T021006Z_epct_2e10c043_FC.tif")

#add extent to the geotiff
bb <- extent(118.6, 119.7, -5.3, -4.4) #extent(xmin,xmax,ymin,ymax)
extent(data) <- bb
data <- setExtent(data, bb)

#new extent for cropping (lat, long for Spermonde)
e <- extent(118.8, 119.55, -5.3, -4.6)
#crop landsat by the extent
data <- crop(data, e)

#north-south inverted, so flip data
data <- flip(data)

#plot
plot(data)

#plot TSM
tsm <- data[[1]] #log10 scale, units g/m3, range: 0-100
tsm_10 <- 10 ^ tsm
tsm_10[tsm_10 > 70] <- NA #based on highest turbidity known from deltas (higher means land)
tsm_log <- log10(tsm_10)

plot(tsm_log, col= rev(hcl.colors(100, "Roma")))

#plot Chlorophyll
chl <- data[[3]] #log10 scale, unites mg(chl-a)/m3, range: 0.01-100
plot(chl, col= hcl.colors(100, "TealGrn"))



#----------------TSM average over July 1st to July 31st---------------------------

#set local path
setwd("~downloads/Dataset_OLCI_Sentinel3/cropped_files_sentinel_3")
getwd()

files <- list.files("~downloads/Dataset_OLCI_Sentinel3/cropped_files_sentinel_3",
                    pattern = ".tif")


bb <- extent(13202492, 13324943, -590836.5, -490287.9) #extent(xmin,xmax,ymin,ymax) - ORIGINAL
oo <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs"
res <- c(300.1261, 300.1451) #of first file
dim <- c(335, 408, 136680) #of first file

tsm_stack <- raster()
extent(tsm_stack) <- bb
tsm_stack <- setExtent(tsm_stack, bb)
res(tsm_stack) <- res
dim(tsm_stack) <- dim
crs(tsm_stack) <- oo

chl_stack <- raster()
extent(chl_stack) <- bb
chl_stack <- setExtent(chl_stack, bb)
res(chl_stack) <- res
dim(chl_stack) <- dim
crs(chl_stack) <- oo

for (i in 1:length(files)){
  
  data <- brick(paste("~downloads/Dataset_OLCI_Sentinel3/cropped_files_sentinel_3/",
                      files[[i]], sep = ""))

  tsm <- data[[1]] #log10 scale, units g/m3, range: 0-100
  tsm <- resample(tsm, tsm_stack, method='ngb')
  tsm_stack <- addLayer(tsm_stack, tsm)
  
  chl <- data[[3]] #log10 scale, unites (mg(chl-a)/m3), range: 0.01-100
  chl <- resample(chl, chl_stack, method='ngb')
  chl_stack <- addLayer(chl_stack, chl)
  
}

#plot TSM
tsm_mean <- calc(tsm_stack, fun = median, na.rm = T)
plot(tsm_mean, col= rev(hcl.colors(100, "Roma")))

#plot Chlorophyll
chl_mean <- calc(chl_stack, fun = median, na.rm = T)
plot(chl_mean , col= hcl.colors(100, "TealGrn"))


#----------------land masking--------------------------------

#load geotiff
lmask <- brick("~downloads/Dataset_OLCI_Sentinel3/landmask/OLL1EFR_20190910T015050Z_20190910T015108Z_epct_cf675526_FPRC.tif")


#writeRaster(lmask, '~downloads/Dataset_OLCI_Sentinel3/landmask/forland.tif',overwrite=TRUE)

color <- stack(lmask[[c(18)]])
x <- 6
color[color <= x] <- 0
color[color > x] <- NA
plot(color,col = c("white", "black"))

land_mask <- color

#writeRaster(land_mask, '~downloads/Dataset_OLCI_Sentinel3/landmask/landmask_spermonde_sentinel.tif', overwrite = TRUE)

ee <- extent(118.9, 119.6, -5.2, -4.5) #extent(xmin,xmax,ymin,ymax) - zoom Spermonde

#apply land mask to tsm mean
land_mask <- resample(land_mask, tsm_mean, method='ngb')
tsm_masked <- mask(tsm_mean, mask = land_mask)
tsm_masked_10 <- 10 ^ tsm_masked
tsm_masked_10[tsm_masked_10 > 70] <- NA #based on highest turbidity known from deltas (higher means land)
tsm_masked_log <- log10(tsm_masked_10)
tsm_masked_log <- projectRaster(tsm_masked_log, crs = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
tsm_masked_log <- crop(tsm_masked_log, ee)
plot(tsm_masked_log, col= rev(hcl.colors(100, "Roma")))

#apply land mask to chl mean
land_mask <- resample(land_mask, chl_mean, method='ngb')
chl_masked <- mask(chl_mean, mask = land_mask)
chl_masked <- projectRaster(chl_masked, crs = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
chl_masked <- crop(chl_masked, ee)
plot(chl_masked , col= hcl.colors(100, "TealGrn"))

#writeRaster(chl_masked, '~downloads/Dataset_OLCI_Sentinel3/final/chl_masked.tif', overwrite = TRUE)

#----------------Fig. 1a,b --- draw circles around the islands---------------------

map_tsm <- raster("~downloads/Dataset_OLCI_Sentinel3final/tsm_masked_log.tif")
map_chl <- raster("~downloads/Dataset_OLCI_Sentinel3/final/chl_masked.tif")

data <- read.csv("~downloads/Dataset_OLCI_Sentinel3/coordinates_islands.csv")

#function to create circles   
make_circle <- function (x, y, r, n = 360, proj4str) 
{
  pts <- seq(0, 2 * pi, length.out = n)
  xy <- cbind(x + r * sin(pts), y + r * cos(pts))
  sl <- SpatialLines(list(Lines(list(Line(xy)), "line")))
  if (!missing(proj4str)) 
    sp::proj4string(sl) <- proj4str
  return(sl)
}

spc <- make_circle(x= data$lon[1],y = data$lat[1],r=0.015)
circle1 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle1)
L1 <- Lines(list(circle_L), ID = data$ID[1])
circle1$island <- data$ID[1]

spc <- make_circle(x= data$lon[2],y = data$lat[2],r=0.015)
circle2 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle2)
L2 <- Lines(list(circle_L), ID = data$ID[2])
circle2$island <- data$ID[2]

spc <- make_circle(x= data$lon[3],y = data$lat[3],r=0.015)
circle3 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle3)
L3 <- Lines(list(circle_L), ID = data$ID[3])
circle3$island <- data$ID[3]

spc <- make_circle(x= data$lon[4],y = data$lat[4],r=0.015)
circle4 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle4)
L4 <- Lines(list(circle_L), ID = data$ID[4])
circle4$island <- data$ID[4]

spc <- make_circle(x= data$lon[5],y = data$lat[5],r=0.015)
circle5 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle5)
L5 <- Lines(list(circle_L), ID = data$ID[5])
circle5$island <- data$ID[5]

spc <- make_circle(x= data$lon[6],y = data$lat[6],r=0.015)
circle6 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle6)
L6 <- Lines(list(circle_L), ID = data$ID[6])
circle6$island <- data$ID[6]

spc <- make_circle(x= data$lon[7],y = data$lat[7],r=0.015)
circle7 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle7)
L7 <- Lines(list(circle_L), ID = data$ID[7])
circle7$island <- data$ID[7]

spc <- make_circle(x= data$lon[8],y = data$lat[8],r=0.015)
circle8 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle8)
L8 <- Lines(list(circle_L), ID = data$ID[8])
circle8$island <- data$ID[8]

spc <- make_circle(x= data$lon[9],y = data$lat[9],r=0.015)
circle9 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle9)
L9 <- Lines(list(circle_L), ID = data$ID[9])
circle9$island <- data$ID[9]

spc <- make_circle(x= data$lon[10],y = data$lat[10],r=0.015)
circle10 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle10)
L10 <- Lines(list(circle_L), ID = data$ID[10])
circle10$island <- data$ID[10]

spc <- make_circle(x= data$lon[11],y = data$lat[11],r=0.015)
circle11 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle11)
L11 <- Lines(list(circle_L), ID = data$ID[11])
circle11$island <- data$ID[11]

spc <- make_circle(x= data$lon[12],y = data$lat[12],r=0.015)
circle12 <- as.data.frame(geom(spc)[,c(4,5)])
circle_L <- Line(circle12)
L12 <- Lines(list(circle_L), ID = data$ID[12])
circle12$island <- data$ID[12]

all_circles <- SpatialLines(list(L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12))
coordinates_circles <- rbind(circle1,circle2,circle3,circle4,circle5,circle6,circle7,circle8,circle9,circle10,circle11,circle12)


#TSM
plot(map_tsm , col= rev(hcl.colors(100, "Roma")))
points(data = data, x = data$lon, y = data$lat, pch = 16, col = "red")
text(data = data,  x = data$lon, y = data$lat, labels = data$ID, pos = 4, offset = 0.7)
plot(all_circles, add = TRUE, col = "red")

#Chl
plot(map_chl , col= hcl.colors(100, "TealGrn"))
points(data = data, x = data$lon, y = data$lat, pch = 16, col = "red")
text(data = data,  x = data$lon, y = data$lat, labels = data$ID, pos = 4, offset = 0.7)
plot(all_circles, add = TRUE, col = "red")



#----------------Extract data value from coordinates------------------------

coordinates(coordinates_circles)= ~ x + y

#_____TSM________________
rasValue <- extract(map_tsm, coordinates_circles)
combinePointValue <- as.data.frame(cbind(coordinates_circles,rasValue))
colnames(combinePointValue) <- c("ID", "tsm_log", "x", "y")

data_adj_tsm <- combinePointValue


#____________Chl___________________________________

rasValue <- extract(map_chl, coordinates_circles)
combinePointValue <- as.data.frame(cbind(coordinates_circles,rasValue))
colnames(combinePointValue) <- c("ID", "chl_log", "x", "y")

data_adj_chl <- combinePointValue


#----------------Supp. Fig S1 --- check violin plot of the points--------------------------


#remove data points above the 90 (and 50) percentile to reduce bottom reflectance
data_adj_tsm$tsm <- 10^data_adj_tsm$tsm_log 
data_adj_chl$chl <- 10^data_adj_chl$chl_log 

filtered_tsm <-data_adj_tsm %>% group_by(ID) %>% mutate(reduced_tsm = ifelse(tsm > quantile(tsm, 0.90, na.rm=TRUE), NA, tsm))
filtered_tsm$reduced_tsm[filtered_tsm$ID == 'Langkai' & filtered_tsm$tsm[filtered_tsm$ID == 'Langkai'] > quantile(filtered_tsm$tsm[filtered_tsm$ID == 'Langkai'], 0.50, na.rm=TRUE)] <- NA
filtered_tsm$reduced_tsm[filtered_tsm$ID == 'Kapoposang Masdar point' & filtered_tsm$tsm[filtered_tsm$ID == 'Kapoposang Masdar point'] > quantile(filtered_tsm$tsm[filtered_tsm$ID == 'Langkai'], 0.50, na.rm=TRUE)] <- NA
filtered_tsm$color <- 'black'
filtered_tsm$color[is.na(filtered_tsm$reduced_tsm)] <- 'red'

filtered_chl <-data_adj_chl %>% group_by(ID) %>% mutate(reduced_chl = ifelse(chl > quantile(chl, 0.90, na.rm=TRUE), NA, chl))
filtered_chl$reduced_chl[filtered_chl$ID == 'Langkai' & filtered_chl$chl[filtered_chl$ID == 'Langkai'] > quantile(filtered_chl$chl[filtered_chl$ID == 'Langkai'], 0.50, na.rm=TRUE)] <- NA
filtered_chl$reduced_chl[filtered_chl$ID == 'Kapoposang Masdar point' & filtered_chl$chl[filtered_chl$ID == 'Kapoposang Masdar point'] > quantile(filtered_chl$chl[filtered_chl$ID == 'Langkai'], 0.50, na.rm=TRUE)] <- NA
filtered_chl$color <- 'black'
filtered_chl$color[is.na(filtered_chl$reduced_chl)] <- 'red'

ID_ord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
            "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
            "Kapoposang Masdar point")

p1 <- ggplot(filtered_tsm[filtered_tsm$color == 'black',], aes(x = factor(ID, rev(ID_ord)), y = tsm)) +
  gghalves::geom_half_boxplot() +
  gghalves::geom_half_point(color = 'black', alpha = 0.3) +
  gghalves::geom_half_point(data = filtered_tsm[filtered_tsm$color == 'red',], aes(x = factor(ID, rev(ID_ord)), y = tsm), color = 'red', alpha = 0.3) +
  scale_y_continuous(breaks=seq(-1, 50, 1)) +
  scale_y_log10() + theme_classic()  +
  labs(tag = 'A')

p2 <- ggplot(filtered_chl[filtered_chl$color == 'black',], aes(x = factor(ID, rev(ID_ord)), y = chl)) +
  gghalves::geom_half_boxplot() +
  gghalves::geom_half_point(color = 'black', alpha = 0.3) +
  gghalves::geom_half_point(data = filtered_chl[filtered_chl$color == 'red',], aes(x = factor(ID, rev(ID_ord)), y = chl), color = 'red', alpha = 0.3) +
  scale_y_continuous(breaks=seq(-1, 50, 1)) +
  scale_y_log10() + theme_classic() +
  labs(tag = 'B')

layout <- "
AB
"

p1 + p2 + plot_layout(design = layout)


summary_data  <- list(data, filtered_mean_chl[,c(1,2)],filtered_mean_tsm[,c(1,2)])     
summary_data  <- summary_data %>% purrr::reduce(full_join, by='ID')

colnames(summary_data) <- c("island", "longitude", "latitude", "mean_chl_filtered", "mean_tsm_filtered")

write.csv(summary_data, "~downloads/Dataset/filtered_tsm_chl_mean_islands.csv", row.names = FALSE)



#----------------Fig. 1c --- TSM vs Chl-a in Spermonde ---------------------
new_df <- summary_data

new_df$dist_coast_km[new_df$island == "Langkai"] <- 42
new_df$dist_coast_km[new_df$island == "Samalona"] <- 7
new_df$dist_coast_km[new_df$island == "Badi"] <- 23
new_df$dist_coast_km[new_df$island == "K. Keke"] <- 15
new_df$dist_coast_km[new_df$island == "Langkadea"] <- 14
new_df$dist_coast_km[new_df$island == "Karanrang"] <- 14
new_df$dist_coast_km[new_df$island == "Barang Baringan"] <- 5
new_df$dist_coast_km[new_df$island == "Lae Lae"] <- 1
new_df$dist_coast_km[new_df$island == "Lumulumu"] <- 31
new_df$dist_coast_km[new_df$island == "Polewali"] <- 11
new_df$dist_coast_km[new_df$island == "Padjenekang"] <- 19
new_df$dist_coast_km[new_df$island == "Kapoposang Masdar point"] <- 62
new_df$dist_coast_km[new_df$island == "Barang Lompo"] <- 17
new_df$dist_coast_km[new_df$island == "Karang kassi"] <- 39
new_df$dist_coast_km[new_df$island == "Bonetambung"] <- 22

ggplot(new_df, aes(x = dist_coast_km)) +
  geom_point(aes(y = mean_chl_filtered), col = "aquamarine4") +
  stat_smooth(aes(y = mean_chl_filtered), method = "lm", formula = y ~ log(x), geom = "smooth", col = "aquamarine4", se = FALSE) +
  geom_point(aes(y = mean_tsm_filtered), col = "salmon4") +
  stat_smooth(aes(y = mean_tsm_filtered), method = "lm", formula = y ~ log(x), geom = "smooth", col = "salmon4", se = FALSE) +
  theme_classic() +
  theme(
    axis.title.y = element_text(color = "aquamarine4"),
    axis.title.y.right = element_text(color = "salmon4"))+
  xlab("Distance to coast (km)")+
  scale_y_log10(name = "Chlorophyll-a (mg/m3)",
                sec.axis = sec_axis(trans=~., name="Total suspended matter (g/m3)"))

cor.test(new_df$dist_coast_km, new_df$mean_chl_filtered, method="kendall")
cor.test(new_df$dist_coast_km, new_df$mean_chl_filtered, method="spearman")
cor.test(new_df$dist_coast_km, new_df$mean_tsm_filtered, method="kendall")
cor.test(new_df$dist_coast_km, new_df$mean_tsm_filtered, method="spearman")
