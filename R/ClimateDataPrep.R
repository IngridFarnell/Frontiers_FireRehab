library(raster)
library(spatialEco)
library(sf)
library(data.table)
##### DEM #####
####### read in DEM - not projected #######
dem <- raster("E:/Spatial Data/DEM/BC_dem.tif") #this is the raster for the whole province 

####################################
## ClimWNA - this is how I got the dem for ClimateNA
Study_plots <- read_sf("./Inputs/SpatialFiles/FRplots.shp")
Study_plots_latLong <- st_transform(Study_plots,crs="+proj=longlat +datum=WGS84 +no_defs")
plot_geoms <- do.call(rbind, st_geometry(Study_plots_latLong))
plot_dat <- as.data.table(Study_plots_latLong)
plot_dat[,geometry:=NULL]
plot_dat<- cbind(plot_dat,plot_dat,plot_geoms)
plot_dat<- setNames(plot_dat,c("ID1","ID2","long","lat"))
Plot_elev <- raster::extract(dem,Study_plots_latLong,method="simple")
plot_dat <- plot_dat[,.(ID1,ID2,lat,long,elev=Plot_elev)]
write.csv(plot_dat,"./Inputs/FR_plots_Clim.csv",row.names=FALSE)
