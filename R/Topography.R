#Landscape topography data
library(raster)
library(spatialEco)
library(sf)
library(data.table)

####### read in DEM#######
dem_studyFires <- raster("./Inputs/SpatialFiles/dem_sa_nad1.tif") #projected DEM
dem <- raster("./Inputs/SpatialFiles/dem_sa.tif") #unprojected dem raster

####################################
#create topography layers
####################################
#slope
DEMslope <- raster::terrain(dem_studyFires, opt=c("slope"))
writeRaster(DEMslope,"./Inputs/SpatialFiles/DEMslope.tif",overwrite=TRUE)

#Topographic position index
DEMtpi <- tpi(dem_studyFires, win="circle", scale=100) #not sure what the scale is
writeRaster(DEMtpi,"./Inputs/SpatialFiles/DEMtpi.tif")

#heat load index (based on McCune and Keon 2002)
DEMhli <- hli(dem, check=TRUE,force.hemisphere = "northern")
DEMhlinad <- projectRaster(DEMhli,crs="+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
writeRaster(DEMhlinad,"./Inputs/SpatialFiles/DEMhli.tif", overwrite=TRUE)
