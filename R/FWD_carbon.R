# FWD Carbon
# Ingrid Farnell
# Feb 11, 2021

# Libraries
library(dplyr)

###########################
## FWD - individual piece transects 
###########################
names(line)[names(line) == "Plot"] <- "PlotID"
names(line)[names(line) == "Line.."] <- "Transect"
names(line)[names(line) == "Horizontal.distance..adj..line.length"] <- "hor.dist"

# Simplify line
line <- line%>%
  dplyr::select(PlotID, Transect, hor.dist) # selecting columns to merge

# Transect line plot summary
line.plot <- line %>%
  group_by(PlotID) %>%
  summarize_at(c("hor.dist"), sum)

# fwd only measured on 10m of each line (=20m both lines), so subtract remaing 80 m from total length
line.plot$fwd.dist <- line.plot$hor.dist - 80

# Sum tally for each diameter class for each plot
fwd.plot <- fwd %>%
  group_by(PlotID, Diam_class) %>%
  summarise_at(c("Tally"), sum)
  
# Merge line and fwd data
fwd.line.plot <- left_join(fwd.plot, line.plot, by = c("PlotID"))

fwd.line.plot <- as.data.frame(fwd.line.plot)
fwd.line.plot$Tally <- as.numeric(as.character(fwd.line.plot$Tally))

##################################################################

# To convert volume (calculated using the VanWagner formula) to carbon you have to 
# calculate the volume for each decay class, then convert that to live biomass then dead biomass

# However, we only tallied the pieces by diameter class and did not record species or decay class
# Volume estimates will be multiplied by the bulk density of FWD and a decay-reduction factor to get biomass


# Volume equation for FWD
# CWD volume (m3/ha) = pi^2/8L  *  n * QMD^2          

# Where: 	L = length of total transect (horizontal distance (HD) in m) 
#                                       HD = SD / Square root of [1 + (% slope / 100)2]
# n = number of pieces in the diameter class
# QMD = quadratic mean diameter at the intersection aka. mean diameter of the diameter class

################################################################

#------- VOLUME ------------#
# Create QMD column (mean of diameter class)
qmdFN <- function(Diam_class){
  if(is.na(Diam_class)){
    print(paste("Diam_class is not found"))
    Diam <- NA
  } else {
    if (Diam_class == "1.1-2.5"){
      Diam <- (2.5+1.1)/2
      } else if (Diam_class == "2.6-5"){
        Diam <- (5+2.6)/2
      } else if (Diam_class == "5.1-7.5"){
        Diam <- (7.5+5.1)/2
      } else {
      print(paste("Diam_class",Species,"not found"))
      Diam <- NA
    }
  }
  return(Diam)
}

# Assign QMD to new column
qmd <- vector()
for(i in 1:nrow(fwd.line.plot)){
  qmd[i] <- qmdFN(Diam_class = fwd.line.plot[i, "Diam_class"])
}
fwd.line.plot$QMD <- qmd


# Create volumn calculation function
fwdVolFN <- function(L, n, QMD){
  fwd_vol <- ((pi^2/(8*L)) * n *(QMD^2))
  return(fwd_vol)
}

# Calculate fwd volume in each size class
fwd.line.plot$volume_ha <- fwdVolFN(L = fwd.line.plot$fwd.dist, 
                                    n = fwd.line.plot$Tally, 
                                    QMD = fwd.line.plot$QMD)

#------------ CARBON ---------------#
# fwdCARBON (T/ha) = volume(m3/ha)
#                     x Live wood density(g/cm3) # use unknown species (Harmon et al. 2008)
#                     x Decay reduction factor for each size class (Harmon and Fasth website)
#                     x CarbonConcentration # use 50% (Harmon and Fasth website)


# Assign carbon to new column
c <- vector()
for(i in 1:nrow(fwd.line.plot)){
  c[i] <- fwdCarbonFN(volume = fwd.line.plot[i,"volume_ha"], 
                      Diam_class = fwd.line.plot[i, "Diam_class"])
}
fwd.line.plot$carbon <- c

# Sum total carbon in each plot
fwd.carbon.plot <- fwd.line.plot %>%
  group_by(PlotID) %>%
  summarize_at(c("volume_ha", "carbon"), sum)


