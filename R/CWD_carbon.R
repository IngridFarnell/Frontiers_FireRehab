# Converting CWD volume to carbon
# Ingrid Farnell
# Feb 5, 2021

# Libraries
library(dplyr)

# Import transect horizontal distances
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

##################################################################

# To convert volume (calculated using the VanWagner formula) to carbon you have to 
# calculate the volume for each decay class and species, then convert that to live biomass then dead biomass

# There are differnet conversion factors for each decay class and species. I will calculate the volume for
# each species and decay class in each plot - then convert that to biomass and carbon - then sum to 
# get the total plot carbon

# Volume equation
# CWD volume (m3/ha) = pi^2/8L  *  sum[D2]          

# Where: 	L = length of total transect (horizontal distance (HD) in m) 
#                                       HD = SD / Square root of [1 + (% slope / 100)2]
# D = diameter of each piece of CWD (cm)

################################################################

# Since the volume is per hectare have to sum diameter to plot, but first square individual diameters
cwd$D2 <- cwd$Diam_cm^2

# Plot summary by Decay class
cwd.plot <- cwd %>%
  group_by(PlotID, Decay_class, Species) %>%
  summarize_at(c("D2"), sum)

# Merge line.plot with cwd.plot
cwd.line.plot <- left_join(cwd.plot, line.plot, by = c("PlotID"))

# Calculate plot volume(m3/ha) for each decay class
cwd.line.plot$volume_ha<- pi^2/(8* cwd.line.plot$hor.dist) * cwd.line.plot$D2 

# Make volume NA = 0
cwd.line.plot$volume_ha[is.na(cwd.line.plot$volume_ha)] <- 0

# export table - to use from this point on
#write.csv(cwd.line.plot, "C:/Users/Farne/Documents/Borealis_Ecological_Services/BVRC_20-37_TreeRecruitmentAfterFire/Data/Exported/CWD_volume_Sp_DC.csv")


# cwdCARBON (T/ha) = volume(m3/ha)
#                     x structural reduction factor # decay class specific (Fraver et al. 2013)
#                     x Absolute density(g/cm3) # species and decay class specific (Harmon et al. 2008)
#                     x CarbonConcentration # species and decay class specific (Harmon et al. 2013)

cwd.line.plot <- as.data.frame(cwd.line.plot)
cwd.line.plot$Decay_class <- as.factor(cwd.line.plot$Decay_class)

# Assign carbon to new column
t <- vector()
for(i in 1:nrow(cwd.line.plot)){
  t[i] <- cwdCarbonFN(volume_ha = cwd.line.plot[i,"volume_ha"], 
                      Decay_class = cwd.line.plot[i,"Decay_class"], 
                      Species = cwd.line.plot[i,"Species"])
}
cwd.line.plot$carbon <- t
cwd.line.plot$carbon[is.na(cwd.line.plot$carbon)] <- 0

# Sum total carbon in each plot
cwd.carbon.plot <- cwd.line.plot %>%
  group_by(PlotID) %>%
  summarize_at(c("volume_ha", "carbon"), sum)
