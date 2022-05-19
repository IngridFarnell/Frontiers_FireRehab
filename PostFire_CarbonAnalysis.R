# Alana Clason, Ingrid Farnell, Erica Lilles
# January 2022

##########################################################################################################
# This script is for cleaning data, calculating carbon, analysing data, and producing figures and tables #
# for Clason, Farnell and Lilles, 2022. Carbon 5 to 60 years after fire: Planting trees does not         #
# compensate for losses in dead wood stores                                                              #
##########################################################################################################

## Load libraries
library(sf)
library(data.table)
library(raster)
library(rsq)
library(ggplot2)
library(lme4)
library(ggsci)
source("./R/CarbonFunctions.R")

## set pathways 
DatPath <- "./Inputs/"

## Import data
#Plot
FR_treatments <- fread(paste0(DatPath,"FR_Treatments.csv")) # has field plot assessed treatments and fire year
#decided to change FR08 to NP because while it was a plantation, it was burned at low severity and not planted after the fire.
#Tree data:
A1trees <- fread(paste0(DatPath,"A1trees.csv"))
B1trees <- fread(paste0(DatPath,"B1trees.csv"))
Regen <- fread(paste0(DatPath,"Regen.csv"))

#Soils data:
Soils <- fread(paste0(DatPath,"Soils.csv"))

#Woody debris:
cwd <- read.csv(paste0(DatPath,"FireRehabData_CWD.csv"),header = T,stringsAsFactors = T)
fwd <- read.csv(paste0(DatPath,"FireRehabData_FWD.csv"),header = T,stringsAsFactors = T)
line <- read.csv(paste0(DatPath,"FireRehabData_TransectDistance.csv"),header = T,stringsAsFactors = T) 

#Plot treatment cleaning
FR_treatments[,ID:=as.factor(ID)][,Planted:=as.factor(Planted)][,CumBurnSev:= as.factor(CumBurnSev)]
FR_treatments$CumBurnSev <- factor(FR_treatments$CumBurnSev, levels=c("L","LM","M","MH","H"))
FR_treatments[,FIRE_YR_cat:=ifelse(FIRE_YEAR<=1979,"1960-1979",
                                   ifelse(FIRE_YEAR>1979 & FIRE_YEAR<=1989,"1980-1989",
                                          ifelse(FIRE_YEAR>1989 & FIRE_YEAR<=1999,"1990-1999",
                                                 ifelse(FIRE_YEAR>1999 & FIRE_YEAR<=2010,"2000-2009",
                                                        "2010+"))))]
FR_treatments[,TimeSinceFire := 2020 - FIRE_YEAR]

#Climate
AnnClim <- fread(paste0(DatPath,"FR_plots_Clim_Normal_1961_1990Y.csv"))
FR_treatments <- merge(FR_treatments, AnnClim[,.(ID1,MAT, MWMT, MCMT, TD, MAP, MSP, AHM, SHM, DD_0, DD5, DD_18,
                                                 DD18, NFFD, bFFP, eFFP, FFP, PAS, EMT, EXT, MAR, Eref, CMD, RH,
                                                 CMI, DD1040)],by.x="ID",by.y="ID1")

hli_rast <- raster("./Inputs/Rasters/DEMhli.tif")
slope_rast <- raster("./Inputs/Rasters/DEMslope.tif")
Study_plots <- read_sf("./Inputs/Shapefiles/FRplots.shp")
Can_dNBR <- raster("./Inputs/Rasters/CanLaBS_dNBR1000_v0.tif")



###################################################################################################### 
############################################ Data Prep ###############################################
######################################################################################################

################################## Tree Carbon ################################## 
# A1 plot = 5.64m radius = 100m2, B1 = 11.28m radius= 400 m2, A2 = 3.99m radius = 50m2
# convert to Mg/ha: Kg/m2 x 10: 
# (Kg C in the tree)/plot are(m2) x (1Mg/1000kg) x (10000m2/ha) == (Kg C in the tree)/plot are(m2) x 10
# also here, we account for the sub-sampling of some plots
A1trees[,Species := as.factor(Species)]
B1trees[,Species := as.factor(Species)]
#A1 trees biomass
t <- vector()
for(i in 1:nrow(A1trees)){
  t[i] <- TreeCarbonFN(Species = A1trees[i,Species], DBH= A1trees[i,DBH], 
                       HT= A1trees[i,Height], Tree_class=A1trees[i,Tree_class])
}
A1trees[,AreaSearchM2 := ifelse(`Sub-plot`=="A1",100, ifelse(`Sub-plot`=="A2",50,400))]
A1trees[Notes=="only 1/4 of A1",AreaSearchM2:=100/4]
A1trees[Notes=="only 1/2 of A1",AreaSearchM2:=100/2]
A1trees[Notes=="only 1/4 of A2",AreaSearchM2:=50/4]
A1trees[,MgPerHa_fac:=10/AreaSearchM2] #calculating the factor to mulitply to get Mg/Ha
A1trees[,CarbonPerHa:=t*MgPerHa_fac]
A1trees[is.na(CarbonPerHa)]#No A1 trees, biomass==0
A1trees[is.na(CarbonPerHa), CarbonPerHa:=0]
A1trees[CarbonPerHa ==0]
#B1 trees biomass
r <- vector()
for(i in 1:nrow(B1trees)){
  r[i] <- TreeCarbonFN(Species = B1trees[i,Species], DBH= B1trees[i,DBH], 
                       HT= B1trees[i,Height], Tree_class=B1trees[i,Tree_class])
}
B1trees[,AreaSearchM2 := ifelse(`Sub-plot`=="A1",100, ifelse(`Sub-plot`=="A2",50,400))]
B1trees[,MgPerHa_fac:=10/AreaSearchM2] #calculating the factor to mulitply to get Mg/Ha
B1trees[,CarbonPerHa:=r*MgPerHa_fac]
B1trees[is.na(CarbonPerHa), CarbonPerHa:=0]

#merge A1 and B1 together
A1B1trees <- rbind(A1trees[,.(PlotID,`Sub-plot`,Species,DBH,Height,Tree_class,AreaSearchM2,CarbonPerHa)], 
                   B1trees[,.(PlotID,`Sub-plot`,Species,DBH,Height,Tree_class,AreaSearchM2,CarbonPerHa)])
A1B1trees[,PHF:=10000/AreaSearchM2] #accounting for smaller search areas
#Live trees by plot
LiveTree_plots <- merge(FR_treatments[,.(ID)],A1B1trees[Tree_class<3, sum(CarbonPerHa),by="PlotID"],
                        by.x="ID", by.y="PlotID", all=TRUE)
LiveTree_plots$V1[is.na(LiveTree_plots$V1)] <-0
setnames(LiveTree_plots,"V1","LiveTreeCperHa")

#Dead trees by plot
DeadTree_plots <- merge(FR_treatments[,.(ID)],A1B1trees[Tree_class>=3, sum(CarbonPerHa),by="PlotID"],
                        by.x="ID", by.y="PlotID", all=TRUE)
DeadTree_plots$V1[is.na(DeadTree_plots$V1)] <-0
setnames(DeadTree_plots,"V1","DeadTreeCperHa")
live_dead_CperHa <- merge(LiveTree_plots,DeadTree_plots,by="ID")

FR_treatments <- merge(FR_treatments,live_dead_CperHa,by="ID")

################################## Regeneration carbon ##################################
Regen[,Height_class:=`Height_class(cm)`][,`Height_class(cm)`:=NULL]
Regen$Diam_est <- 1 # for 31-130 cm height class
Regen$Diam_est[Regen$Height_class == "0-30"] <- 0.1
h <- vector()
for(i in 1:nrow(Regen)){
  h[i] <- RegenCarbonFN_Ung(Species = Regen[i,Species], Diam_est = 0.1, 
                            Height_class = Regen[i,Height_class], Health = Regen[i,`Live/Dead`] )
}
Regen[,AreaSearchM2:=ifelse(`Sub-Plot`=="A1",100,50)]
Regen[,MgPerHa_fac:=10/AreaSearchM2]
Regen[,CarbonPerHa := h*Tally*MgPerHa_fac] 
Regen[CarbonPerHa==0]
Regen[is.na(CarbonPerHa)]
PlotRegen <- Regen[,.(Regen_MGHa = sum(CarbonPerHa)),by=c("PlotID","Live/Dead")]
PlotRegen <- dcast(PlotRegen, PlotID ~ `Live/Dead`, value.var="Regen_MGHa")
setnames(PlotRegen,c("L","D"),c("Regen_L_MGHa","Regen_D_MGHa"))

FR_treatments <- merge(FR_treatments,PlotRegen,by.x="ID",by.y="PlotID")


################################## Root carbon ##################################
FR_treatments[,LiveRootC:=((LiveTreeCperHa*0.222)+(Regen_L_MGHa*0.222)),by="ID"] #root =  0.222 x live biomass
FR_treatments[,DeadRootC:=((DeadTreeCperHa*exp(-0.1044*TimeSinceFire))+ #decaying the tree roots by time since fire
                             (Regen_D_MGHa*0.222)),by="ID"] #root =  0.222 x  dead biomass


################################## Mineral soil carbon ##################################
Soils[,C_pro := MinSoil_C_PC/100]

#calculating bulk density
Soils[,TotalSample_wgt_kg := (MinSoilAll_DryWgt-MinSoil_TinWgt+ #this includes the coarse fragments
                                Root_wgt_inMin+Black_C_wgt_inMin+
                                Litter_Wgt_inMin)/1000]
Soils[,CoarseFrag_wgt_kg := (Coarse_frag_weight_W_tin - MinSoil_TinWgt)/1000]
#Bulk density calculations:
Soils[,OrgFrag_wgt_kg := (Root_wgt_inMin + Black_C_wgt_inMin + Litter_Wgt_inMin)/1000]
Soils[,FineFract_wgt_kg := (MinSoilAll_DryWgt - Coarse_frag_weight_W_tin)/1000 - OrgFrag_wgt_kg]
Soils[,CoarseAndOrgFragM3M3 := (CoarseFrag_wgt_kg/2650 + OrgFrag_wgt_kg/500)/ (BulkDensity_ml/1000000) ]
## Fine fraction BulkDensity ####
#fine fraction bulk density defined in Bulmer and Simpson 2010 (Soil Compaction Reduced the Growth of Lodgepole Pine and Douglas-fi r Seedlings in Raised Beds after Two Growing Seasons) Fine fraction
Soils[,BDcalc_finefract := FineFract_wgt_kg/((BulkDensity_ml/1000000)-
                                               (CoarseFrag_wgt_kg/2650)-(OrgFrag_wgt_kg/500))]#kg/m3
FR_minSOC_finefract_kg_m2 <- Min_SOC(Soc = Soils[,C_pro], BD = Soils[,BDcalc_finefract],
                                     depth = Soils[,BDDepth], CoarseFrags = Soils[,CoarseAndOrgFragM3M3])
## Soil Organic Carbon ####
#using 0.5 biomass to C conversion for litter and roots
#using a 0.75 char to C mass conversion from Donato et al. 2009 (Quantifying char in postfire woody detritus inventories)
Soils[,FR_BlackC_kg := (Black_C_wgt_inMin*0.75)/1000] #removed roots, and added litter to litter
Soils[,FR_BlackC_kg_m2:= FR_BlackC_kg/(BulkDensity_ml/1000000) * BDDepth]
Soils[,FR_SOM:= C_pro/0.47]
Soils[,BDcalc_finefract_g_cm3 := BDcalc_finefract*1000/1e+06] #calculate finefraction

#Test the predicted (Perie and Ouimet) BD vs observed BD
Soils[,BDest_finefract_g_cm3 := -1.977+4.105*FR_SOM - 1.229*log(FR_SOM) - 0.103*log(FR_SOM)^2] #estimate fine fraction

## Estimating missing data in plot 63 
#use estimated fine fraction from Perie and Ouimet
Soils[ID=="FR63",BDcalc_finefract := BDest_finefract_g_cm3*1000]
#back calculate likely volume of hole for 63 to estimate other C volume
Soils[ID=="FR63", BulkDensity_ml := (FineFract_wgt_kg/(BDcalc_finefract + (CoarseFrag_wgt_kg/2650) -
                                                         (OrgFrag_wgt_kg/500)))*1000000]
Soils[ID=="FR63",FR_BlackC_kg_m2 := FR_BlackC_kg/(BulkDensity_ml/1000000) *0.2] #calculate black C Kg/m2
Soils[ID=="FR63",CoarseAndOrgFragM3M3 := (CoarseFrag_wgt_kg/2650 + OrgFrag_wgt_kg/500)/ (BulkDensity_ml/1000000) ]
FR_minSOC_finefract_kg_m2[63] <- Min_SOC(Soc = Soils[ID=="FR63",C_pro], 
                                         BD = Soils[ID=="FR63",BDcalc_finefract],
                                         depth = 0.2, 
                                         CoarseFrags = Soils[ID=="FR63",CoarseAndOrgFragM3M3])
FR_minSOC_kg_m2 <-FR_minSOC_finefract_kg_m2 + Soils[,FR_BlackC_kg_m2] #fine fraction SOC + BlackCarbon C
FR_minSOC_Mg_ha <- (FR_minSOC_kg_m2/1000)*10000 #convert to Mg per hectare
FR_treatments[,MinSoilC_Mgha := FR_minSOC_Mg_ha]

# Add textures
Soils[,CLAY := sum(MinSoil_Clay_PC_hyd,MinSoil_Clay_PC_h202_hyd,na.rm = TRUE), by="ID"]
Soils[,SILT := sum(MinSoil_Silt_PC_hyd, MinSoil_Silt_PC_h202_hyd,na.rm = TRUE), by="ID"]
Soils[,SAND := sum(MinSoil_Sand_PC_hyd,MinSoil_Sand_PC_h202_hyd,na.rm = TRUE), by="ID"]
Soils[,TotalTexture := CLAY + SILT + SAND]

### Add C:N ratio
Soils[,CN := MinSoil_C_PC/MinSoil_N_PC]

#combine this into the plot treatments
FR_treatments <- merge(FR_treatments, Soils[,.(ID,CN,CLAY)],by="ID")
FR_treatments[,CoarseFrag_wgt_kg := Soils[,CoarseFrag_wgt_kg]]
FR_treatments[,CN:=scale(CN, center=TRUE,scale=TRUE)]
FR_treatments[,CoarseFrag_wgt_kg:=scale(CoarseFrag_wgt_kg, center=TRUE,scale=TRUE)]
FR_treatments[,CLAY:=scale(CLAY, center=TRUE,scale=TRUE)]

################################## Litter carbon ##################################

FR_treatments[,Litter_C_g := Soils[,(Litter_DryWgt*(Litter_C_PC/100)) + Litter_Wgt_inMin*0.5]]
FR_treatments[,Litter_C_MgHa := Litter_C_g/16]

################################## Forest floor carbon ##################################
Soils[,ForestFloor_C_g := sum((ForFloor_DryWgt*ForFloor_C_PC/100), Wood_chunks_inFF*0.5,
                              Black_C_wgt_inFF*0.75,na.rm=TRUE), by="ID"]
#tiny amt of FF in FR03 didn't make it to the lab, so just use 0.5
Soils[ID=="FR03",ForestFloor_C_g := ForFloor_DryWgt*0.5+Wood_chunks_inFF*0.5+Black_C_wgt_inFF*0.75 ]
#the only true NA (missing sample) is FR02
Soils[ID=="FR02",ForestFloor_C_g := NA]
FR_treatments[,ForestFloor_C_g := Soils[,ForestFloor_C_g]]
FR_treatments[,ForestFl_MgHa :=ForestFloor_C_g/16] #converting to Mg/Ha

################################## Coarse woody debris carbon ##################################
source("./R/CWD_carbon.R")
cwd.carbon.plot <- as.data.table(cwd.carbon.plot)
plot.desc<- FR_treatments
# First column imported weird - rename it
names(plot.desc)[names(plot.desc) == "ID"] <- "PlotID"
# Merge dataframes via PlotID
cwd.desc <- merge(cwd.carbon.plot, plot.desc, by = "PlotID")
cwd.desc$FIRE_YR_cat <- as.factor(cwd.desc$FIRE_YR_cat)
FR_treatments[,CWD_C:=cwd.carbon.plot$carbon]

################################## Fine woody debris carbon ##################################
source("./R/FWD_carbon.R")
fwd.carbon.plot <- as.data.table(fwd.carbon.plot)
plot.desc<- FR_treatments
# First column imported weird - rename it
names(plot.desc)[names(plot.desc) == "ID"] <- "PlotID"
# Merge dataframes via PlotID
fwd.desc <- merge(fwd.carbon.plot, plot.desc, by = "PlotID") #not working now
fwd.desc <- as.data.frame(fwd.desc)
fwd.desc$FIRE_YR_cat <- as.factor(fwd.desc$FIRE_YR_cat)
FR_treatments[,FWD_C:=fwd.carbon.plot$carbon]

################################## Carbon pool summaries ##################################
CarbonPlots <- melt(FR_treatments, id.vars = c("ID", "Date", "HarvHist","Planted","CumBurnSev","CumBurnSevCat",
                                               "Aspect","Slope_PC", "SlopePos","FIRE_NUMBE", "FIRE_NAME", 
                                               "FIRE_YEAR", "area","FIRE_YR_cat"),
                    measure.vars = c("LiveTreeCperHa","Regen_L_MGHa","Regen_D_MGHa","DeadTreeCperHa",
                                     "MinSoilC_Mgha","Litter_C_MgHa",
                                     "ForestFl_MgHa","CWD_C","FWD_C","LiveRootC"),
                    variable.name = "CarbonSource",
                    value.name = "CarbonMgHa")

FR_treatments <- merge(FR_treatments,CarbonPlots[,.(TotalCarbon = sum(na.omit(CarbonMgHa))), by="ID"],by="ID")
FR_treatments <- merge(FR_treatments,CarbonPlots[CarbonSource!="MinSoilC_Mgha",
                                                 .(TotalCarbon_NotMin = sum(na.omit(CarbonMgHa))), by="ID"],by="ID")

FR_treatments <- merge(FR_treatments,CarbonPlots[CarbonSource=="DeadTreeCperHa"|
                                                   CarbonSource=="Regen_D_MGHa"|
                                                   CarbonSource=="CWD_C"|
                                                   CarbonSource=="FWD_C",
                                                 .(DeadCarbon = sum(na.omit(CarbonMgHa))),
                                                 by="ID"],by="ID")
FR_treatments <- merge(FR_treatments,CarbonPlots[CarbonSource=="LiveTreeCperHa"|
                                                   CarbonSource=="Regen_L_MGHa"|
                                                   CarbonSource=="LiveRootC",
                                                 .(LiveCarbon = sum(na.omit(CarbonMgHa))),
                                                 by="ID"],by="ID")
FR_treatments <- merge(FR_treatments,CarbonPlots[CarbonSource=="ForestFl_MgHa"|
                                                   CarbonSource=="Litter_C_MgHa",
                                                 .(ForestFloorTotC = sum(na.omit(CarbonMgHa))),
                                                 by="ID"],by="ID")

FR_tables <- melt(FR_treatments,id.vars = c("ID","Planted","TimeSinceFire"),
                  measure.vars = c("LiveTreeCperHa","Regen_L_MGHa","LiveRootC","DeadTreeCperHa","Regen_D_MGHa",
                                   "CWD_C","FWD_C"),
                  variable.name = "CarbonPool",
                  value.name = "CarbonMgHa")

######################################## Additional Predictors #####################################

################################## Climate ##################################
FR_treatments[,DD5 := scale(DD5, center=TRUE, scale=TRUE)]
FR_treatments[,MAP := scale(MAP, center=TRUE, scale=TRUE)]

################################## Topography ##################################
Plot_hli <- raster::extract(hli_rast,Study_plots,method="simple")
Plot_slope <- raster::extract(slope_rast,Study_plots,method="simple")

FR_treatments[,hli := scale(Plot_hli,center=TRUE,scale=TRUE)]
FR_treatments[,slope_ra := scale(Plot_slope,center=TRUE,scale=TRUE)] 

################################## Fire severity ##################################
Study_plots_tr <- st_transform(Study_plots,crs=crs(Can_dNBR))
Plot_dNBR <- raster::extract(Can_dNBR,Study_plots_tr,method="simple")
#how many are NA after 1985:
FR_treatments[,dNBR:=Plot_dNBR]
FR_treatments[,dNBR := ifelse(!is.na(Plot_dNBR),Plot_dNBR,ifelse(FIRE_YEAR>1985,1,Plot_dNBR))]
FR_treatments[,dNBR_nsc := dNBR]
FR_treatments[,dNBR_Class := ifelse(dNBR_nsc<=100,"Unburned",
                                    ifelse(dNBR_nsc>=100 & dNBR_nsc <270,"L",
                                           ifelse(dNBR_nsc>=270 & dNBR_nsc<440, "LM",
                                                  ifelse(dNBR_nsc>=440 & dNBR_nsc<660, "MH",
                                                         "High"))))]
FR_treatments[,dNBR := scale(dNBR, center=TRUE,scale=TRUE)]



##################### Table 2 ###################################################
FR_tables[,TSF_cat:=ifelse(TimeSinceFire >50,55,
                           ifelse(TimeSinceFire >40,45,
                                  ifelse(TimeSinceFire > 30, 35,
                                         ifelse(TimeSinceFire>20,25,
                                                ifelse(TimeSinceFire>10,15,
                                                       5)))))]
FR_pool_summary <- FR_tables[,.(mn_CarbonMgHa=mean(CarbonMgHa),sd_CarbonMgHa=sd(CarbonMgHa)),
                             by=c("CarbonPool","Planted","TSF_cat")]
###################################################################################################### 
############################################ Analyses ################################################
######################################################################################################
Mod <- list()
Mod[[1]] <- lmer(TotalCarbon ~ Planted-1 + TimeSinceFire + DD5 + MAP + CN + CoarseFrag_wgt_kg +
                   CLAY+hli + slope_ra +(1|FIRE_NAME),data=FR_treatments)
Mod[[2]] <- lmer(MinSoilC_Mgha ~ Planted-1 + TimeSinceFire  + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli +
                   slope_ra +(1|FIRE_NAME), data=FR_treatments)
Mod[[3]] <- lmer(LiveCarbon ~ Planted-1 + TimeSinceFire + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli + 
                   slope_ra +(1|FIRE_NAME), data=FR_treatments)
Mod[[4]] <- lmer(DeadCarbon ~ Planted-1 + TimeSinceFire  + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli +
                   slope_ra+(1|FIRE_NAME),data=FR_treatments)
Mod[[5]] <- lmer(ForestFloorTotC ~ Planted-1 + TimeSinceFire  + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli +
                   slope_ra +(1|FIRE_NAME),data=FR_treatments)

ModNoPlant <- list()
ModNoPlant[[1]] <- lmer(TotalCarbon ~  TimeSinceFire + DD5 + MAP + CN + CoarseFrag_wgt_kg +
                          CLAY+hli + slope_ra +(1|FIRE_NAME),data=FR_treatments)
ModNoPlant[[2]] <- lmer(MinSoilC_Mgha ~ TimeSinceFire  + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli +
                          slope_ra +(1|FIRE_NAME), data=FR_treatments)
ModNoPlant[[3]] <- lmer(LiveCarbon ~ TimeSinceFire + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli + 
                          slope_ra +(1|FIRE_NAME), data=FR_treatments)
ModNoPlant[[4]] <- lmer(DeadCarbon ~ TimeSinceFire  + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli +
                          slope_ra+(1|FIRE_NAME),data=FR_treatments)
ModNoPlant[[5]] <- lmer(ForestFloorTotC ~ TimeSinceFire  + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli +
                          slope_ra +(1|FIRE_NAME),data=FR_treatments)
#Total, Mineral soil, live, dead, forestfloor
lapply(Mod,AIC)
lapply(ModNoPlant,AIC)
lapply(Mod,rsq.lmm)
lapply(Mod,summary)
lapply(ModNoPlant,summary)
lapply(ModNoPlant,rsq.lmm)

########## With fire severity ############
Mod <- list()
Mod[[1]] <- lmer(TotalCarbon ~ Planted-1 + TimeSinceFire  + dNBR+ DD5 + MAP + CN + CoarseFrag_wgt_kg +
                   CLAY+hli + slope_ra +(1|FIRE_NAME),data=FR_treatments)
Mod[[2]] <- lmer(MinSoilC_Mgha ~ Planted-1 + TimeSinceFire  + dNBR + DD5 + MAP + CN + CoarseFrag_wgt_kg + 
                   CLAY+hli +slope_ra +(1|FIRE_NAME), data=FR_treatments)
Mod[[3]] <- lmer(LiveCarbon ~ Planted-1 + TimeSinceFire  + dNBR+ DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli + 
                   slope_ra +(1|FIRE_NAME), data=FR_treatments)
Mod[[4]] <- lmer(DeadCarbon ~ Planted-1 + TimeSinceFire  + dNBR + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli +
                   slope_ra+(1|FIRE_NAME),data=FR_treatments)
Mod[[5]] <- lmer(ForestFloorTotC ~ Planted-1 + TimeSinceFire + dNBR  + DD5 + MAP + CN + CoarseFrag_wgt_kg +
                   CLAY+hli + slope_ra +(1|FIRE_NAME),data=FR_treatments)

ModNoPlant <- list()
ModNoPlant[[1]] <- lmer(TotalCarbon ~  TimeSinceFire  + dNBR+ DD5 + MAP + CN + CoarseFrag_wgt_kg +
                          CLAY+hli + slope_ra +(1|FIRE_NAME),data=FR_treatments)
ModNoPlant[[2]] <- lmer(MinSoilC_Mgha ~ TimeSinceFire + dNBR + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli +
                          slope_ra +(1|FIRE_NAME), data=FR_treatments)
ModNoPlant[[3]] <- lmer(LiveCarbon ~ TimeSinceFire + dNBR + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli + 
                          slope_ra +(1|FIRE_NAME), data=FR_treatments)
ModNoPlant[[4]] <- lmer(DeadCarbon ~ TimeSinceFire + dNBR  + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli +
                          slope_ra+(1|FIRE_NAME),data=FR_treatments)
ModNoPlant[[5]] <- lmer(ForestFloorTotC ~ TimeSinceFire + dNBR + DD5 + MAP + CN + CoarseFrag_wgt_kg + CLAY+hli +
                          slope_ra +(1|FIRE_NAME),data=FR_treatments)

lapply(Mod,AIC)
lapply(ModNoPlant,AIC)
lapply(Mod,summary)
lapply(Mod,rsq.lmm)
lapply(ModNoPlant,summary)
lapply(ModNoPlant,rsq.lmm)


#####test the difference in variance between planted and not planted:
source("./R/DensityFunctions.R")
PlotRegen <- RegenDensity(Regen)
PlotRegenSum <- PlotRegen[,sum(SPH),by=c("PlotID,Species")] #combine 0-0.3 and 0.3-1.3
setnames(PlotRegenSum,"V1","SPH")
FR_Regen <- merge(FR_treatments,PlotRegenSum,by.x="ID",by.y="PlotID")
FR_Regen[,TSF_cat:=as.factor(ifelse(TimeSinceFire >50,55,
                                    ifelse(TimeSinceFire >40,45,
                                           ifelse(TimeSinceFire > 30, 35,
                                                  ifelse(TimeSinceFire>20,25,
                                                         ifelse(TimeSinceFire>10,15,
                                                                5))))))]
#Young stands:
regen_summary <- FR_Regen[TimeSinceFire<30, sum(SPH), by=c("ID","Planted","TimeSinceFire","TSF_cat")]
setnames(regen_summary,"V1","regenStems")
car::leveneTest(regenStems ~ Planted, data = regen_summary)
t.test(regen_summary[Planted=="P"]$regenStems,regen_summary[Planted=="NP"]$regenStems)

FR_regen_test <- merge(FR_treatments, dcast(PlotRegenSum, PlotID~Species,
                                            value.var = "SPH"),by.x="ID",by.y="PlotID")
#using the kruskal wallis test because of large number of ties that make te Wilcox rank test difficult
wilcox.test(Pl~Planted, data=FR_regen_test, paired=FALSE)
wilcox.test(Bl~Planted, data=FR_regen_test, paired=FALSE)
wilcox.test(Ac~Planted, data=FR_regen_test, paired=FALSE)
wilcox.test(Sx~Planted, data=FR_regen_test, paired=FALSE)


###################################################################################################### 
############################################ Figures #################################################
######################################################################################################

FR_pools <- FR_treatments[,.(ID,Planted,TimeSinceFire,TotalCarbon,LiveCarbon,DeadCarbon,ForestFloorTotC,MinSoilC_Mgha)]
FR_pools[,Planted:=ifelse(Planted=="NP","Not planted", "Planted")]
FR_pools[,TSF_cat:=ifelse(TimeSinceFire >50,55,
                          ifelse(TimeSinceFire >40,45,
                                 ifelse(TimeSinceFire > 30, 35,
                                        ifelse(TimeSinceFire>20,25,
                                               ifelse(TimeSinceFire>10,15,
                                                      5)))))]
FR_pools_m <- melt(FR_pools,id.vars = c("ID","Planted","TSF_cat"),
                   measure.vars = c("LiveCarbon","DeadCarbon","ForestFloorTotC","MinSoilC_Mgha"),
                   variable.name = "CarbonPool",
                   value.name = "CarbonMgHa")
FR_means <- FR_pools_m[,.(mn_CarbonMgHa=mean(CarbonMgHa),sd_CarbonMgHa=sd(CarbonMgHa)),
                       by=c("CarbonPool","Planted","TSF_cat")]

##Figure 2a
FR_pools <- FR_treatments[,.(ID,Planted,TimeSinceFire,TotalCarbon,LiveCarbon,DeadCarbon,ForestFloorTotC,MinSoilC_Mgha)]
FR_pools[,Treatment:=ifelse(Planted=="NP","Not planted", "Planted")]
FR_pools_m <- melt(FR_pools,id.vars = c("ID","Treatment","TimeSinceFire"),
                   measure.vars = c("TotalCarbon","LiveCarbon","DeadCarbon","ForestFloorTotC","MinSoilC_Mgha"),
                   variable.name = "CarbonPool",
                   value.name = "CarbonMgHa")
supp.labs <- c("Total","Live","Dead","Forest floor","Mineral soil")
names(supp.labs) <- c("TotalCarbon","LiveCarbon","DeadCarbon","ForestFloorTotC","MinSoilC_Mgha")

ggplot(FR_pools_m,aes(x=TimeSinceFire, y=CarbonMgHa))+
  geom_point(size=2,aes(colour=Treatment))+
  geom_smooth(method="lm",aes(colour=Treatment,fill=Treatment))+
  scale_color_npg()+
  scale_fill_npg()+
  theme_minimal() +
  ylab(expression("Carbon Mg" ~ ha^-1))+
  xlab("Time since fire (years)")+
  #theme(legend.position = "bottom")+
  theme(legend.justification = c(1, 0), legend.position = c(1, 0))+
  facet_wrap("CarbonPool",labeller=labeller(CarbonPool=supp.labs))+
  theme(strip.text.x = element_text(face="bold"))
  #theme(strip.text.x = element_text(face="bold"), text=element_text(size=21)) #for presentations
  ggsave(filename = "Fig2a.jpg",path = "./Outputs/", device='jpeg', dpi=300, bg="white")

## Figure 2b
ggplot(FR_means)+
  geom_area(aes(x=TSF_cat,y=mn_CarbonMgHa,fill=CarbonPool))+
  facet_grid(rows="Planted")+
  scale_fill_npg(name="Carbon Pool",labels=c("Live","Dead","Forest floor","Mineral soil"))+
  scale_color_npg()+
  ylab(expression("Carbon Mg" ~ ha^-1))+
  xlab("Time since fire (years)")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank())
ggsave(filename = "Fig2b.jpg",path = "./Outputs/", device='jpeg', dpi=300, bg="white")

####### Figure 3
Dead_part <- FR_treatments[,.(ID,Planted,TimeSinceFire, dNBR_Class,DeadTreeCperHa, Regen_D_MGHa, CWD_C, FWD_C)]
Dead_part2 <-  melt(Dead_part,id.vars = c("ID","Planted","TimeSinceFire","dNBR_Class"),
                    measure.vars = c("DeadTreeCperHa","Regen_D_MGHa","CWD_C","FWD_C"),
                    variable.name = "DeadCarbonPool",
                    value.name = "CarbonMgHa")
Dead_part2[,TSF_cat:=ifelse(TimeSinceFire >50,55,
                            ifelse(TimeSinceFire >40,45,
                                   ifelse(TimeSinceFire > 30, 35,
                                          ifelse(TimeSinceFire>20,25,
                                                 ifelse(TimeSinceFire>10,15,
                                                        5)))))]
Dead_means <- Dead_part2[,.(mn_CarbonMgHa=mean(CarbonMgHa)),by=c("DeadCarbonPool","TSF_cat","Planted")]
supp.labs <- c("Not planted","Planted")
names(supp.labs) <- c("NP","P")

### 3a
ggplot(Dead_means)+
  geom_area(aes(x=TSF_cat,y=mn_CarbonMgHa,fill=DeadCarbonPool))+
  facet_grid(rows="Planted",labeller=labeller(Planted=supp.labs))+
  scale_fill_npg(name="Dead Carbon Pool",labels=c("Snags","Dead seedlings","CWD","FWD", "Dead roots"))+
  scale_color_npg()+
  ylab(expression("Carbon Mg" ~ ha^-1))+
  xlab("Time since fire (years)")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank())
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
   #     strip.background = element_blank(),strip.text.x = element_text(face="bold"), 
    #    text=element_text(size=21)) #for presentations

ggsave(filename = "Fig3a.jpg",path = "./Outputs/", device='jpeg', dpi=300, bg="white")

### 3b
Dead_means <- Dead_part2[!is.na(dNBR_Class),.(mn_CarbonMgHa=mean(CarbonMgHa)),
                         by=c("DeadCarbonPool","dNBR_Class","Planted")]
Dead_means[,dNBR_Class:=as.factor(dNBR_Class)]
Dead_means[,dNBR_Class:=ifelse(dNBR_Class=="L","Low",
                               ifelse(dNBR_Class=="LM","Low-Moderate",
                                      ifelse(dNBR_Class=="MH","Moderate-High",
                                             ifelse(dNBR_Class=="Unburned","Unburned",
                                                    "High"))))]
Dead_means[,dNBR_Class:=factor(dNBR_Class, levels = c("Unburned","Low","Low-Moderate",
                                                      "Moderate-High","High"))]

supp.labs <- c("Not planted","Planted")
names(supp.labs) <- c("NP","P")
ggplot(Dead_means)+
  geom_col(aes(x=dNBR_Class,y=mn_CarbonMgHa,fill=DeadCarbonPool))+
  facet_grid(rows="Planted",labeller=labeller(Planted=supp.labs))+
  scale_fill_npg(name="Dead Carbon Pool",labels=c("Snags","Dead seedlings","CWD","FWD", "Dead roots"))+
  scale_color_npg()+
  ylab(expression("Carbon Mg" ~ ha^-1))+
  xlab("Burn severity")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank())
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
   #     strip.background = element_blank(),strip.text.x = element_text(face="bold"),
    #    text=element_text(size=21)) #for presentations

ggsave(filename = "Fig3b.jpg",path = "./Outputs/", device='jpeg', dpi=300, bg="white")

## Figure 5
FR_Regen[,TSF_cat:=ifelse(TimeSinceFire<=10,"\u2264 10",ifelse(TimeSinceFire<=20,"11-20",
                                                         ifelse(TimeSinceFire<=30,"21-30",
                                                                ifelse(TimeSinceFire<=40,"31-40",
                                                                       ifelse(TimeSinceFire<=5,
                                                                              "41-50","51-60")))))]

FR_tables[,TSF_cat:=ifelse(TimeSinceFire >50,55,
                           ifelse(TimeSinceFire >40,45,
                                  ifelse(TimeSinceFire > 30, 35,
                                         ifelse(TimeSinceFire>20,25,
                                                ifelse(TimeSinceFire>10,15,
                                                       5)))))]
supp.labs <- c("Not planted","Planted")
names(supp.labs) <- c("NP","P")
ggplot(FR_Regen[SPH>0])+
  geom_boxplot(aes(x=TSF_cat,y=SPH, fill=Planted))+
  scale_fill_npg()+
  scale_y_log10()+
  ylab(expression("Density (seedlings" ~ ha^-1 ~ "log scale)"))+
  xlab("Time since fire (years) category")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank())+
  facet_wrap(~Planted,labeller=labeller(Planted=supp.labs))+
  theme(legend.position = "none")+
  theme(strip.text.x = element_text(face="bold"))
ggsave(filename = "Fig5.jpg",path = "./Outputs/", device='jpeg', dpi=300, bg="white")


#### Figure 6
supp.labs <- c("Not planted","Planted")
names(supp.labs) <- c("NP","P")
FR_Regen[,Species:= ifelse(Species=="Lw","Western larch",
                           ifelse(Species=="Ac","Cottonwood",
                                  ifelse(Species=="At","Trembling aspen",
                                         ifelse(Species=="Bl","Subalpine fir",
                                                ifelse(Species =="Sx","Hybrid spruce",
                                                       ifelse(Species=="Pl","Lodgepole pine",
                                                              ifelse(Species=="Fd","Douglas-fir",
                                                                     ifelse(Species=="Ep","Paper birch",Species))))))))]
FR_Regen[,ConDec:=ifelse(Species=="Larch"|Species=="Fir"|Species=="Spruce"|Species=="Pine"|
                           Species=="DouglasFir","Conifer","Deciduous")]

ggplot(FR_Regen[SPH>0],aes(x=TimeSinceFire, y=SPH, color=Planted))+
  geom_point(size=4, position="jitter", alpha=0.5)+
  scale_y_log10()+
  scale_color_npg()+
  ylab(expression("Seedlings" ~ ha^-1 ~ "(log scale)"))+
  xlab("Time since fire")+
  facet_wrap(~Species,ncol=3)+
  #theme(panel.border=element_rect(colour="black",size=1))+
  theme_bw() +
  theme(legend.position = "bottom")+
  theme(strip.text.x = element_text(face="bold"))
ggsave(filename = "Fig6.jpg",path = "./Outputs/", device='jpeg', dpi=300, bg="white")

###################################################################################################### 
#################################100 years since fire ################################################
######################################################################################################
Trees_SO <- fread("./Inputs/SORTIEruns/SE_Trees_SO_1_75.csv") #live trees

Trees_SO_live <- Trees_SO[!is.na(CarbonPerHa)]
rm(Trees_SO)
gc()
FR_Sor_trees <- Trees_SO_live[,.(LiveCperHa = sum(CarbonPerHa)),by=c("plotID","timestep")]
plotID <- c(paste0("FR0",seq(1,9,by=1)),paste0("FR",seq(10,75,by=1)))
FR_Sor_trees <- merge(FR_Sor_trees,FR_treatments[ID %in% plotID,.(ID,TimeSinceFire,Planted)],
                      by.x="plotID",by.y="ID")
FR_Sor_trees[,TSF:= timestep+TimeSinceFire]
FR_Sor_trees <- melt(FR_Sor_trees, id.vars = c("plotID", "TSF","Planted"),
                     measure.vars = c("LiveCperHa"),
                     variable.name = "CarbonPerHa",
                     value.name = "CarbonMgHa")

library(likelihood)
#Live Carbon functions 
Models <- list()
Models[[1]] <- function(parA,parB,TimeSinceFire){parA+(parB*TimeSinceFire)}   #linear function
Models[[2]] <- function(parA,parB,TimeSinceFire){(parA*TimeSinceFire)/(parB+TimeSinceFire)} #Michel-Menton function
Models[[3]] <- function(parA,parB,TimeSinceFire){parA*(1-exp(-parB^TimeSinceFire))}#Monomolecular function
Models[[4]] <- function(parA,parB,TimeSinceFire){(parA*TimeSinceFire)^2/(parB^2+TimeSinceFire^2)} #Holling type 3
Models[[5]] <- function(parA,parB,TimeSinceFire){parA+parB*log(TimeSinceFire)} #log model
Models[[6]] <- function(parA,parB,parC,TimeSinceFire){parA+(parB*TimeSinceFire)+(parC*TimeSinceFire)^2} #polynomial
Models[[7]] <- function(parA, parB, parC, TimeSinceFire){parA*(1 - exp(-1* parB * TimeSinceFire))^parC} #CR

var <- list(TimeSinceFire = "TSF",x="CarbonMgHa",mean="predicted",log=TRUE)
Planted_fits <- list()
NotPlanted_fits <- list()
for(j in 1:length(Models)){
  if(j==1|j==2|j==3|j==4|j==5){
    par <- list(parA = 0.1, parB = 0.1,sd=0.5)
    par_lo <- list(parA = -1000, parB = -1000, sd=0)
    par_hi <- list (parA = 4000, parB = 1000, sd=100)
  }else if(j==6){
    par <- list(parA = 0.1, parB = 0.1,parC=0.1, sd=0.5)
    par_lo <- list(parA = -1000, parB = -1000, parC= -1000, sd=0)
    par_hi <- list (parA = 4000, parB = 1000, parC = 1000, sd=100)
  }else{
    par <- list(parA = 0.1, parB = 0.1,parC=0.1, sd=0.5)
    par_lo <- list(parA = 0, parB = -10, parC=-10, sd=0)
    par_hi <- list (parA = 500, parB = 10, parC = 10, sd=100)
  }
  Planted_fits[[j]] <- anneal(model = Models[[j]], par = par, var = var, source_data = FR_Sor_trees[Planted=="P"],
                              par_lo=par_lo, par_hi=par_hi, pdf = dnorm, dep_var="CarbonMgHa", initial_temp = 5, 
                              temp_red = 0.85, max_iter=50000, hessian = FALSE)
  
  NotPlanted_fits[[j]] <- anneal(model = Models[[j]], par = par, var = var, source_data = FR_Sor_trees[Planted=="NP"],
                                 par_lo=par_lo, par_hi=par_hi, pdf = dnorm, dep_var="CarbonMgHa", initial_temp = 5, 
                                 temp_red = 0.85, max_iter=5000, hessian = FALSE)
}

AIC_table <- as.data.table(cbind(Model=c(paste0("Model",seq(1,length(Models)))),
                                 Planted=0,NotPlanted=0,Planted=0,NotPlanted=0,Planted=0,
                                 NotPlanted=0,Planted=0,NotPlanted=0))
for(ii in 1:length(Models)){
  if(ii <6){
    AIC_table[ii,2] <- round(Planted_fits[[ii]]$aic,0)
    AIC_table[ii,3] <- round(NotPlanted_fits[[ii]]$aic,0)
    AIC_table[ii,4] <- round(Planted_fits[[ii]]$best_pars$parA,3)
    AIC_table[ii,5] <- round(NotPlanted_fits[[ii]]$best_pars$parA,3)
    AIC_table[ii,6] <- round(Planted_fits[[ii]]$best_pars$parB,3)
    AIC_table[ii,7] <- round(NotPlanted_fits[[ii]]$best_pars$parB,3)
  }else{
    AIC_table[ii,2] <- round(Planted_fits[[ii]]$aic,0)
    AIC_table[ii,3] <- round(NotPlanted_fits[[ii]]$aic,0)
    AIC_table[ii,4] <- round(Planted_fits[[ii]]$best_pars$parA,3)
    AIC_table[ii,5] <- round(NotPlanted_fits[[ii]]$best_pars$parA,3)
    AIC_table[ii,6] <- round(Planted_fits[[ii]]$best_pars$parB,3)
    AIC_table[ii,7] <- round(NotPlanted_fits[[ii]]$best_pars$parB,3)
    AIC_table[ii,8] <- round(Planted_fits[[ii]]$best_pars$parC,3)
    AIC_table[ii,9] <- round(NotPlanted_fits[[ii]]$best_pars$parC,3)
  }
}
AIC_table
NotPlanted_fits[[1]]$R2
NotPlanted_fits[[4]]$R2
NotPlanted_fits[[5]]$R2
NotPlanted_fits[[6]]$R2
NotPlanted_fits[[7]]$R2
Planted_fits[[1]]$R2
Planted_fits[[4]]$R2
Planted_fits[[5]]$R2
Planted_fits[[6]]$R2
Planted_fits[[7]]$R2

supp.labs <- c("Not planted","Planted")
names(supp.labs) <- c("NP","P")

a <- ggplot(FR_Sor_trees[Planted=="NP"], aes(x=TSF,y=CarbonMgHa))+
  geom_point(colour="#E64B35B2", alpha = 0.1)+
  geom_line(aes(y=NotPlanted_fits[[1]]$best_pars$parA +
                  (NotPlanted_fits[[1]]$best_pars$parB*TSF)),colour="black",size=1, linetype =2)+
  geom_line(aes(y=(NotPlanted_fits[[7]]$best_pars$parA*
                     (1-exp(-1*NotPlanted_fits[[7]]$best_pars$parB*TSF))^NotPlanted_fits[[7]]$best_pars$parC)),
            colour="black",size=1, linetype=1)+
  scale_color_npg()+
  scale_fill_npg()+
  theme_minimal() +
  ylab(expression("Carbon Mg" ~ ha^-1))+
  xlab("Time since fire (years)")+
  ylim(0,200)+
  theme(strip.text.x = element_text(face="bold"))

b <- ggplot(FR_Sor_trees[Planted=="P"], aes(x=TSF,y=CarbonMgHa))+
  geom_point(colour="#4DBBD5B2", alpha=0.1)+
  geom_line(aes(y=Planted_fits[[1]]$best_pars$parA +
                  (Planted_fits[[1]]$best_pars$parB*TSF)),colour="black",size=1, linetype=2)+
  geom_line(aes(y=(Planted_fits[[7]]$best_pars$parA*
                     (1-exp(-1*Planted_fits[[7]]$best_pars$parB*TSF))^Planted_fits[[7]]$best_pars$parC)),
            colour="black",size=1, linetype=1)+
  scale_color_npg()+
  scale_fill_npg()+
  theme_minimal() +
  ylab(expression("Carbon Mg" ~ ha^-1))+
  xlab("Time since fire (years)")+
  ylim(0,200)+
  theme(strip.text.x = element_text(face="bold"))
ggarrange(a,b)
ggsave(filename = "Fig4.jpg",path = "./Outputs/", device='jpeg', dpi=300, bg="white")


###################################################################################################### 
############################################ Map #####################################################
######################################################################################################
library(sf)
library(raster)
library(data.table)

## set pathways 
DatPath <- "./Inputs/"

## Import data
studyFires <- read_sf(paste0(DatPath,"Shapefiles/StudyFires.shp"))
HistoricFires <- read_sf(paste0(DatPath,"Shapefiles/HistoricFires_StudyArea.shp"))
StudyArea <- read_sf(paste0(DatPath,"Shapefiles/StudyArea.shp"))

## Fires since 1960
fires_dt <- as.data.table(HistoricFires)
fires_since60 <- fires_dt[FIRE_YEAR>1959 & FIRE_YEAR<2016]
# area
sum(st_area(st_as_sf(fires_since60)))/10000 # 584,904 ha
nrow(fires_since60)
write_sf(fires_since60,paste0(DatPath,"Shapefiles/AllFires_1960_2015.shp"))




