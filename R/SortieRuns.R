# This creates, runs and calculates carbon for SORTIE runs (100 years and 250 years)
library(data.table)
library(xml2)
library(stringr)

source("../r-SORTIE/ParseXML.R")
source("../r-SORTIE/Functions.R")
source("../r-SORTIE/ReplaceInfo.R")
source("../r-SORTIE/SORTIE-HelperFunctions.R")
source("./R/CarbonFunctions.R")
source("./R/DensityFunctions.R")

####################################################################################
################## Create the initial conditions ##################################
####################################################################################
##data
A1trees <- fread("./Inputs/A1trees.csv")
B1trees <- fread("./Inputs/B1trees.csv")
Regen <- fread("./Inputs/Regen.csv")
FR_treatments <- fread("./Inputs/FR_treatments.csv") # has field plot assessed treatments and fire year
FR_treatments[,TimeSinceFire := 2020 - FIRE_YEAR] # Calculate the time since fire when measurements happened

#create plot starting conditions
InitRootName <- "Init.Dens_"
InitTreesPath <- "./Inputs/SORTIEInits/"

dbhclassSize <- 2 #define diameter class size and number of classes
diamClasses <- seq(0,40, by=dbhclassSize)

#Define the run type - Options are "FireReturn" (100 years) or "OldGrowth" (250 years)
Run_type <- "OldGrowth"

#Get trees ready for SORTIE
PlotTrees <- TreeDensity(A1trees,B1trees)
PlotTrees[,SO_sp:=ifelse(Species=="Lw","Western_Larch", #grow larch as spruce
                         ifelse(Species=="Fd","Douglas_Fir", #grow Doug fir as spruce
                                ifelse(Species=="Bl","Subalpine_Fir",
                                       ifelse(Species =="Sx","Interior_Spruce",
                                              ifelse(Species=="Pl","Lodgepole_Pine",
                                                     ifelse(Species=="At","Trembling_Aspen",
                                                            ifelse(Species=="Ac","Black_Cottonwood",
                                                                   ifelse(Species=="Ep","Paper_Birch",Species))))))))]

PlotTrees[,Init.Values := as.factor(paste0(InitRootName, formatC(as.numeric(DBH_bin),
                                                                 width=3,format='f',digits=1)))]
#Get regen ready for SORTIE
Regen[,Height_class:=`Height_class(cm)`][,`Height_class(cm)`:=NULL]
PlotRegen <- RegenDensity(Regen) #run RegenDensity function
PlotRegen[,SO_sp:=ifelse(Species=="Lw","Western_Larch", #grow larch as spruce
                         ifelse(Species=="Fd","Douglas_Fir", #grow Doug fir as spruce
                                ifelse(Species=="Bl","Subalpine_Fir",
                                       ifelse(Species =="Sx","Interior_Spruce",
                                              ifelse(Species=="Pl","Lodgepole_Pine",
                                                     ifelse(Species=="At","Trembling_Aspen",
                                                            ifelse(Species=="Ac","Black_Cottonwood",
                                                                   ifelse(Species=="Ep","Paper_Birch",Species))))))))]
PlotRegen[,Init.Values := as.factor(paste0("Init_Sdl_Hgt_",SdlHgt))]

#Create the initial conditions csvs with regen, trees, and plot run information
PlotList <- c(paste0("FR0",seq(1,9,by=1)),paste0("FR",seq(10,75,by=1)))
BlankSortie <- data.table(Init.Values = as.factor(paste0(InitRootName, formatC(as.numeric(diamClasses[diamClasses>0]),
                                                                               width=3,format='f',digits=1))),
                          "Interior_Spruce"=0,"Lodgepole_Pine"=0,"Subalpine_Fir"=0,"Trembling_Aspen"=0,
                          "Western_Larch"=0,"Douglas_Fir"=0,"Black_Cottonwood"=0,"Paper_Birch"=0)

for(ii in 1:length(PlotList)){
  #Trees
  if(nrow(PlotTrees[PlotID==PlotList[ii]])!=0){
    a <- dcast(PlotTrees[PlotID==PlotList[ii],.(SO_sp,SPH,Init.Values)],Init.Values ~ SO_sp, value.var = "SPH")
    b <- rbind(BlankSortie[!Init.Values %in% a$Init.Values],a,fill=TRUE)
    for (i in seq_along(b)) set(b, i=which(is.na(b[[i]])), j=i, value=0)
  } else {
    b <- BlankSortie
  }
  c <- rbind(data.table(Init.Values="na"),b,fill=TRUE)
  #Regen
  if(nrow(PlotRegen[PlotID==PlotList[ii]])!=0){
    d <- dcast(PlotRegen[PlotID==PlotList[ii],.(SO_sp,SPH,Init.Values)],Init.Values ~ SO_sp, value.var = "SPH")
    #because regen already has 0s for every species, don't need to merge with blank
  }else{
    d <- data.table(Init.Values = c("Init_Sdl_Hgt_1","Init_Sdl_Hgt_2"),
                    "Interior_Spruce"=0,"Lodgepole_Pine"=0,"Subalpine_Fir"=0,"Trembling_Aspen"=0,
                    "Western_Larch"=0,"Douglas_Fir"=0,"Black_Cottonwood"=0,"Paper_Birch"=0)
  }
  e <- rbind(c,d)
  if(Run_type=="FireReturn"){
    Run_Time <- FR_treatments[ID==PlotList[ii],100-TimeSinceFire] #using 100 as mean fire return interval
    #write out the initial conditions as csv with output file paths pasted in
    write.csv(rbind(e,data.table(Init.Values = c("Output","ShortOutput","Timesteps"), 
                                 Interior_Spruce=c("E:\\SORTIE_runs\\","E:\\SORTIE_runs\\",Run_Time)),fill=TRUE),
              paste0(InitTreesPath,PlotList[ii],"Inits.csv"),quote=TRUE,row.names=FALSE) 
  }else{
    Run_Time <- FR_treatments[ID==PlotList[ii],250-TimeSinceFire]
    #write out the initial conditions as csv with output file paths pasted in
    write.csv(rbind(e,data.table(Init.Values = c("Output","ShortOutput","Timesteps"), 
                                 Interior_Spruce=c("E:\\SORTIE_runs\\","E:\\SORTIE_runs\\",Run_Time)),fill=TRUE),
              paste0(InitTreesPath,PlotList[ii],"Inits250.csv"),quote=TRUE,row.names=FALSE)
  }
}

####################################################################################
################## Create the SORTIE parameter files ##################################
####################################################################################
#File paths
xmls_path <- "./Inputs/SORTIEInits/"
param_path <- "./Inputs/SORTIEInits/"

#create list of starting condition files
if(Run_type == "FireReturn"){
  lstFiles <- data.table(type=c(0,rep(1,75)),name=c("SBS-01.xml",
                                                    c(paste0(c(paste0("FR0",seq(1,9,by=1)),paste0("FR",seq(10,75,by=1))),
                                                             "Inits.csv"))))
}else{
  lstFiles <- data.table(type=c(0,rep(1,75)),name=c("SBS-01.xml",
                                                    c(paste0(c(paste0("FR0",seq(1,9,by=1)),paste0("FR",seq(10,75,by=1))),
                                                             "Inits250.csv"))))
}

VariableNames <- read.csv("./Inputs/VariableNames.csv", header=TRUE, strip.white = TRUE)

xmlList <- c()
paramList1 <- vector("list",5)
maxtype <- 0
for (i in 1:nrow(lstFiles)) {
  fn <- as.character(trimws(lstFiles$name[i]))
  itype <- lstFiles$type[i]
  
  if (itype == 0) {
    xmlList <- c(xmlList,fn)
  }
  else {
    paramList1[[itype]] <- c(paramList1[[itype]],list(fn))
  }
  if (itype > maxtype) {maxtype <- itype}
}
numtype <- c()
for (iii in 1:5) {
  numtype <- c(numtype,length(paramList1[[iii]]))
}

ListOfFiles <- c()
for (ix in 1:length(xmlList)) { #start loop over xml files
  
  #read the given xml
  res <- read_xml(paste0(xmls_path,xmlList[ix]))
  #write the xml to a file again (this will put in the missing line breaks)
  write_xml(res, "temp.xml")
  
  #read the newly printed file, this time as lines of text
  tmp <- readLines("temp.xml", encoding="UTF-8")
  xml1 <- gsub("\\\\", "//",tmp)    #reverse the slash marks
  
  #make a vector that contains the length of each file type
  for (ip in 1:numtype[1]) {
    for (ip2 in 1:max(1,numtype[2])) {
      for (ip3 in 1:max(1,numtype[3])) {
        for (ip4 in 1:max(1,numtype[4])) {
          for (ip5 in 1:max(1,numtype[5])) {
            ip_vals <- c(ip,ip2,ip3,ip4,ip5)
            newname <- ""
            newname <- paste(substr(xmlList[ix],1,nchar(xmlList[ix])-4),"-",substr(paramList1[[1]][ip],1,nchar(paramList1[[1]][ip])-4),sep="")
            #newname <- paste(substr(xmlList[ix],1,nchar(xmlList[ix])-4),sep="")    #remove type 1 files from the naming
            for (iii in 2:5) {
              if (numtype[iii] >0) {
                newname <- paste(newname,"-",substr(paramList1[[iii]][ip_vals[iii]],1,nchar(paramList1[[iii]][ip_vals[iii]])-4),sep="")
              }
            }
            
            #for each of the files, prepare it, and process it
            # note: we have to do all five files each time because we don't know which of the files might have the output directories (which need 'newname')
            
            for (iii in 1:5) {
              if (numtype[iii] > 0) {
                #print(paste("MakeFiles",iii,ip_vals[iii]))
                #print(paramList1[[iii]][ip_vals[iii]])
                xml2 <- ModifyFile(paste0(param_path,paramList1[[iii]][ip_vals[iii]]),xml1)
              } else {
                xml2 <- xml1
              }
              xml1 <- xml2
            } 
            
            xml2 <- gsub("//", "\\\\", xml2)    #turn any forward slashes into back into double backwards slashes
            #write the new file 
            newname <- paste(newname,".xml",sep="")
            writeLines(xml2,paste0(xmls_path,newname))
            ListOfFiles <- c(ListOfFiles, newname)    #store the newly created file in a list so it can be run automatically later
          }
        }
      }
    }
  }
}

##############################################################################
##### read in SORTIE Outputs as csvs of extracted plots from SORTIE extraction
#of full plots and then clip to get results on 20x20
##############################################################################
sortie_out_path <- "D:/SORTIEruns/FireCarbon/extracted/"
Run_type <- "OldGrowth"
#Checking files
dt_table <- data.table()
plotID250 <- paste0(c(paste0("FR0",seq(1,9,by=1)),paste0("FR",seq(10,75,by=1))))
##Checking that all files are there
for(pl in 1:length(plotID250)){
  PlID <- plotID250[pl]
  FireYr <- FR_treatments[ID==plotID250[pl],TimeSinceFire]
  Yrs2ext <- c(seq(0,250,by=1)-FireYr)
  Yrs2ext <- Yrs2ext[Yrs2ext>=0]
  for(i in 1:length(Yrs2ext)){
    if(file.exists(paste0(sortie_out_path,"ext_SBS-01-",plotID250[pl],"Inits250_det_",Yrs2ext[i]))==TRUE){
      next
    }else{
      print(paste0(sortie_out_path,"ext_SBS-01-",plotID250[pl],"Inits250_det_",Yrs2ext[i]," does not exist"))
    }
  }
}

if(Run_type=="OldGrowth"){
  ###### 250 year runs ########
  dt_table <- data.table()
  plotID250 <- paste0(c(paste0("FR0",seq(1,9,by=1)),paste0("FR",seq(10,75,by=1))))
  for(pl in 1:length(plotID250)){
    PlID <- plotID250[pl]
    FireYr <- FR_treatments[ID==plotID250[pl],TimeSinceFire]
    Yrs2ext <- c(seq(0,250,by=10)-FireYr) #extract every 10 years
    Yrs2ext <- Yrs2ext[Yrs2ext>=0]
    for(i in 1:length(Yrs2ext)){
      dt <- fread(paste0(sortie_out_path,"ext_SBS-01-",plotID250[pl],"Inits250_det_",Yrs2ext[i]), 
                  sep="\t", header=T,na.strings = "--", skip=1)
      dt[,':='(timestep = Yrs2ext[i],plotID = PlID)]
      dt <- dt[X >50 & X <150 & Y>50 & Y <150] #clip 1ha in the middle of the stand
      dt_table <- rbind(dt_table,dt,fill=TRUE)
    }
    print(paste("Plot",PlID,"done"))
  }
  dt_table <- dt_table[,SO_sp:=ifelse(Species=="Western_Larch", "Lw", #grow larch as spruce
                                      ifelse(Species=="Douglas_Fir","Fd", #grow Doug fir as spruce
                                             ifelse(Species=="Subalpine_Fir","Bl",
                                                    ifelse(Species =="Interior_Spruce","Sx",
                                                           ifelse(Species=="Lodgepole_Pine","Pl",
                                                                  ifelse(Species=="Trembling_Aspen","At",
                                                                         ifelse(Species=="Black_Cottonwood",
                                                                                "Ac",
                                                                                ifelse(Species=="Paper_Birch",
                                                                                       "Ep",Species))))))))]
  dt_table[,Tree_class:=ifelse(Type=="Seedling"|Type=="Sapling"|Type=="Adult","L","D")] 
  
  ###Tree data from SORTIE
  Trees_SO <- dt_table[!is.na(DBH) & Height >1.3|!is.na(DBH) & is.na(Height)]
  Trees_SO[,Tree_class := ifelse(Tree_class=="L",2,4)] #will have to change when have dead trees initating
  
  g <- vector() #rJust live trees
  for(i in 1:nrow(Trees_SO)){
    g[i] <- TreeCarbonFN(Species = Trees_SO[i,SO_sp], DBH= Trees_SO[i,DBH], 
                         HT= Trees_SO[i,Height], Tree_class=Trees_SO[i,Tree_class])
  }
  # xkg/(100x100m) = xkg/1ha * 1Mg/1000kg = x Mg/1000
  Trees_SO[,CarbonPerHa := g/1000]
  #can't calculate dead carbon here because no heights - need to use the dbh only equations, which aren't in our function yet
  #just write out live trees
  write.csv(Trees_SO[Tree_class==2],"D:/SORTIEruns/SE_Trees_SO_250.csv")

}else{
  ###### 100 year runs ########
  
  dt_table <- data.table()
  for(pl in 1:length(plotID)){
    PlID <- plotID[pl]
    yrs <- seq(0,length(list.files("E:/SORTIE_runs/SortieExtract/",pattern =PlID))-1)
    for(i in 1:length(yrs)){
      dt <- fread(paste0("E:/SORTIE_runs/SortieExtract/Ext_SBS-01-",plotID[pl],"Inits_det_",yrs[i]), 
                  sep="\t", header=T,na.strings = "--", skip=1)
      dt[,':='(timestep = yrs[i],plotID = PlID)]
      dt <- dt[X >50 & X <150 & Y>50 & Y <150] #clip 1ha in the middle of the stand
      dt_table <- rbind(dt_table,dt,fill=TRUE)
    }
    print(paste("Plot",PlID,"done"))
  }
  
  
  dt_table <- dt_table[,SO_sp:=ifelse(Species=="Western_Larch", "Lw", #grow larch as spruce
                                      ifelse(Species=="Douglas_Fir","Fd", #grow Doug fir as spruce
                                             ifelse(Species=="Subalpine_Fir","Bl",
                                                    ifelse(Species =="Interior_Spruce","Sx",
                                                           ifelse(Species=="Lodgepole_Pine","Pl",
                                                                  ifelse(Species=="Trembling_Aspen","At",
                                                                         ifelse(Species=="Black_Cottonwood",
                                                                                "Ac",
                                                                                ifelse(Species=="Paper_Birch",
                                                                                       "Ep",Species))))))))]
  dt_table[,Tree_class:=ifelse(Type=="Seedling"|Type=="Sapling"|Type=="Adult","L","D")] 
  
  ###Tree data from SORTIE
  Trees_SO <- dt_table[!is.na(DBH) & Height >1.3|!is.na(DBH) & is.na(Height)]
  Trees_SO[,Tree_class := ifelse(Tree_class=="L",2,4)] #will have to change when have dead trees initating
  
  g <- vector() #just live trees
  for(i in 1:nrow(Trees_SO)){
    g[i] <- TreeCarbonFN(Species = Trees_SO[i,SO_sp], DBH= Trees_SO[i,DBH], 
                         HT= Trees_SO[i,Height], Tree_class=Trees_SO[i,Tree_class])
  }
  # xkg/(100x100m) = xkg/1ha * 1Mg/1000kg = x Mg/1000
  Trees_SO[,CarbonPerHa := g/1000]
  #just write out every 10 years:
  write.csv(Trees_SO[Tree_class==2],"E:/SORTIE_runs/SortieExtract/SE_Trees_SO_1_75.csv")
}

