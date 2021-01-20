#######################################################################################################################
#### Thylacine extinction inference using TTSRD - BW Brook, 2020 ####
rm(list=ls()); options(scipen=999,warn=-1); libs <- c("parallel","sp","rgdal","raster", "smooth","maptools","pbapply")
lapply(libs, require, character.only=T)
source('thy_ede_functions.r') # put in the working directory (all data files are assumed to be in working directory)
#######################################################################################################################

#######################################################################################################################
#### Load data, add regionalisation to data frame, prepare spatial files, etc. ####
thyly <- read.csv("ttsrd.csv") # read TTSRD input data and verify
map.grid <- read.csv("tas_grid.csv") # Pre-defined map.grid for Tas at 0.01 degree scale
tas.stack <- stack("tas.stack.12.grd") # Read pre-created stack of SDM raster layers for Tas
IBRA.shape <- shapefile("IBRA_Tas.shp") # Load Tasmanian Interim Bioregionalisation of Australia polygons
IBRA <- list(A="Furneaux",B="Ben Lomond",C="Tasmanian Northern Midlands",D="Tasmanian South East",E="Tasmanian West",
             F="Tasmanian Northern Slopes",G="Tasmanian Southern Ranges",H="Tasmanian Central Highlands",I="King")

hist(thyly$year,breaks=15,main="Thylacine sightings in Tasmania by decade",xlab="",ylab="Number of observations)")
str(thyly); table(thyly$type) # verify/summarise TTSRD
thyly <- add.ibra(thyly) # overlay sightings on IBRA regionalisation, fix beach locations, etc.
map.grid <- map.grid[-c(which(map.grid$Lat>-40.64)),] # delete King and Flinders Island groups
#######################################################################################################################

#######################################################################################################################
#### Set reliability vector (scenario) or select sensitivity analysis ####
es <- set.rel(d=thyly, rel="XA"); head(es) # set sighting probability according to type and return data for analysis
#es <- droplevels(es[which(es$edec=="PS"),]) # Physical only
#es <- droplevels(es[which(es$edec=="PS" | es$edec=="EO" | es$edec=="EI"),]) # Physical + Expert
#es <- droplevels(es[which(es$rating==5 & es$edec!="OI"),]) # PS,EO,EI,OO but rated 5 only
#es <- droplevels(es[which(es$rating>2 & es$edec!="OI"),]) # PS,EO,EI,OO but rated 3-5
table(es$edec,es$rating)
#######################################################################################################################

#######################################################################################################################
### Pre-process and display sighting record data ####
dat <- dmy.cprob(data.frame(year=es$year,prob=es$prob)) # pre-process data
ma.pd(dat,n=20) # plot moving-average of 1-rel as indication of persistence, based on a window of n sightings
#######################################################################################################################

#######################################################################################################################
#### Set EDE parameters ####
ede <- "ole" # Current methods are "LAD", "rw64", so93", "mc06" and "ole" ('f' takes sighting record and returns TE)
rep <- 1e5 # fast for all methods except OLE, suggest 1e3 for fast, 1e5 for slow (more precise)
obs.year <- 2060 # use 2019 for current probability of extinction, or large number (e.g., 3000) to get full UCI
#######################################################################################################################

#######################################################################################################################
#### EDE using all TTSRD records (Tasmania wide) ####
max(es$year); bbj.2019(dat, iter=rep, ey=obs.year, m=ede) # BBJ18 method: m = pre-defined EDE function to use (default: "LAD")
#str(es); table(es$edec,es$rating)
#######################################################################################################################

#######################################################################################################################
#### Spatially discrete analysis based on IBRA with following bioregion options: ####
IBRA.select <- c("A") # see above list for IBRA code, A to I
esd <- es[which(es$IBRA==IBRA[[IBRA.select]]),]
dat.esd <- dmy.cprob(data.frame(year=esd$year,prob=esd$prob))
IBRA[IBRA.select]; length(which(esd$prob>0)); bbj.2019(dat.esd, iter=rep, ey=obs.year, m=ede)

#######################################################################################################################
IS <- which(es$IBRA==IBRA["A"] | es$IBRA==IBRA["B"]) # NE: Ben Lomond, Furneaux = c("A","B")
IS <- which(es$IBRA==IBRA["F"] | es$IBRA==IBRA["I"]) # NW: Northern Slopes, King = c("F","I")
IS <- which(es$IBRA==IBRA["C"] | es$IBRA==IBRA["D"]) # SE: Northern Midlands, South East = c("C","D")
IS <- which(es$IBRA==IBRA["E"] | es$IBRA==IBRA["G"] | es$IBRA==IBRA["H"]) # WHA: West, Southern Ranges, Central Highlands = c("E","G","H")
esd <- es[IS,]; dat.esd <- dmy.cprob(data.frame(year=esd$year,prob=esd$prob))
max(esd$year); length(which(esd$prob>0)); bbj.2019(dat.esd, iter=rep, ey=2020, m="mc06")
#######################################################################################################################

#######################################################################################################################
#### Spatial extinction map (continuous) ####
thy.ss <- es[-which(is.na(es$lat)|es$prob==0),] # set sighting prob., remove zero-prob records and those without spatial data

cl <- makeCluster(detectCores()-1); clusterExport(cl, c("spDistsN1","jr.2014","dmy.cprob","obs.year"))
map.grid$TE <- pbsapply(cl=cl, 1:dim(map.grid)[1], te.grid, ts=thy.ss, grd=map.grid, m="so93", result="UCI"); stopCluster(cl); rm(cl)

rmap <- cbind(map.grid$Lon, map.grid$Lat, pmin(2060,map.grid$TE)); colnames(rmap) <- c("Longitude", "Latitude", "TE")
plot(rasterFromXYZ(rmap),col=rev(terrain.colors(30,0.8)),zlim=c(min(min(map.grid$TE),1920),min(max(map.grid$TE),2060)),
     xlab="Lon E",ylab="Lat S") # plot gridded map, colorized by TE, tune n to have 2 col per decade
write.csv(rmap,"TE-expert-MTE.csv",row.names=F)
#######################################################################################################################

#######################################################################################################################
#### Mapping of aggregated sightings for TTSRD ####
thyly.sub <- thyly[which((thyly$edec=="PS" | thyly$edec=="EO") & thyly$year>=1910 & thyly$year<=2020),]

stn <- 3; tas.stack[[stn]] # layer 3 is vegetation type
thyly.sub$unique <- paste0(thyly.sub$lon,"_",thyly.sub$lat)
un.rec <- unique(thyly.sub$unique)
aggr.thyly <- aggregate(x=thyly.sub$unique, by=list(unique.values=thyly.sub$unique), FUN=length)
coord <- strsplit(aggr.thyly$unique.values, "_")
for(i in 1:length(coord)) {
  aggr.thyly$lon[i] <- coord[[i]][1]
  aggr.thyly$lat[i] <- coord[[i]][2]
}
aggr.thyly <- aggr.thyly[-length(coord),] # remove NAs
image(tas.stack[[stn]],col=rev(terrain.colors(64)),asp=1.25,main="TTSRD (PS/EO) 1910-2020",xlab="Lon (E)",ylab="Lat (S)")
points(aggr.thyly[,3:4],col="red",cex=as.numeric(aggr.thyly$x/2),lwd=2) # zoom the map to make it look clearest
#######################################################################################################################

#######################################################################################################################
#### Extra: Simulation of MTE based on random (within pre-defined range) sighting reliabilities ####
thy <- thyly[which(thyly$IBRA==IBRA[[IBRA.select]]),] # or thy <- thly for full dataset

mte.vec <- rep(NA,rep) # Could improve via a function with mc.apply
for(i in 1:rep) {
  es <- set.ran(thy)
  mte.vec[i] <- bbj.sim(dmy.cprob(data.frame(year=es$year,prob=es$prob),plot=F),m=ede)
  if(i%%round(s/10)==0) print(i)
}
summary(mte.vec); quantile(mte.vec,0.95) # quantile(mte.vec[which(is.na(mte.vec))],0.95)
hist(mte.vec,freq=F,xlab="Year",main="SCEN: Distribution of TE")
#######################################################################################################################
