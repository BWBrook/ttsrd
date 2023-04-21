#######################################################################################################################
#### Thylacine extinction inference using TTSRD - BW Brook, 2020 ####
rm(list=ls()); options(scipen=999,warn=-1) 
libs <- c("parallel","sp","pbapply","gstat","MASS","raster","ggplot2"); lapply(libs, require, character.only=T)
source('thy_ede_spfunc.r') # put in the working directory (all data files are assumed to be in working directory)
#######################################################################################################################

#### Load data, add regionalisation to data frame, prepare spatial files, etc. ####
thyly <- read.csv("ttsrd.csv") # read TTSRD input data and verify
map.grid <- read.csv("tas_grid.csv") # Pre-defined map.grid for Tas at 0.01 degree scale
map.grid <- map.grid[-c(which(map.grid$Lat>-40.64)),] # delete King and Flinders Island groups
es <- set.rel(d=thyly, rel="A"); head(es) # set sighting probability according to type and return data for analysis

#######################################################################################################################
#### Spatial extinction map (continuous) ####

# Optimize spatial decay parameters based on match to required distance-frequency relationship
d <- c(10,20,30,50,75) # distance markers in km from sighting
rf <- c(1,0.5,0.25,0.05,0.01) # relative weighting of sighting by distance
# alternative rf vectors: c(1,0.5,0.25,0.05,0.01) [best] c(0.5,0.25,0.12,0.06,0.01) [rapid] or c(1,0.75,0.5,0.25,0.1) [slow]
opt_result <- optim(c(5, 2), objective, method = "BFGS"); opt_result
round(sapply(d, tr.exp.decay, a=opt_result$par[1], b=opt_result$par[2]),3) # show weighting by distance

#######################################################################################################################
ede <- "mc06" # Current methods are "LAD", "rw64", so93", "mc06" and "ole" ('f' takes sighting record and returns TE)
obs.year <- 1940 # maximum value to be plotted on the spatial maps
ede.result <- "MTE" # median time to extinction (other option is "UCI" for upper confidence internval)

start.year <- 1910
end.year <- 1929

thy.ss <- es[es$year>=start.year & es$year<=end.year,] # select span of dates to include
thy.ss <- thy.ss[-which(is.na(thy.ss$lat)|thy.ss$prob==0),] # set sighting prob., remove zero-prob records and those without spatial data

cl <- makeCluster(detectCores()-1); clusterExport(cl, c("spDistsN1","jr.2014","dmy.cprob","obs.year"))
map.grid$TE <- pbsapply(cl=cl, 1:dim(map.grid)[1], te.grid, ts=thy.ss, grd=map.grid, 
                        a=opt_result$par[1], b=opt_result$par[2], p=0.01, ede=ede, result=ede.result)
stopCluster(cl); rm(cl) # stop and remove parallel processing cluster

map.grid$TE[is.na(map.grid$TE)] <- start.year # a floor on TE values

rmap <- cbind(map.grid$Lon, map.grid$Lat, pmin(obs.year,map.grid$TE)); colnames(rmap) <- c("Longitude","Latitude","TE")
plot(rasterFromXYZ(rmap),col=rev(terrain.colors(30,0.8)),
     zlim=c(start.year, obs.year), xlab="Lon E",ylab="Lat S") # plot gridded map, colorized by TE

#######################################################################################################################
# Kernel desnity estimate
start.year <- 1910
end.year <- 1929
min.rec.prob <- 0.2 # mininum probability for a record to be included (using expert only)
# 1910–1929, 1930–1949, 1950–1969, and 1970–1989. 
thy.ss <- es[es$year>=start.year & es$year<=end.year,] # select span of dates to include
thy.ss <- thy.ss[-which(is.na(thy.ss$lat)|thy.ss$prob==0),] # set sighting prob., remove zero-prob records and those without spatial data
thy.ss <- thy.ss[which(thy.ss$prob>min.rec.prob),]

dens <- kde2d(thy.ss$lon, thy.ss$lat) # Estimate the kernel density using kde2d
dens_df <- data.frame(expand.grid(lon = kde_result$x, lat = kde_result$y), z = as.vector(kde_result$z))

# Create a SpatialPointsDataFrame objects from the data
coords <- data.frame(lon = dens_df$lon, lat = dens_df$lat)
spdf <- SpatialPointsDataFrame(coords, data.frame(value = dens_df$z))
pred.coords <- data.frame(lon=map.grid$Lon, lat=map.grid$Lat)
pred.spdf <- SpatialPointsDataFrame(pred.coords, data.frame(value=map.grid$TE))

# Inverse-Distance Weighting interpolation to interpolate values to the prediction points
interp <- krige(value ~ 1, spdf, pred.spdf) # interpolate to map grid using kridging
ipred <- interp$var1.pred
iprs <- (ipred - min(ipred)) / (max(ipred) - min(ipred)) # rescale to 0-1 range

# Plot results
rmap <- cbind(interp@coords[,1], interp@coords[,2], 0.7); colnames(rmap) <- c("Longitude","Latitude","KDE")
plot(rasterFromXYZ(rmap),col=rev(terrain.colors(30,0.8)), zlim=c(min(iprs), max(iprs)),
     xlab="Lon E",ylab="Lat S", main="1910-1929", legend=FALSE) # plot gridded map, colorized by TE
contour(dens, nlevels = 10, add = TRUE, col = "red", lwd = 1) # apply the KDE contours

# Raster plot
ggplot() +
  geom_raster(data = map.grid, aes(x = Lon, y = Lat, fill = TE), interpolate = TRUE) +
  scale_fill_gradient(low = "white", high = "red", name = "TE") +
  
  # KDE contours
  geom_contour(data = dens_df, aes(x = lon, y = lat, z = z), color = "blue") +
  
  # Customizing the plot
  theme_minimal() +
  coord_equal() +
  labs(title = "KDE Contours and TE Values", x = "Longitude", y = "Latitude")

#######################################################################################################################
