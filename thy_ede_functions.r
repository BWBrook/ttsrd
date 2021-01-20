#### Functions for data processing and EDE estimation, BBJ method ####

#######################################################################################################################
#### Set reliability vector, PO (base) is default ####
set.rel <- function(d, rel="A") {
  d <- d[,-c(4,7,10:16)]
  d$prob <- NA

  # Ratings: P = 3-5, S = 2-5, O = 1-5, I = 1-5 with: max(thyly[which(thyly$edec=="O"),8])
  # For es, 5 keeps set rating, 4-1 cause fixed or proportional (multiplier) decreases
  # Reference reductions: 5=5/5=1, 4=4/5=0.8, 3=3/5=0.6, 2=2/5=0.4, 1=1/5=0.2
  p.div <- 5
  e.div <- 5
  o.div <- 5
  # incr or decr if want  rating multiplier to be shallower or steeper, e.g. 4 gives 1,1,0.75,0.5,0.25
  pr <- switch(rel,
        # Probabilities assigned to each type:
        #      PS EO    EI   OO    OI
        "P"= c(1, 0.00, 0.0, 0.00, 0.00),
        "EH"=c(1, 0.50, 0.2, 0.00, 0.00),
        "EO"=c(0, 0.25, 0.1, 0.00, 0.00),
        "E"= c(1, 0.25, 0.1, 0.00, 0.00),
        "A"= c(1, 0.25, 0.1, 0.05, 0.01),
        "XA"=c(1, 0.05, 0.01, 0.005, 0.001))
  # P = physical records only, EH = physical/expert high, E = physical/expert default, A = all record types included

  d[d$edec=="PS",8] <- pmin(1,pr[1]*(d[d$edec=="PS",6]/p.div))
  d[d$edec=="EO",8] <- pmin(1,pr[2]*(d[d$edec=="EO",6]/e.div))
  d[d$edec=="EI",8] <- pmin(1,pr[3]*(d[d$edec=="EI",6]/e.div))
  d[d$edec=="OO",8] <- pmin(1,pr[4]*(d[d$edec=="OO",6]/o.div))
  d[d$edec=="OI",8] <- pmin(1,pr[5]*(d[d$edec=="OI",6]/o.div))
  return(d)
}
#######################################################################################################################

#######################################################################################################################
#### Set random reliability vector, based on a priori table ####
set.ran <- function(d) {
  d <- d[,-c(4,7,10:16)]
  d$prob <- NA

  d[d$edec=="PS" & d$rating==5,9] <- runif(dim(d[d$edec=="PS" & d$rating==5,])[1],    1,     1) # PS5
  d[d$edec=="PS" & d$rating==4,9] <- runif(dim(d[d$edec=="PS" & d$rating==4,])[1],  0.9,     1) # PS4
  d[d$edec=="PS" & d$rating==3,9] <- runif(dim(d[d$edec=="PS" & d$rating==3,])[1],  0.8,   0.9) # PS3
  d[d$edec=="EO" & d$rating==5,9] <- runif(dim(d[d$edec=="EO" & d$rating==5,])[1],  0.2,   0.3) # EO5
  d[d$edec=="EO" & d$rating==4,9] <- runif(dim(d[d$edec=="EO" & d$rating==4,])[1],  0.1,   0.2) # EO4
  d[d$edec=="EO" & d$rating==3,9] <- runif(dim(d[d$edec=="EO" & d$rating==3,])[1],    0,   0.1) # EO3
  d[d$edec=="EI" & d$rating==5,9] <- runif(dim(d[d$edec=="EI" & d$rating==5,])[1],    0,   0.1) # EI5
  d[d$edec=="EI" & d$rating==4,9] <- runif(dim(d[d$edec=="EI" & d$rating==4,])[1],    0,  0.08) # EI4
  d[d$edec=="EI" & d$rating==3,9] <- runif(dim(d[d$edec=="EI" & d$rating==3,])[1],    0,  0.05) # EI3
  d[d$edec=="OO" & d$rating==5,9] <- runif(dim(d[d$edec=="OO" & d$rating==5,])[1],    0,   0.1) # OO5
  d[d$edec=="OO" & d$rating==4,9] <- runif(dim(d[d$edec=="OO" & d$rating==4,])[1],    0,  0.08) # OO4
  d[d$edec=="OO" & d$rating==3,9] <- runif(dim(d[d$edec=="OO" & d$rating==3,])[1],    0,  0.06) # OO3
  d[d$edec=="OO" & d$rating==2,9] <- runif(dim(d[d$edec=="OO" & d$rating==2,])[1],    0,  0.04) # OO2
  d[d$edec=="OO" & d$rating==1,9] <- runif(dim(d[d$edec=="OO" & d$rating==1,])[1],    0,  0.02) # OO1
  d[d$edec=="OI" & d$rating==4,9] <- runif(dim(d[d$edec=="OI" & d$rating==4,])[1],    0,  0.03) # OI4
  d[d$edec=="OI" & d$rating==3,9] <- runif(dim(d[d$edec=="OI" & d$rating==3,])[1],    0,  0.02) # OI3
  d[d$edec=="OI" & d$rating==2,9] <- runif(dim(d[d$edec=="OI" & d$rating==2,])[1],    0,  0.01) # OI2
  d[d$edec=="OI" & d$rating==1,9] <- runif(dim(d[d$edec=="OI" & d$rating==1,])[1],    0, 0.005) # OI1

  return(d)
}
#######################################################################################################################

#######################################################################################################################
add.ibra <- function(d) {
  # Fix beach locations that miss IBRA polygons by shifting jittering just inland
  d[125,5:6] <- c(-40.993404, 145.742551) # ID 143
  d[153,5:6] <- c(-43.227800, 145.973654) # ID 171
  d[211,5:6] <- c(-42.059861, 145.261377) # ID 241
  d[253,5:6] <- c(-42.437799, 145.255131) # ID 303
  d[329,5:6] <- c(-43.109792, 145.729889) # ID 390
  d[369,5:6] <- c(-42.855020, 145.439628) # ID 437
  d[427,5:6] <- c(-43.229396, 145.975300) # ID 502
  d[469,5:6] <- c(-42.060872, 145.260024) # ID 550
  d[541,5:6] <- c(-42.350612, 145.481646) # ID 631
  d[570,5:6] <- c(-43.229396, 145.975300) # ID 670
  d[728,5:6] <- c(-42.060872, 145.260024) # ID 832
  d[748,5:6] <- c(-41.228730, 144.688209) # ID 854
  d[761,5:6] <- c(-42.060872, 145.260024) # ID 867
  d[816,5:6] <- c(-42.437799, 145.255131) # ID 929
  d[829,5:6] <- c(-43.229396, 145.975300) # ID 942
  d[876,5:6] <- c(-41.000683, 147.208343) # ID 992
  IBRA <- shapefile("IBRA_Tas.shp") # Load Tasmanian Interim Bioregionalisation of Australia polygons
  pts<-data.frame(x=d$lon,y=d$lat); pts[which(is.na(pts$x)),]<-c(0,0)
  coordinates(pts)<-~x+y; proj4string(pts)<-proj4string(IBRA) # geog points
  plot(IBRA,main="Tasmania IBRA Regions"); points(pts, pch=20, col='red')
  d$IBRA <- unlist(over(pts,IBRA)) # overlay IBRA region on sightings, for spatial analysis
  #d[which(is.na(d$IBRA) & !is.na(d$lat)),] # 16 records getting NA for IBRA due to being on beach - fix by hand edit
  #print(es[which(is.na(es$IBRA)),])
  return(d)
} # 16 records getting NA for IBRA due to on beach - fix by hand, shifting geog
#######################################################################################################################

#######################################################################################################################
#### Combine repeated sightings in a year and return dataframe ####
dmy.cprob <- function(dd,plot=T) {
  dd <- dd[order(dd$year),] # ensure data is ordered by year
  for(i in unique(dd$year)) {
    row.sel <- which(dd$year==i)
    if(length(row.sel) > 1) {
      dd$prob[row.sel[1]] <- 1-prod(1-dd[row.sel,"prob"]) # calculate combined probability for the year
      dd <- dd[-row.sel[-1],] # eliminate all but the first row for each year
    }
  }
  # Tablulate cumulative combinatorial prob. by year, versus last confirmed sighting (similar to model 0)
  dd$cum.prob <- NA
  for(i in 1:length(dd$prob)) {
    dd$cum.prob[i] <- round(1-prod(1-dd$prob[i:length(dd$prob)]),3)
  } # calculate the cumulative probability of persistence, from last 'confirmed' year to last observation
  MTE <- ifelse(max(dd$cum.prob<=0.5),dd$year[min(which(dd$cum.prob<=0.5))],NA) # estimate median TE based on cum.prob

  if(plot==T) {
    plot(dd$cum.prob~dd$year, main="Prob. Persistence based on uncertain records",ylim=c(0,1),
         xlab="Year",ylab="cumulative probability")
    lines(dd$year,dd$cum.prob); abline(v=c(min(dd$year),max(dd$year)),col="blue",lty=2); abline(v=MTE,col="red",lty=1)
    print(dd); print(paste("MTE =",MTE))
  } # Plot the annual vector of combinatorial probabilities

  return(dd)
}
#######################################################################################################################

#######################################################################################################################
#### Plot moving-average of PE based on sighting reliability ####
ma.pd <- function(dd,n=5) {
  cx <- c(0,cumsum(dd$prob)); yr <- dd$year[n:length(dd$year)]
  rsum <- (cx[(n+1):length(cx)]-cx[1:(length(cx)-n)])/n # calculate moving average of probabilities
  plot((1-rsum)~yr, main="Moving average PE based on uncertain records",ylim=c(0,1), xlab="Year",ylab="MA probability of extinction")
  lines(yr,(1-rsum)); abline(v=c(min(yr),max(yr)),col="blue",lty=2)
}
#######################################################################################################################

#######################################################################################################################
#### Jaric & Roberts 2014 method applied to Solow 1993 EDE ####
jr.2014 <- function(dd,ey=2018,a=0.05,m="LAD") {

  ts <- dd$year-dd$year[1]; rs <- dd$prob; n <- length(ts) # initialise vars
  tf <- ts[1]*rs[1] # Note: if prob(rs[1]) = 1 then tf==ts[1]
  for(i in n:2) tf <- tf + ts[i]*rs[i]*prod(1-rs[(i-1):1])
  t0 <- tf+dd$year[1] # most likely first sighting year, in original units

  ts <- ts[-min(which(rs>0))]; rs <- rs[-min(which(rs>0))]; n <- n-1
  tl <- ts[n]*rs[n] # Note: if prob(t.n) = 1 then tl==t.n
  for(i in 1:(n-1)) tl <- tl + ts[i]*rs[i]*prod(1-rs[(i+1):n])

  tr <- tl-tf # most likely sighting period
  r <- sum(rs) # most likely number of observations
  To <- ey-t0 # most likley total observation period

  switch(m,
   "so93"={mod <- "Solow 1993"; MTE <- (r+1)/r*tr
      UCI <- tr/(a^(1/r)); p <- (tr/To)^r},
   "mc06"={mod <- "McInerney et al. 06"; MTE <- ceiling(tr+log(0.5,1-(r/tr)))
      UCI <- ceiling(tr+log(a,1-(r/tr))); p <- (1-(r/tr))^(To-tr)},
   "ss89"={mod <- "Strauss & Sadler 1989"; MTE <- (r*tr)/(r-1)
      UCI <- tr+(a^(-1/(r-1))-1)*tr; p <- ((To-tr)/tr+1)^-(r-1)},
   {mod <- "Most likely endpoint (LAD)"; MTE <- tr; UCI <- NA; p <- NA})

  return(list(Method=mod, MTE=round(t0+MTE),
              UCI=ifelse(is.na(UCI),NA,round(t0+UCI)), p=p))
}
#######################################################################################################################

#######################################################################################################################
#### EDE functions for BBJ18 Method, returning TE (additions here) ####
rw64 <- function(ts) {
  n <- length(ts)
  return(ts[n]+(ts[n]-ts[n-1]))
} # Robson & Whitlock 1964, returning TE

so93 <- function(ts) {
  sr <- ts-ts[1]; n <- length(sr)-1
  return(ts[1]+(n+1)/n*sr[n+1])
} # Solow 1993 (+ Strauss & Sadler 1989), returning TE

mc06 <- function(ts) {
  sr <- ts-ts[1]; n <- length(sr)-1
  return(ts[1]+ceiling(sr[n+1]+log(0.5,1-(n/sr[n+1]))))
} # McInerney et al. 2006, returning TE

ole <- function(ts) {
  gam.fit <- function(i,j,v) (gamma(2*v+i)*gamma(v+j))/(gamma(v+i)*gamma(j))
  sights <- rev(sort(ts)); k <- length(sights)
  v <- (1/(k-1))*sum(log((sights[1]-sights[k])/(sights[1]-sights[2:(k-1)])))
  lambda <- outer(1:k,1:k,gam.fit,v=v); lambda <- ifelse(lower.tri(lambda),lambda,t(lambda))
  e <- matrix(rep(1,k),ncol=1)
  a <- try(as.vector(solve(t(e)%*%solve(lambda)%*%e))*solve(lambda)%*%e,silent=TRUE)
  return(try(round(sum(t(a)%*%sights)),silent=TRUE))
} # Roberts & Solow 2003 OLE, returning TE
#######################################################################################################################

#######################################################################################################################
#### BBJ 2019 computational method applied to LAD, Solow 1993 and OLE ####
bbj.2019 <- function(dd,iter=10000,ey=2018,m="LAD",plot=T,OS="W") {

  ifelse(exists(m), ede <- match.fun(m), ede <- function(ts) max(ts)+1)

  r.rec <- function(d) {
    sr <- d$year[which(d$prob>runif(length(d$year)))] # select random SR using reliability
    ext <- if(length(sr)>2) ede(sr[order(sr)])
    return(ifelse(length(ext)==0,NA,ext))
  } # generates a random sighting record and returns TE

  if(OS=="W") {
    cl <- makeCluster(detectCores()-1); dd # run loop with parallel computing
    ext.yr <- parSapply(cl=cl, 1:iter, function(x) r.rec(dd)) # vector of TE from random subsamples
    stopCluster(cl) # stop parallel computing
  } # parallel processing for Windows OS
    else {
      ext.yr <- mclapply(1:iter,mc.cores=detectCores()-1, function(x) r.rec(dd)) # Linux multi-core option
      ext.yr <- do.call(rbind,ext.yr)
    } # parallel processing for Linux

  ext.yr <- na.omit(as.numeric(ext.yr)) # eliminate bad runs
  ext.yr <- na.omit(ext.yr); min.y <- round(min(ext.yr)) # exclude non-convergence, calc. first ext yr
  if(sd(ext.yr)==0) return(list(MTE=round(ext.yr[1]),UCI=round(ext.yr[1]),
                                PP=ifelse(round(ext.yr[1])>=ey,1,0))) # for only physical records only

  freq.ext <- hist(ext.yr, breaks=round(max(ext.yr))-min.y, plot=F) # use hist to calc. freq per year
  bbj <- list(); bbj$Method <- m; x <- min.y:ey; pp.vec <- rep(NA,length(x)); index <- 1 # initialisation of variables
  for(i in min.y:ey) {
    pp.vec[index] <- 1 - (sum(freq.ext$counts[which(freq.ext$mids<i)])/length(ext.yr))
    index <- index + 1
  } # tally cumulative proportion of TE by year since first TE

  if(round(mean(ext.yr))>=ey) {
    bbj$MTE <- which(pp.vec>0.5)
    print("Median TE reported because mean TE is later than end year")
  } else bbj$MTE <- round(mean(ext.yr)) # Mean of distribution of TE from r.rec sampling
  bbj$LCI <- ifelse(x[max(which(pp.vec>0.975))]==ey,NA,x[max(which(pp.vec>0.975))]) # Lower 2.5th percentile of distribution
  bbj$UCI <- ifelse(x[max(which(pp.vec>0.025))]==ey,NA,x[max(which(pp.vec>0.025))]) # Upper 2.5th percentile of distribution
  bbj$end_year <- pp.vec[length(x)]; names(bbj$end_year) <- paste("PP",ey) # proportion of simulations with TE>ey

  if(plot==T) {
    par(mfrow=c(1,2))
    hist(ext.yr,main="Probability Density of Extinction", xlab="Year",freq=F) 
    abline(v=bbj$MTE,col="blue",lty=2); abline(v=bbj$UCI,col="red",lty=3)
    plot(pp.vec~x, main="Cumulative Probability of Persistence", xlab="Year", ylab="Probability",xlim=c(min(x),ey-1),type='l')
    lines(x,pp.vec); abline(v=bbj$MTE,col="blue",lty=2); abline(v=bbj$UCI,col="red",lty=3)
    par(mfrow=c(1,1))
  }
  return(bbj)
}
#######################################################################################################################

#######################################################################################################################
bbj.sim <- function(dd,m="LAD") {

  ifelse(exists(m), ede <- match.fun(m), ede <- function(ts) max(ts)+1)

  r.rec <- function(d) {
    sr <- d$year[which(d$prob>runif(length(d$year)))] # select random SR using reliability
    ext <- if(length(sr)>2) ede(sr[order(sr)])
    return(ifelse(length(ext)==0,NA,ext))
  } # generates a random sighting record and returns TE

  return(round(r.rec(dd)))
}
#######################################################################################################################

#######################################################################################################################
#### MTE by year method using year-cutoff or jackknifing ####
bbj.mte_yr <- function(dd,iter=10000,ey=2018,m="LAD",jk=F,srm=0,plot=T) {
  dd <- dd[which(dd$prob>=srm),]; print(dd) # exclude sightings below minimum sighting reliability
  n <- dim(dd)[1]; if(n<3) return("Too few observations!")

  if(jk==F) { # assess impact of sequentially assuming all later years are false
    last.P <- max(which(dd$prob==max(dd$prob))) # start with last 'confirmed' record
    x <- dd$year[last.P:n]
    cull.pp <- rep(NA,length(x))
    for(i in last.P:n) {
      cull.pp[i-last.P+1] <- bbj.2018(dd[1:i,],iter,ey,m,plot=F)$MTE
      if(i%%5==0) print(paste("Last Data Year =", dd$year[i], "MTE =", cull.pp[i-last.P+1]))
    }
    if(plot==T) {
      plot(cull.pp~x, main="MTE by Year", xlab="Last Data Year", ylab="MTE at Year")
      lines(x[order(x)],cull.pp[order(x)]); abline(0,1,lty=3)
    }
    return(data.frame(year=x,MTE=cull.pp))
  }

  else { # assess impact of dropping each year in turn (jackknife)
    cull.pp <- rep(NA,n-1)
    print(paste(n,"observations, testing sensitivity of MTE to dropping each in turn"))
    for(i in 1:(n-1)) {
      cull.pp[i] <- bbj.2018(dd[-i,],iter,ey,m,plot=F)$MTE
      if(i%%5==0) print(i)
    }
    if(plot==T) hist(cull.pp, main="Jackknife freq. of predicted TE", xlab="Year")
    return(list(MTE=round(median(cull.pp)),UCI=round(quantile(cull.pp,0.95))))
  }
}
#######################################################################################################################

#######################################################################################################################
#### Make a frequency time series ####
make.ts <- function(d=thyly,timespan=1910:2020) {
  x <- data.frame(year=timespan,count=0)
  y <- aggregate(x=d,by=list(d$year),FUN=length)[1:2]; colnames(y)<-c("year","count")
  for(i in unique(y$year)) x[x$year==i,2] <- y[y$year==i,2]
  return(x)
} 
#######################################################################################################################

#######################################################################################################################
#### Return the TE of a grid call in a spatial map ####
te.grid <- function(x, ts, grd, a=2.15, b=1.074, p=0.01, m="so93", result="MTE"){
  r <- ifelse(result=="UCI",3,2) # report MTE (default) or UCI
  ts$v <- spDistsN1(as.matrix(ts[,c(5,4)]),as.matrix(grd[x,1:2]), longlat=TRUE) # dist of records from grid cell
  ts$v <- pmin(1,a*b^-ts$v) # 1 for all distances up to 10 km, declining expon to 0.01 at 75 km (for default parameters)
  ts <- ts[which(ts$v>p),] # truncate all records further than d km (~75 km for default parameters)
  te <- jr.2014(dmy.cprob(data.frame(year=ts$year,prob=ts$prob*ts$v),plot=F),ey=obs.year,m="so93")[[r]] # LAD, so93 or mc06
  te <- ifelse(length(te)>0,te,NA)  # to avoid crashing if NA value returned by jr.2014
  return(te)
} 
#######################################################################################################################
