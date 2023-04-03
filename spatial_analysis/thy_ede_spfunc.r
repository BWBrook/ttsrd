#### Functions for data processing spatial BBJ method ####

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
        "XA"=c(1, 0.05, 0.01, 0.005, 0.001),
        "S5"=c(1, 0.1, 0, 0, 0)) 
  # P = physical records only, EH = physical/expert high, E = physical/expert default, A = all record types included

  d[d$edec=="PS","prob"] <- pmin(1,pr[1]*(d[d$edec=="PS","rating"]/p.div))
  d[d$edec=="EO","prob"] <- pmin(1,pr[2]*(d[d$edec=="EO","rating"]/e.div))
  d[d$edec=="EI","prob"] <- pmin(1,pr[3]*(d[d$edec=="EI","rating"]/e.div))
  d[d$edec=="OO","prob"] <- pmin(1,pr[4]*(d[d$edec=="OO","rating"]/o.div))
  d[d$edec=="OI","prob"] <- pmin(1,pr[5]*(d[d$edec=="OI","rating"]/o.div))
  return(d)
}

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

#### Define distance-decay function ####
tr.exp.decay <- function(a=2, b=1, d=10) min(1, a*b^-d)

#### Set objective function to minimize deviation from required distance-frequency relationship ####
objective <- function(params) {
  
  ted.opt <- function(a, b, x, f) {
    pr.freq <- sapply(x, tr.exp.decay, a=a, b=b)
    return(sum((f - pr.freq)^2))
  } 
  
  a <- params[1]
  b <- params[2]
  return(ted.opt(a, b, d, rf))
}

#### EDE functions for BBJ18 Method, returning TE (additions here) ####
so93 <- function(ts) {
  sr <- ts-ts[1]; n <- length(sr)-1
  return(ts[1]+(n+1)/n*sr[n+1])
} # Solow 1993 (+ Strauss & Sadler 1989), returning TE

mc06 <- function(ts) {
  sr <- ts-ts[1]; n <- length(sr)-1
  return(ts[1]+ceiling(sr[n+1]+log(0.5,1-(n/sr[n+1]))))
} # McInerney et al. 2006, returning TE

#### Jaric & Roberts 2014 method applied to McInerney 2006, Solow 1993 EDE, or LAD ####
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
   {mod <- "Most likely endpoint (LAD)"; MTE <- tr; UCI <- NA; p <- NA})

  return(list(Method=mod, MTE=round(t0+MTE),
              UCI=ifelse(is.na(UCI),NA,round(t0+UCI)), p=p))
}

#### Return the TE of a grid call in a spatial map ####
te.grid <- function(x, ts, grd, a=2.15, b=1.074, p=0.01, ede="mc06", result="MTE"){
  r <- ifelse(result=="UCI",3,2) # report MTE (default) or UCI
  ts$v <- spDistsN1(as.matrix(ts[,c(5,4)]),as.matrix(grd[x,1:2]), longlat=TRUE) # dist of records from grid cell
  ts$v <- pmin(1,a*b^-ts$v) # 1 for all distances up to 10 km, declining expon to 0.01 at 75 km (for default parameters)
  ts <- ts[which(ts$v>p),] # truncate all records further than d km (~75 km for default parameters)
  te <- jr.2014(dmy.cprob(data.frame(year=ts$year,prob=ts$prob*ts$v),plot=F),ey=obs.year,m=ede)[[r]] # LAD, so93 or mc06
  te <- ifelse(length(te)>0,te,NA)  # to avoid crashing if NA value returned by jr.2014
  return(te)
} 
