PSUMv2 <- function(longitude, latitude, emissions, sperr, countyMap, projection, resolution, num){
  lat <- latitude
  long <- longitude
  emit <- emissions
  err <- sperr
  county <- countyMap
  res <- resolution
  decplace <- decimalplaces(res)
  num <- num
  
  #define constants
  e2 <- 0.00669437999014;
  a <- 6378137; #radius of the Earth
  
  # number of km per degree latitude as a function of latitude 
  dlat <- pi*a*(1-e2)/(180*(1-e2*(sin(pi/180*lat))^2)^(3/2))/1000;
  # number of km per degree longitude as a function of latitude
  dlong <- pi*a*(cos(pi/180*lat))/(180*(1-e2*(sin(pi/180*lat))^2)^(1/2))/1000;
  
  # assume that each component contributes equally in Euclidean norm
  corx <- (err/abs(dlong))^2;
  cory <- (err/abs(dlat))^2;
  
  
  mu=c(long,lat); # mean for space - use unrounded numbers
  corrs=matrix(c(corx*0.6366,0,0,cory*0.6366),2,2); 
  
  samplePoints <- rmvnorm(n=num,mean=mu,corrs);
  
  print('finished conversions')
  xt <- samplePoints[,1] #x-coordinates of set of sample points
  yt <- samplePoints[,2] #y-coordinates of set of sample points
  sp <- SpatialPoints(coords=samplePoints,proj4string=CRS(projection))
  
  
  #Set up the data table
  xtRnd <- round(xt,decplace)
  ytRnd <- round(yt,decplace)
  
  xs<-seq(min(xtRnd),max(xtRnd),by=10^(-decplace))
  ys<-seq(min(ytRnd),max(ytRnd),by=10^(-decplace))
 
  grid <- expand.grid(xs,ys) 
  z <- v <- rep(0,nrow(grid));
  
  DF <- data.table(cbind(grid,z,v))
  setnames(DF,colnames(DF),c('x','y','z','v'))
  setkey(DF,x,y)
  
  print('set up data table')  
  
  outBounds <- is.na(over(sp,county)$COUNTY) #whether or not each sample point is outside the bounds of the county

  for (i in 1:num) {  
    if (!outBounds[i]){
      ind <- which(round(DF$x,decplace)==xtRnd[i] & round(DF$y,decplace)==ytRnd[i])  
      DF$z[ind] <- DF$z[ind]+emit;
    }
  }

  print('finished computing z')
  
  # compute average
  counter <- num - sum(outBounds)
  DF$z <- DF$z/counter;
  if (nrow(DF)==1){
    DF$v <- 0
  } else
  {
    DF$v <- emit*sqrt(DF$z/emit)
    ind <- which(round(DF$x,decplace)==round(long,decplace) & round(DF$y,decplace)==round(lat,decplace))
    DF$v[ind] <- emit*sqrt(round(1,8)-round(DF$z[ind]/emit,8))
  }

  
  setnames(DF,colnames(DF),c('x','y','z','v'))
  setkey(DF,x,y)
  return(DF)
}