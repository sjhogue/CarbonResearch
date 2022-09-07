#Handy functions (used in runPSUMv2.R)

#Convert eGRID data to input format
converteGRID <- function(data){
  ED <-data[,c(1,2,3,7)]
  colnames(ED) <- c('State','Lat','Long','Emissions')
  
  CD <- read.csv('eGRIDUSCounties.csv',stringsAsFactors=FALSE);
  colnames(CD)<-c('x','y','cnty','state');
  
  coordsNA <- (nchar(ED$Long)==2 | nchar(ED$Lat)==2)
  
  County <- NULL;
  removed <- NULL;
  for (j in 1:nrow(ED))
  {
    if (!coordsNA[j])
    {
      #latMatch gives indices in the county data file of the latitudes that are the same as that of the first point source in ED
      latMatch <- which(ED$Lat[j]==CD$x);
      #gives the index of the entry in CD which is the same as point source j in ED
      match <- latMatch[which(ED$Long[j]==CD$y[latMatch])];
      if (length(match)==0) 
      {
        County<-c(County,'NA')
      }else County <- c(County,CD$cnty[match[1]]);
    } else
    {
      print('NA geographic coordinates')
      removed <- c(removed,j)
    }
  }
  if (length(removed)!=0)
  {
    write.csv(ED[removed,],paste(ED$State[1],'NACoordinates.csv',sep=''))
    ED <- ED[-removed,]
  }
    
  ED <- cbind(ED,County)  
  return(ED)  
}

#string matching - not using, may take out
closestMatch = function(string, stringVector){
  maxD <- max(nchar(stringVector))
  dist <- stringdist(string,stringVector,method='lv',maxDist = maxD)  
  return(which(dist==min(na.omit(dist))))  
}

#Not using currently
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

#convert from degrees to radians. Coordinate can be a number, list of numbers, or table/dataframe
toRad <- function(coordinate)
{
  if (class(coordinate)=="numeric") 
  {
    return((coordinate*pi)/180)
  } else
  {
    radTable <- NULL
    for (j in 1:nrow(coordinate))
    {
      radians <- rep(0,length(coordinate[j,]))
      for (i in 1:length(coordinate[j,]))
      {
        radians[i] <- (coordinate[j,i]*pi)/180
      }
      radTable <- rbind(radTable,radians)   
    }
    rownames(radTable) <- seq(1,nrow(radTable))
    return(radTable);
  }
}

#Compute average distance from the centroid of map to the edge of map
avgDist <- function(map,numPts)
{ 
  polys <- do.call(rbind,lapply(map@polygons[[1]]@Polygons,function(x)c(x@labpt,x@area)))
  polys <- data.frame(polys)
  colnames(polys) <- c("long","lat","area")
  polys$area <- with(polys,area/sum(area))
  cen <- with(polys,c(x=sum(long*area),y=sum(lat*area))) 
  #create boundary box
  bbox  <- bbox(map)
  
  #adding points.
  x   <- seq(bbox[1,1],bbox[1,2],length=numPts)   # longitude
  y   <- seq(bbox[2,1],bbox[2,2],length=numPts)   # latitude
  all <- SpatialPoints(expand.grid(x,y),proj4string=CRS(proj4string(map)))
  pts <- gIntersection(map,all)                # points inside the polygons
  coords <- pts@coords
  
  #uncomment below if you wish to plot
  #plot(map)
  #points(all,col='black',pch='.')
  #points(coords,col='red',pch='.')
  
  #convert to radians
  cenRad <- toRad(cen)  
  ptsRad <- toRad(coords)

  q <- (sin(abs(ptsRad[,2]-cenRad[2])/2))^2+cos(cenRad[2])*cos(ptsRad[,2])*(sin(abs(ptsRad[,1]-cenRad[1])/2))^2
  w <- 2*6371*asin(sqrt(q))
  return(mean(w))
}

getCentroid <- function(map)
{
  polys <- do.call(rbind,lapply(map@polygons[[1]]@Polygons,function(x)c(x@labpt,x@area)))
  polys <- data.frame(polys)
  colnames(polys) <- c("long","lat","area")
  polys$area <- with(polys,area/sum(area))
  centroid <- with(polys,c(x=sum(long*area),y=sum(lat*area)))
  return(centroid)
}

#get coordinates of polygon data x
allcoordinates = function(x) {
  ret = NULL
  polys = x@polygons
  for(i in 1:length(polys)) {
    pp = polys[[i]]@Polygons
    for (j in 1:length(pp)){
      ret = rbind(ret, coordinates(pp[[j]]))
    }
  }
  return(ret)
}

#finds the nearest county if a point is not in the state or is in the water
nearestExternalCounty <- function(point,map)
{
  minD <- NULL
  for (i in 1:length(map)){
    q <- allcoordinates(map[i,])
    d <- spDistsN1(as.matrix(q), point, longlat = TRUE)
    minD <- c(minD,min(d))
  }
  minInd <- which(minD==min(minD))
  return(map[minInd,])
}