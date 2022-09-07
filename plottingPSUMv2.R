stateName <- 'NC'
projection <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 

#Read in county shape file data from National Atlas, saved in working directory as 'countyp010.shp'
USCounties <- readOGR(dsn='CountyMaps',layer="countyp010",stringsAsFactors=FALSE)

#Points are in WGS84 but map data is NAD83, so convert map to same coordinate system
counties84 <- spTransform(USCounties,CRS(projection))


#Read in old points
oPts <- read.csv(paste('converted',stateName,'eGRID.csv',sep=''))
oPts <- oPts[oPts$State==stateName,]
oPtsDF <- SpatialPointsDataFrame(coords=cbind(oPts$Long,oPts$Lat),data=oPts,proj4string=CRS(projection))

nPts <- read.csv(paste(stateName,'Pointstest.csv',sep=''))[,2:3]
nPtsDF <- SpatialPointsDataFrame(coords=nPts,data=data.frame(nPts),proj4string=CRS(projection))


data <- read.csv(paste(stateName,'PSUMtest.csv',sep=''))
data$x <- data$x
data$y <- data$y
grid <- SpatialPointsDataFrame(coords=cbind(data$x,data$y),data=data.frame(data$z),proj4string=CRS(projection))
gridded(grid) <- TRUE
r <- raster(grid)

plot(r)
plot(counties84[which(counties84$STATE==stateName),],add=TRUE)
points(oPtsDF,col='black',pch=20)
points(nPtsDF,col='red',pch=20)

#to see the original points if they fell outside of the previous plot:
plot(oPtsDF)
plot(counties84[which(counties84$STATE==stateName),],add=TRUE)