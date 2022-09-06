system('python readGrid.py gridcar.2009')

a <- read.asciigrid('cdiacOld.txt') #read in cdiac data

ras <- raster(a) #convert to raster

stateName <- 'IA'
data <- read.csv(paste(stateName,'PSUMtest.csv',sep=''))
data$x <- data$x - 0.5
data$y <- data$y - 0.5
grid <- SpatialPointsDataFrame(coords=cbind(data$x,data$y),data=data.frame(data$z),proj4string=CRS(projection))
gridded(grid) <- TRUE
r <- raster(grid)
r <- r*(3/11) #convert to units of C from CO2

c2 <- crop(ras,r) #reduce size of cdiac data to the size of r

sum2 = c2 + r #combine 
