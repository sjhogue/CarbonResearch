########################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~To run~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# - Will need shapefile map data saved in CountyMaps called 'countyp010', 'functionsPSUMv2.R',                         #
#   and 'PSUMv2.R' saved in working directory                                                                          #
# - Will also need a .csv file containing the point source data including coordinates, county name,                    #  
#   state name (can be retrieved automatically - see below), emissions (in annual metric tons of CO2)                  #
# - Data must have columns named State, County, Emissions, Lat, Long with these exact names in any order.              #
#   Additional columns do not matter.                                                                                  #
# - Use list.files('.') to see what's in your working directory. Change using 'setwd('Path/to/DesiredDirectory')'      #              #
# - Will need mean locational error for dataset (mean linear distance between reported and actual locations) in km     #
########################################################################################################################

#USER INPUT
err <- 0.842049 #mean locational error for data
res <- 1 #desired resolution in degrees. This cannot be changed in this version.
num <- 10000 #number of points used in the Monte Carlo simulation (more increases runtime but also accuracy)
numPts <- 100 #square root of the number of points distributed over each county when computing average distances
subfolder <- '../../pointsource/Data' #name of subfolder. If not in one, leave as is.
eGRID <- TRUE #If true will automatically convert eGRID data in the format we have it to the format it needs to be in. 

#input proj4string of your .csv data here. MUST be in proj4string format.
#http://www.remotesensing.org/geotiff/proj_list/ look here for information on how to format these
#For example NAD83 is +proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0
projection <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 
data_separator <- ',' #delimeter used by your data. input '\t' for tab separated data, ' ' for white space, ',' for .csv files etc




#filename<-list.files(subfolder)[k] #this is useful for simplifying running the code through a large number of files


#Load Packages
library(spdep)

#Geography packages
library(rgdal) #readOGR
library(stringdist)
library(rgeos) #avgDist (gIntersection)
library(maptools)
library(sp)
library(maps)

#Stats packages
library(mvtnorm) #rmvnorm
require(akima)
library(stats)

#Data structures packages
library(data.table)
library(raster)
require(plyr)
require(grDevices)
 

##### Read in all necessary files #####

#Load file with Monte Carlo function
source('PSUM1deg.R')

#Other useful functions:
source('functionsPSUM1deg.R')


#Read in county shape file data from National Atlas, saved in working directory as 'countyp010.shp'
USCounties <- readOGR(dsn='CountyMaps',layer="countyp010g",stringsAsFactors=FALSE)
names(USCounties@data)[names(USCounties@data)=='NAME'] <- "COUNTY"

#Points are in WGS84 but map data is NAD83, so convert map to same coordinate system
counties84 <- spTransform(USCounties,CRS(projection))

options(warn=2)

for(iter in 38:49)
{
  filename<-list.files(subfolder)[iter]
  #File with columns for latitude, longitude, state abbrev, county name, emissions
  if (filename=='DCeGRID.csv') next
  PS<-read.table(paste(subfolder,'/',filename,sep=''),sep=data_separator, header=TRUE,stringsAsFactors=FALSE)
  
  if (eGRID)
  {
    PS <- converteGRID(PS)
    write.csv(PS,paste('converted',filename,sep=''))
  }
  
  if (is.null(PS$County))
  {
    County <- rep('NA',nrow(PS))
    PS <- cbind(PS,County)
  }
  
  PS <- cbind(seq(1,nrow(PS)),PS,rep(err,nrow(PS)))
  colnames(PS) <- c('ID',colnames(PS)[3:ncol(PS)-1],'Error')
  
  
  states <- unique(PS$State)
  nStates <- length(states) #store for more convenient referal
  
  
  #### Loop through all states in data ####
  
  skipped <- NULL
  wrongCountry <- NULL
  for (i in 1:nStates)
  {
    PS2 <- PS[which(PS$State==states[i]),] #select point sources in state i
    n <- nrow(PS2) #number of point sources in state i
    
    #name of state as found in map data
    stateName <- unique(counties84$STATE[closestMatch(as.character(states[i]),counties84$STATE)])
    
    #Trim the transformed map to only contain county data for the relevant state
    stateMap <- counties84[counties84$STATE==stateName & !is.na(counties84$STATE),]
    
    #to get rid of water borders on maps
    if (stateName != 'DC')
    {
      outerBorder <- which(nchar(stateMap$COUNTY)==2)
      if (length(outerBorder) != 0) stateMap<- stateMap[-outerBorder,]
    }

    
    #nbs <- poly2nb(state) #generate list of neighbors of every polygon in 'stateMap'
    
    
    #### Loop through all point sources listed under each state ####
    tmp <- NULL
    coordsCounty<- NULL
    for (j in 1:n)
    {    
      coordsNA <- FALSE
      #Check if has emissions value
      if (is.na(PS2$Emissions[j])) 
      {
        print('skipping point source with no emissions')
        skipped <- c(skipped, PS2$ID[j])
        next
      }
      
      #check if there is something entered for the county name
      if(nchar(as.character(PS2$County[j]))==0) 
      {
        namedCountyIndex <- NA
        namedCountyName <- NA
      } else 
      {
        namedCountyIndex <- closestMatch(as.character(PS2$County[j]),stateMap$COUNTY)[1]
        namedCountyName <- stateMap$COUNTY[namedCountyIndex]
      }
      
      #check if coordinates are NA
      if (is.na(PS2$Lat[j]) | is.na(PS2$Long[j]))
      {
        coordsNA <- TRUE
        actualCountyName <- NA
        actualCountyIndex <- NA
        if (is.na(namedCountyIndex)){
          print('skipping point source with not enough information')
          skipped <- c(skipped, PS2$ID[j])
          next
        } else errStr <- 'bigger'
      } else #if not, determine which county the coordinates fall into
      {
        #Convert point source data to SpatialPoints object so it can be compared to county data  
        spPS <- SpatialPoints(coords=cbind(PS2$Long[j],PS2$Lat[j]),proj4string=CRS(projection))
        actualCounty <- c(over(spPS,stateMap)$COUNTY,over(spPS,stateMap)$COUNTYP010); #do I still need P010?
        actualCountyName <- actualCounty[1]
        actualCountyIndex <- which(actualCountyName==stateMap$COUNTY)[1]
      }
      
      #if coordinates are given for the point source
      if (!coordsNA){
        if (is.na(actualCountyIndex)) 
        {
          #check if in other state
          locTest <- over(spPS,counties84)
          if (is.na(locTest$STATE))
          {
            print('point source in wrong country')
            wrongCountry <- c(wrongCountry, PS2$ID[j])
          } else if (locTest$STATE != stateName)
          {
            realState <- counties84[counties84$STATE==locTest$STATE & !is.na(counties84$STATE),]
            nameTest <- unique(realState$COUNTY[amatch(as.character(PS2$County[j]),realState$COUNTY,method='lv',maxDist=7)])
            if (is.na(nameTest))
            {
              errStr <- 'bigger'
              countyName <- namedCountyName
            }
            else if (locTest$COUNTY == nameTest) 
            {
              print('skipping point source in wrong state')
              skipped <- c(skipped, PS2$ID[j])
              next
            }
          }
          
          if (is.na(namedCountyIndex)) errStr <- 'biggest'
          else
          {
            errStr <- 'bigger'
            countyName <- namedCountyName
          }
        } else if (is.na(namedCountyIndex)) 
        {
          errStr <- 'normal'
          countyName <- actualCountyName
        } else
        {
          if (namedCountyIndex!=actualCountyIndex)
          {
            errStr <- 'bigger'
          } else errStr <- 'normal'
          countyName <- namedCountyName        
        }   
      }
      
      if (errStr == 'normal')
      {
        map <- stateMap[which(stateMap$COUNTY==countyName),]
      } else if (errStr == 'bigger') 
      {
        map <- stateMap[which(stateMap$COUNTY==countyName),]
        cen <- getCentroid(map)
        PS2$Long[j] <- cen[1]
        PS2$Lat[j] <- cen[2]
        PS2$Error[j] <- avgDist(map,numPts)
      } else if (errStr == 'biggest') 
      {
        map <- nearestExternalCounty(cbind(PS2$Long[j],PS2$Lat[j]),stateMap)
        
        cen <- getCentroid(map)
        PS2$Long[j] <- cen[1]
        PS2$Lat[j] <- cen[2]
        PS2$Error[j] <- avgDist(map,numPts)
      }
      run <- PSUMv2(PS2$Long[j],PS2$Lat[j],PS2$Emissions[j],PS2$Error[j],map,projection,res,num)
      tmp <- rbind(tmp,run)
      coordsCounty <- c(coordsCounty,actualCountyName) #store for output
      print(paste('done with powerplant',j,sep=' '))
    }
    
    summedz <- aggregate(tmp$z,by=list(tmp$x,tmp$y),FUN=sum)
    summedv <- aggregate(tmp$v^2,by=list(tmp$x,tmp$y),FUN=sum)
    setnames(summedv,colnames(summedv),c('x','y','v'))
    summed <- cbind(summedz,sqrt(summedv$v))
    setnames(summed,colnames(summed),c('x','y','z','v'))
    zeroRows <- which(summed$z==0 & summed$v==0)
    if (length(zeroRows)!=0) summed<-summed[-zeroRows,]
    
    print(paste(states[i],'PSUMtest.csv',sep=''))
    write.csv(summed,paste(states[i],'PSUMtest.csv',sep=''))
    
    data <- cbind(PS2$Long,PS2$Lat,PS2$Error)
    ind<-NULL
    if (length(skipped) !=0 ) 
    {
      skippedDF <- NULL
      for (i in 1:length(skipped))
      {
        ind <- c(ind,which(PS2$ID==skipped[i]))
        skippedDF <- rbind(skippedDF,PS2[ind,])
      }
      write.csv(skippedDF,paste(states[i],'SkippedTest.csv',sep=''))
      data <- data[-ind,]
      write.csv(cbind(PS2[-ind,],coordsCounty),'checkingStuff.csv')
    } else write.csv(cbind(PS2,coordsCounty),'checkingStuff.csv')
    
    write.csv(data,paste(states[i],'Pointstest.csv',sep=''))
    
  }  
  
}


