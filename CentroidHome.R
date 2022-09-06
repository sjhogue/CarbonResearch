library(ggplot2)
library(rgdal)    # for readOGR(...)
library(rgeos)    # for gIntersection(...)
library(stats)

county<- "Orange County"
state<-"NC"

setwd("~/Dropbox/Carbon/pointsource")
map   <- readOGR(dsn=".",layer="countyp010")
NC    <- map[map$COUNTY==county & map$STATE==state & !is.na(map$COUNTY),]
NC.df <- fortify(NC)

bbox  <- bbox(NC)
x   <- seq(bbox[1,1],bbox[1,2],length=100)   # longitude
y   <- seq(bbox[2,1],bbox[2,2],length=100)   # latitude
all <- SpatialPoints(expand.grid(x,y),proj4string=CRS(proj4string(NC)))
pts <- gIntersection(NC,all)                # points inside the polygons
pts <- data.frame(pts@coords)               # ggplot wants a data.frame
centroid <- data.frame(x=mean(pts$x),y=mean(pts$y))

#guess we don't really need to see it now that we know it works
#ggplot(NC.df)+
#  geom_path(aes(x=long,y=lat, group=group),    colour="grey50")+      
#  geom_polygon(aes(x=long,y=lat, group=group), fill="lightgreen")+
#  geom_point(data=pts, aes(x,y), colour="blue")+
#  geom_point(data=centroid, aes(x,y), colour="red", size=5)+
#  coord_fixed()

r<-nrow(pts)
nrow(pts)
v<-matrix(centroid,nrow=r,ncol=2,byrow=TRUE)
p<-matrix(pts,nrow=r,ncol=2,byrow=TRUE)
dist(as.matrix(p))#v<-data.frame(v)
#nrow(v)
#dimnames(v) <- list(NULL, NULL)
#dimnames(p) <- list(NULL, NULL)

##########################################

euclidean_distance <- function(p,v){
  sqrt(sum((p - v)^2))
}
euclidean_distance(p,v)

##########################################

sqrt((p-v)^2)
p[,1]-v[,1]
p*v
