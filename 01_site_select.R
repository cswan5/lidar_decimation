## Select sites from each lidar point cloud with point densities > 150 pts/m2

library(lidR)
library(raster)
library(tidyverse)

## List files
las.list <- list.files("./Lidar/Lidar_subsets/",pattern="\\.las",full.names = TRUE)[-3]

for(file in 1:length(las.list)){
  ## Get DTM name
  dtm.name <- gsub("\\.las","_dtm.tif",
                   gsub("Lidar_subsets","DTM",las.list[file]))
  
  ## Set out file name
  out.norm <- gsub("\\.las","_norm.las",
                   gsub("Lidar_subsets","Normalized",las.list[file]))
  
  ## Read in las file
  las <- readLAS(las.list[file])
  
  ## Read in DTM
  dtm <- raster(dtm.name)

  ## Move all points in Assateague lidar dataset to 0
  if(str_detect(las.list[file],"assa")){
    las$Z <- las$Z - min(las$Z)
  }
  
  ## Filter out any duplicate points
  las <- filter_duplicates(las)

  ## Normalize point cloud
  las.norm <- normalize_height(las,dtm)
  
  ## Filter out gross statistical outliers
  las.norm<-filter_poi(las.norm,Z>(0)&Z<mean(Z)+6*sd(Z))
  
  ## Write out normalized point cloud
  writeLAS(las.norm,out.norm)
} 


## Grid sizes
sizes <- c(25,20,15,10)

norm.lst <- list.files("./Lidar/Normalized",pattern="\\.las",full.names = TRUE)

## Empty list to hold results
las.dens <- vector(mode="list",length=length(norm.lst))

## Make a grid of point density for each pixel, turn rasters into data frames,
## filter out points that have < 150 pts/m2
for(i in 1:3){
  
  las <- readLAS(norm.lst[i])
  for(j in 1:length(sizes)){
    las.dens[[i]][[j]]<-as.data.frame(grid_density(las,res=sizes[j]),xy=TRUE) %>% 
      filter(density >= 150)
  }
  
}

## Create a list of centroids for the 25 x 25 m plots at each site
centroids <- list(las.dens[[1]][[1]], las.dens[[2]][[1]],las.dens[[3]][[1]])

## Check that plots have data in them
for(i in 1:length(norm.lst)){
  las <- readLAS(norm.lst[i])
  
  buff = 25
  ## Make buffer polygons around the filtered points
  yPlus <- centroids[[i]]$y+(buff/2)
  xPlus <- centroids[[i]]$x+(buff/2)
  yMinus <- centroids[[i]]$y-(buff/2)
  xMinus <- centroids[[i]]$x-(buff/2)
  
  # calculate polygon coordinates for each plot centroid. 
  square=cbind(xMinus,yPlus,  # NW corner
               xPlus, yPlus,  # NE corner
               xPlus,yMinus,  # SE corner
               xMinus,yMinus, # SW corner
               xMinus,yPlus)  # NW corner again - close ploygon
  
  ID = 1:nrow(centroids[[i]])
  
  buff.polys <- SpatialPolygons(mapply(function(poly, id) {
    xy <- matrix(poly, ncol=2, byrow=TRUE)
    Polygons(list(Polygon(xy)), ID=id)
  }, 
  split(square, row(square)), ID),
  proj4string=CRS(st_crs(las)$proj4string))
  
  chm <- rasterize_canopy(las,p2r(),res=1)
  plot(chm)
  plot(buff.polys,add=TRUE)
}

## Keep only plots that have data
assa1 <- c(5:7,12,13,17:22,24:29,32:35,38,39)
assa2 <- c(5:8)
krc09 <- c(1,3:14,16:29,31:46,48:63,65:80,82:97,99:114,116:129,
           131:139,141:143,145:166)

centroids <- list(centroids[[1]][assa1,],
                  centroids[[2]][assa2,],
                  centroids[[3]][krc09,])

## Clip las to different plot sizes
for(i in 3:length(norm.lst))  {
  ctg <- readLAScatalog(norm.lst[i])
  
  file.name.lng <- str_split(norm.lst[i],pattern="/")[[1]][4]
  file.name <- gsub("_norm\\.las",'',file.name.lng)
  
  for(size in 1:length(sizes)){
    buff = sizes[size]
    ## Make buffer polygons around the filtered points
    yPlus <- centroids[[i]]$y+(buff/2)
    xPlus <- centroids[[i]]$x+(buff/2)
    yMinus <- centroids[[i]]$y-(buff/2)
    xMinus <- centroids[[i]]$x-(buff/2)
    
    # calculate polygon coordinates for each plot centroid. 
    square=cbind(xMinus,yPlus,  # NW corner
                 xPlus, yPlus,  # NE corner
                 xPlus,yMinus,  # SE corner
                 xMinus,yMinus, # SW corner
                 xMinus,yPlus)  # NW corner again - close ploygon
    
    ID = rownames(centroids[[i]])
    
    buff.polys <- SpatialPolygons(mapply(function(poly, id) {
      xy <- matrix(poly, ncol=2, byrow=TRUE)
      Polygons(list(Polygon(xy)), ID=id)
    }, 
    split(square, row(square)), ID),
    proj4string=CRS(st_crs(las)$proj4string))
    
    opt_output_files(ctg) <- paste0("./Lidar/Buffers_new/",file.name,"_{XCENTER}_{YCENTER}_buffer_",buff)
    
    ##Clip lidar plots
    las.clp <- clip_roi(ctg,buff.polys)
  }
} 


## Go through lidar plots to ensure they are all high enough density

plot.lst <-  list.files("./Lidar/Buffers_new/",full.names=TRUE)

plot.stats <- data.frame(file = character(0),
                         density = numeric(0))

for(i in 1:length(plot.lst)){
  las <- readLAS(plot.lst[i])
  plot.stats[i,1] <- plot.lst[i]
  plot.stats[i,2] <- density(las)
}

## Filter out plots that have density < 150
las.low.dens <- plot.stats %>% filter(density < 150)
las.high.dens <- plot.stats %>% filter(density >= 150)

## Get just the file names for the low density plots
las.low.dens.lst <- las.low.dens$file

## Get the plot names for plots that should be removed
files.rm <- paste(unique(gsub("buffer_\\d+\\.las","",las.low.dens.lst)),collapse="|")

## Use stringr to remove plots that have low densities
plot.lst.sm <- plot.lst[!str_detect(plot.lst,pattern=files.rm)]

write.csv(plot.lst.sm,"./Lidar/plot_list.csv")
