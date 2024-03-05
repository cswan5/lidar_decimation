library(doParallel)
library(tidyverse)

source("./Code/00_lidar_decimation_functions.R")

## Read in the list of lidar plots
plot.lst <- read.csv("./Lidar/plot_list.csv")$x

## Vector of point densities
dens <- c(175,150,125,100,75,50,40,25,20,15,10,5,4,3,2,1)

## Make a data frame with the names of the lidar tiles, point densities, and name
## of the save file
df <- data.frame(las=rep(plot.lst,each=length(dens)),
                 dens = dens) %>% 
  mutate(out.name = paste0(gsub("\\.las","",
                                gsub("\\./Lidar/Buffers_new/","\\./Lidar/Lidar_metrics/Full/",las)),
                           "_density_",dens,".csv"))

## Set number of cores for parallel processing
nc <- detectCores()-1
registerDoParallel(cores=nc)


foreach(j = 1:nrow(df))%dopar%{
  ##Bring in libraries for parallel processing
  library(lidR)
  library(tidyverse)
  
  ##Refer to decimation functions
  source("./Code/00_lidar_decimation_functions.R")
  
  ## Create a data frame to store the metrics
  las.mets <- data.frame()
  
  ## Initialize i for each iteration, j
  i = 1
  
  ## Get the numbers from the file name
  nums <- str_extract_all(df$las[j],"\\d+")
  
  ## Set the x and y coordinates
  x = paste0(nums[[1]][3],".",nums[[1]][4])
  y = nums[[1]][5]
  
  ## Get the buffer from the file name
  buffer = gsub("\\.las","",
                gsub(".*(?<=\\D)(\\d)","\\1",df$las[j],perl=TRUE))

  while(i < 251){
    ## run the lidar decimation algorithm 250 times
    skip <- FALSE
    
    print(i)
    
    ## Run the lidar decimation algorithm - this skips any iterations that 
    ## have errors because of the random decimation protocol
    tryCatch({las.met <- las.sensitivity(lidar=df[j,1], pts=df[j,2],reps=1)
             las.mets <- bind_rows(las.mets,las.met)}, error = function(e)
             { skip <<- TRUE})
    
    ## If there's an error, redo this iteration
    if(skip) {next}
    ## If no error, move to the next iteration
    else {i = i+1}
  }
  
  ## Add x and y coordinates, buffer, and point density to the output
  las.mets <- las.mets %>%
    mutate(x = x, 
           y = y,
           buffer = buffer,
           density = df$dens[j])
  
  ## Remove row names
  rownames(las.mets) <- c()
  
  ## Write out results
  write.csv(las.mets,df$out.name[j])
}

## Repeat for original data
foreach(a=1:length(plot.lst))%dopar%{
  library(lidR)
  library(tidyverse)
  source("./Code/lidar_decimation_functions.R")
  out.file <- gsub("\\.las","_org_dens_metrics.csv",
                   gsub("Buffers_new","Lidar_metrics/Full",plot.lst[a]))
  x = str_split(plot.lst[a],"_")[[1]][6]
  y = str_split(plot.lst[a],"_")[[1]][7]
  buffer = gsub("\\.las","",
                gsub(".*(?<=\\D)(\\d)","\\1",plot.lst[a],perl=TRUE))
  las <- readLAS(plot.lst[a],select="xyz")
  dens <- density(las)
  tryCatch({
    las.proc(lidar=plot.lst[a]) %>% 
      mutate(density = dens,
             x = x,
             y = y,
             buffer = buffer) %>% 
      write.csv(.,out.file)}, error=function(e){print(e)})
}

