##Sensitivity analysis for lidar metrics - do for one plot at different sizes
##random decimation of n=c(10,25,50,100,500,750,1000) randomizations
source("./Code/00_lidar_decimation_functions.R")

library(lidR)
library(tidyverse)
library(raster)
library(doParallel)
library(viridis)
library(patchwork)

## Read in plot list
plot.lst <- read.csv("./Lidar/plot_list.csv")$x

## Select a random plot from each site (this was based on a random number generator)
sens.lst <- plot.lst[c(13:16,305:308)]

## Make vectors to store the number of repetitions and point densities
reps <- c(10,25,50,100,250,500,750,1000)
dens <- c(175,150,125,100,75,50,40,25,20,15,10,5,4,3,2,1)

## Repeat the file names, densities, and repetitions 16*8*8 times
files <- rep(rep(sens.lst,each=16),8)
dens.lst <- rep(rep(dens,8),8)
reps.lst <- rep(rep(reps,each=16),each=8)

## Create a data frame with files names, point densities, and number of repetitions
df <- data.frame(las = files, dens = dens.lst, reps = reps.lst)

## Register parallel cores
nc <- detectCores()-1
registerDoParallel(cores=nc)

foreach(i=1:nrow(df))%dopar%{
  ## Read in libraries used for this code
  library(lidR)
  library(tidyverse)
  
  ## Refer to code for functions
  source("./Code/lidar_decimation_functions.R")
  
  ## Extract numbers from the lidar file name
  nums <- str_extract_all(df$las[i],"\\d+")
  
  ## Get x and y coordinates, buffers from the numbers
  x = paste0(nums[[1]][3],".",nums[[1]][4])
  y = nums[[1]][5]
  buffer = nums[[1]][6]  
  
  ## Extract file name to create out file name
  file.name <- gsub("\\.las","",
                    gsub("./Lidar/Buffers_new/","",df$las[i]))
  out.file <- paste0("./Lidar/SEM/",file.name,
                     "_density_",df$dens[i],"_reps_",df$reps[i],".csv")
  
  ## Perfor las sensitivity algorithm (repeats the lidar decimation code for the
  ## set number of iterations)
  tryCatch({las.sensitivity(df$las[i],df$dens[i],reps=df$reps[i]) %>% 
      ## Add the density, buffer, and x and y coordinates to the data frame
      mutate(density = df$dens[i],
             x = x,
             y = y,
             buffer = buffer) %>% 
      ## Write out the results
      write.csv(.,out.file)}, error=function(e){print(e)})}


## Code to calculate SEM for the different metrics/different densities/different reps

## List the files generated for the SEM analysis
sem.files <- list.files("./Lidar/SEM/",pattern="*csv",full.names=TRUE)

## List the unique repetitions to use with str_detect
rep.lst <- unique(paste0("_reps_",reps.lst,"\\.csv$"))

## Initiate a data frame
dat <- data.frame()

for(i in 1:length(rep.lst)){
  ## List the files for a specific number of repetitions
  sdm.lst <- sem.files[str_detect(sem.files,rep.lst[i])]
  
  for(j in 1:length(sdm.lst)){
    ## Extract the site name
    site <- str_extract(sdm.lst[j],"krc|assa")
    
    ## Add the site name and number of repetitions to the data
    f <- read.csv(sdm.lst[j]) %>% 
      mutate(site = site, nreps = reps[i])
    
    ## Add to the data frame
    dat <- rbind(dat,f)
  }
}

## Calculate SEM
sem.summ <- dat %>% 
  group_by(nreps=as.factor(nreps),density=as.factor(density),buffer=as.factor(buffer),
           site=site) %>% 
  summarize(across(zmax:mhorshanlad,sem),
            .groups="drop") %>%
  distinct() %>% 
  pivot_longer(cols = -c(density,nreps,site,buffer), names_to = "variable",values_to="SEM") %>% 
  mutate(nreps = as.numeric(as.character(nreps)))

## Write out the results from the SEM calculation
write.csv(sem.summ,"./Lidar/SEM/Collated/sem_summary.csv")
