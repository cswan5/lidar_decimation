## These are the base functions for decimating the lidar and quantifying forest
## structure metrics from the lidar. I will indicate pieces of code that I used
## from other packages with included citations.

## The code between the two lines of hashtags comes from the leafR pacakage. 
## Danilo Roberti Alves de Almeida, Scott Christopher Stark, Carlos Alberto Silva, 
## Caio Hamamura and Ruben Valbuena (2021). leafR: Calculates the Leaf Area Index 
## (LAD) and Other Related Functions. R package version 0.3.5. 
## https://CRAN.R-project.org/package=leafR
################################################################################

#' Count number of points in each Z slice
#'
#' @param Z numeric vector. The heights vector.
#' @param maxZ numeric. The maximum height expected in the whole dataset.
#'
#' @return A [`list`][base::list] of point counts in each Z slice of 1 meter
#'
#' @importFrom data.table data.table
#' @importFrom stats aggregate
#' @export
pointsByZSlice = function(Z, maxZ){
  heightSlices = as.integer(Z) # Round down
  zSlice = data.table::data.table(Z=Z, heightSlices=heightSlices) # Create a data.table (Z, slices))
  sliceCount = stats::aggregate(list(V1=Z), list(heightSlices=heightSlices), length) # Count number of returns by slice
  
  ##############################################
  # Add columns to equalize number of columns
  ##############################################
  colRange = 0:maxZ
  addToList = setdiff(colRange, sliceCount$heightSlices)
  n = length(addToList)
  if (n > 0) {
    bindDt = data.frame(heightSlices = addToList, V1=integer(n))
    sliceCount = rbind(sliceCount, bindDt)
    # Order by height
    sliceCount = sliceCount[order(sliceCount$heightSlices),]
  }
  
  colNames = as.character(sliceCount$heightSlices)
  colNames[1] = "ground_0_1m"
  colNames[-1] = paste0("pulses_", colNames[-1], "_", sliceCount$heightSlices[-1]+1, "m")
  metrics = list()
  metrics[colNames] = sliceCount$V1
  
  return(metrics)
  
} #end function pointsByZSlice

#' Creates a data frame of the 3D voxels information (xyz)
#' with Leaf Area Density values from las file
#'
#' @param normlas.file normalized las file
#' @param grain.size horizontal resolution (suggested 1 meter for lad profiles and 10 meters for LAI maps)
#' @param k coefficient to transform effective LAI to real LAI (k = 1; for effective LAI)
#'
#' @return A [`data.frame`][base::data.frame] of the 3D voxels information (xyz) with Leaf Area Density values
#'
#' @note The values of LAD are not estimated below 1 meter. For the following reasons:
#' ground points influence
#' realtive low sampling
#'
#' @examples
#' # Get the example laz file
#' normlas.file = system.file("extdata", "lidar_example.laz", package="leafR")
#'
#' VOXELS_LAD = lad.voxels(normlas.file,
#'                         grain.size = 2, k=1)
#'
#' @importFrom raster values
#' @importFrom lidR grid_metrics readLAS
#' @importFrom sp coordinates
#' @importFrom stats formula
#' @export
lad.vox = function(las, grain.size = 1, k = 1){
  
  #empty list object that will be fueling with binneds data.frames
  LAD_VOXELS = list()
  Z = NA

  las@data$Z[las@data$Z < 0] = 0
  
  maxZ = floor(max(las@data$Z))
  
  func = formula(paste0("~pointsByZSlice(Z, ", maxZ, ")"))
  t.binneds    = lidR::grid_metrics(las, func, res = grain.size,
                                    start = c(min(las@data$X), max(las@data$Y)))
  t.binneds    = data.frame(sp::coordinates(t.binneds), raster::values(t.binneds))
  names(t.binneds)[1:2] = c("X", "Y")
  
  
  #getting the coordinates X and Y
  #t.binneds$X = coordinates(t.binneds)[,1]
  #t.binneds$Y = coordinates(t.binneds)[,2]
  #t.binneds = as.data.frame(t.binneds) #transforming in a data.frame
  
  #clip product by las files limits
  #t.binneds = t.binneds[t.binneds$X < xmax(.las) &
  #                        t.binneds$X > xmin(.las) &
  #                        t.binneds$Y > ymin(.las) &
  #                        t.binneds$Y < ymax(.las),]
  
  
  #select ground returns
  ground.returns = t.binneds[, grep("ground", names(t.binneds))]
  
  #select columns vegetation above 1m:
  if(nrow(t.binneds) != 1){ #this if is necessary when grain size is the whole plot
    pulses.profile.dz1 = t.binneds[, c(grep("pulses", names(t.binneds)))]
  }else{
    pulses.profile.dz1 = data.frame(matrix(as.numeric(as.character(t.binneds[, c(grep("pulses", names(t.binneds)))])), ncol = length(grep("pulses", names(t.binneds)))))
    names(pulses.profile.dz1) = names(t.binneds)[c(grep("pulses", names(t.binneds)))]
  }
  
  #invert data.frames for the sky be first
  pulses.profile.dz1 = pulses.profile.dz1[,length(pulses.profile.dz1):1] #invert columns
  
  #add grounds returns (0-1m)
  pulses.profile.dz1 = cbind(pulses.profile.dz1, ground.returns)
  rm(ground.returns)
  
  ### total matriz and cumsum.matrix:
  total.pulses.matrix.dz1 = matrix(apply(pulses.profile.dz1, 1, sum), ncol = length(pulses.profile.dz1), nrow = nrow(pulses.profile.dz1))
  cumsum.matrix.dz1 = matrix(apply(pulses.profile.dz1, 1, cumsum), ncol = length(pulses.profile.dz1), nrow = nrow(pulses.profile.dz1), byrow = TRUE)
  
  rm(pulses.profile.dz1)
  
  #Pulses out for each voxel
  pulse.out.dz1 = total.pulses.matrix.dz1 - cumsum.matrix.dz1
  
  #The pulses.out of voxel 1 is the pulses.in of voxel 2 and so on...
  #Therefore, pulse.in is pulse.out without the last line and adding in the
  #first line the total pulses:
  if(nrow(t.binneds) != 1){ #if used when grain size of the whole plot
    pulse.in.dz1 <- cbind(total.pulses.matrix.dz1[,1], pulse.out.dz1[,-c(ncol(pulse.out.dz1))])
  }else{
    pulse.in.dz1 <- c(total.pulses.matrix.dz1[,1], pulse.out.dz1[,-c(ncol(pulse.out.dz1))])
  } #enf if
  
  rm(total.pulses.matrix.dz1, cumsum.matrix.dz1)
  
  # MacArthur-Horn eqquation
  # LAD = ln(S_bottom/S_top)*(1/(dz*K))
  #k value for LAD equation
  dz = 1
  
  LAD.dz1 = log(pulse.in.dz1/pulse.out.dz1) * 1/k * 1/dz
  
  rm(pulse.in.dz1, pulse.out.dz1)
  
  # Remove infinite and NaN values
  #Inf ocorre qndo pulses.out eh zero
  #NaN ocorre qndo pulses.in eh zero
  LAD.dz1[is.infinite(LAD.dz1)] <- NA; LAD.dz1[is.nan(LAD.dz1)] <- NA;
  
  #remove the first 1 meter close to the ground (and the ground too)
  LAD.dz1 = LAD.dz1[, -c(ncol(LAD.dz1))]
  
  #fuel list object
  LAD_VOXELS[["LAD"]] = LAD.dz1
  LAD_VOXELS[["coordenates"]] = t.binneds[,c("X", "Y")]
  
  rm(LAD.dz1, t.binneds)
  
  return(LAD_VOXELS)
}#End function

################################################################################
## End of code from leafR ######################################################

## Function to calculate the mode value of the lidar for each plot
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

## Function to get lidar metrics. Function takes in the following metrics from the
## lidar point cloud: z = elevation, rn = return number, th = threshold (for delineating
## canopy vs non-canopy), class = classification of the point (in this case, ground
## or unclassified)
m <- function(z,rn,th,class){
  first <- rn == 1L ## Define first return
  zfirst <- z[first] ## Get z values for first returns
  nfirst <- length(zfirst) ## Count number of first returns
  nzfirst_canopy <- sum(zfirst>th) ## Sum the number of first returns that are greater thant the canopy threshold value
  probs <- c(0.25,0.5,0.75,0.9,0.95,0.99) ## Set up probability breakpoints
  zq 	  <- as.list(stats::quantile(z, probs)) ## Get z values for each probability quantile
  names(zq) <- paste0("zq", probs*100) ## Associate the z value with its appropriate quantile
  n <- length(z) ## Total number of returns
  metrics <- list(
    zmax <- max(z), ## Maximum elevation
    zmean <- mean(z), ## Mean elevation
    zmode = Mode(z), ## Mode elevation
    zskew <- (sum((z - zmean)^3)/n)/(sum((z - zmean)^2)/n)^(3/2), ##skewness of elevation distribution
    zkurt <-  n * sum((z - zmean)^4)/(sum((z - zmean)^2)^2), ## kurtosis of elevation distribution
    zsd <- stats::sd(z), ## standard deviation of elevations
    zcov <- zsd/zmean, ## covariance of elevations
    zvd <- entropy(z), ## vertical diversity, calculated from the entropy function from the lidR package
    cover <- sum(z>th)/n, ## canopy cover
    cdens <- nzfirst_canopy/nfirst, ## canopy density
    lowveg <- sum(z >= 1 & z <= th)/n, ## low vegetation (between 1 m and canopy threshold)
    pground <- sum(z < 1)/n*100) ## percent of ground points (less than 1 m)
  names(metrics) <- c("zmax","zmean","zmode","zskew","zkurt","zsd","zcov","zvd",
                      "cover","cdens","lowveg","pground") ## rename the metrics list
  
  return(c(metrics,zq)) ## returns quantiles and metrics
  
}

## Shannon LAD  and some related calculations from the leafR package ###########
## Danilo Roberti Alves de Almeida, Scott Christopher Stark, Carlos Alberto Silva, 
## Caio Hamamura and Ruben Valbuena (2021). leafR: Calculates the Leaf Area Index 
## (LAD) and Other Related Functions. R package version 0.3.5. 
## https://CRAN.R-project.org/package=leafR
################################################################################

shan.lad <- function(VOXELS_LAD){
  lad.df <- as.data.frame(VOXELS_LAD$LAD)
  r <- range(lad.df,na.rm=TRUE)
  bin_size <- (max(r) - min(r))/50
  bins <- seq(from=min(r),to=max(r),by=bin_size)
  
  ##Calculations for overall LAD diversity
  lad.cuts <- cut(VOXELS_LAD$LAD,bins)
  lad.cuts <- lad.cuts[!is.na(lad.cuts)]
  lad.summ <- data.frame(lad.bins = summary(lad.cuts)) %>% 
    filter(lad.bins >  0) %>% 
    mutate(lad.prop = lad.bins/sum(lad.bins))
  slad <- -sum(lad.summ$lad.prop*log(lad.summ$lad.prop))
  
  ## Calculation for vertical Shannon LAD diversity (shanladij)
  ladij.cuts <- list()
  for(i in 1:nrow(lad.df)){
    ladij.cuts[[i]] <- cut(as.numeric(lad.df[i,]),bins)
  }
  ladij.cuts <- ladij.cuts[!sapply(ladij.cuts,is.null)]
  ladij.cuts <- lapply(ladij.cuts, function(x) x[!is.na(x)])
  ladij.summ <- lapply(ladij.cuts,summary)
  ladij.summ <- lapply(ladij.summ,data.frame)
  ladij.summ <- lapply(ladij.summ, function(x) x %>% filter(.!=0))
  ladij.summ <- lapply(ladij.summ, function(x) x/sum(x))
  ladij.summ <- ladij.summ[sapply(ladij.summ,nrow)>0]
  shanlad.ij <- data.frame()
  for(i in 1:length(ladij.summ)){
    shanlad.ij[i,1] <- -sum(ladij.summ[[i]]*log(ladij.summ[[i]]))
  }
  
  ## Calculation for horizontal Shannon LAD diversity (shanladh)
  ladh.cuts <- list()
  for(i in 1:ncol(lad.df)){
    ladh.cuts[[i]] <- cut(as.numeric(lad.df[,i]),bins)
  }
  
  ladh.cuts <- ladh.cuts[!sapply(ladh.cuts,is.null)]
  ladh.cuts <- lapply(ladh.cuts, function(x) x[!is.na(x)])
  ladh.summ <- lapply(ladh.cuts,summary)
  ladh.summ <- lapply(ladh.summ,data.frame)
  ladh.summ <- lapply(ladh.summ, function(x) x %>% filter(.!=0))
  ladh.summ <- lapply(ladh.summ, function(x) x/sum(x))
  ladh.summ <- ladh.summ[sapply(ladh.summ,nrow)>0]
  shanlad.h <- data.frame()
  for(i in 1:length(ladh.summ)){
    shanlad.h[1,i] <- -sum(ladh.summ[[i]]*log(ladh.summ[[i]]))
  }
  
  shan.div <- list(shanlad = slad,
                   vershanlad = shanlad.ij,
                   horshanlad = shanlad.h)
  
  return(shan.div)
}

lad.profile = function(VOXELS_LAD, relative = FALSE){
  
  if(relative == TRUE){
    t.lad.profile = apply(VOXELS_LAD$LAD, 2, mean, na.rm = TRUE)
    t.lad.profile = t.lad.profile/sum(t.lad.profile)*100
  }else{
    t.lad.profile = apply(VOXELS_LAD$LAD, 2, mean, na.rm = TRUE)
  }
  
  max_height = ncol(VOXELS_LAD[[1]]) + .5
  
  t.lad.profile = data.frame(height = seq(1.5, max_height), lad = t.lad.profile[length(t.lad.profile):1])
  
  return(t.lad.profile)
  
}#end looping

FHD = function(lad_profile, evenness = TRUE, LAD.threshold = -1){
  
  # applying threshold
  if(LAD.threshold == -1) LAD.threshold <- 1 / length(lad_profile$height)
  lad_profile <- lad_profile[lad_profile$lad >= LAD.threshold,]
  
  lad_profile$lad <- lad_profile$lad / 100
  
  #calculating FHD
  if(evenness){
    FHD = - sum( lad_profile$lad * log(lad_profile$lad) ) / log( length(lad_profile$height) )
  }else{
    FHD = - sum( lad_profile$lad * log(lad_profile$lad) )
  } #end if else
  
  return(FHD)
  
} #end function

GC = function(las, threshold = 1){
  
  las <- filter_poi(las, Z > threshold)
  
  # calculate Gini
  n <- length(las@data$Z)
  x <- sort(las@data$Z)
  G <- 2 * sum(x * 1L:n)/sum(x) - (n + 1L)
  GC <- G/(n - 1L)
  
  return(GC)
  
} #end function

################################################################################
## End of functions from leafR #################################################

## Function to calculate voxel-level metrics from the point cloud where las = lidar
## point cloud originated from readLAS
m.vox <- function(las){
  las <- las[las$Z>=0.5] ## Filter points greater than 0.5 m in elevation
  lad.voxels <- lad.vox(las) ## perform lad.vox function (above)
  lad.prof <- lad.profile(lad.voxels,relative=TRUE) ## perform lad profile function (above)
  lad <- as.data.frame(lad.voxels$LAD) ## change the lad.voxels output to a data frame
  
  ## Calculate the row means of the LAD (calculate LAD mean for each x,y coordinate,
  ## which is the mean of each voxel column (see explainer image in README))
  mlad.ij <- lad %>% summarize(mladij=rowMeans(.,na.rm=TRUE)) %>%  
    na.omit(.)
  
  lad.ij.bar <- mean(mlad.ij[,1]) ## Calculate mean of all LAD, ij
  cvmlad <- sqrt(sum((mlad.ij-lad.ij.bar)^2)/(nrow(mlad.ij)-1))/lad.ij.bar ## Coefficient of variance of mean LAD ij
  
  ## Calculate the mean LAD across horizontal slices of the canopy (see explainer
  ## image in README)
  lad.h.bar <- lad %>% 
    summarize(mladh = colMeans(.,na.rm=TRUE))
  
  ## Initialize a data frame to store results
  LADdiffijh <- data.frame()
  
  ## Calculate the square of the difference between the LAD for a voxel and the mean LAD ij
  for(i in 1:nrow(mlad.ij)){
    for(h in 1:nrow(lad.h.bar)){
      LADdiffijh[i,h] <- (lad[i,h]-mlad.ij[i,])^2
    }
  }
  
  sumLADdiff <- LADdiffijh %>% summarize(sumLADdiff=rowSums(.,na.rm=TRUE)) ## Take the sum of all differences calculated above
  cvladij <- sqrt((1/(nrow(lad.h.bar)-1))*sumLADdiff)/lad.ij.bar ## Calculate the coefficient of variation of LAD ij
 
  ## Initialize a data frame to store the results
  LADdiffhij <- data.frame()
  
  ## Calculate the difference between LAD for each voxel in a horizontal slice
  ## and the mean LAD of all voxels in the same horizontal slice
  for(h in 1:nrow(lad.h.bar)){
    for(i in 1:nrow(mlad.ij)){
      LADdiffhij[h,i] <- (lad[i,h]-lad.h.bar[h,])^2
    }
  }
  
  ## Take the sum of the differences calculated above
  sumLADdiffhij <- LADdiffhij %>% summarize(sumLADdiff=rowSums(.,na.rm=TRUE))
  
  horcvlad <- sqrt((1/(nrow(mlad.ij)-1))*sumLADdiffhij)/lad.h.bar ## Calculate the coefficient of variation of LADh
  
  lad.divs <- shan.lad(lad.voxels) ## Calculate Shannon diversity
  slad <- lad.divs[[1]] ## Retrieve Shannon diversity from the results of the diversity calculation
  mverslad <- mean(lad.divs[[2]][,1]) ## Take the mean of the Shannon diversity across columns of voxels
  mhorslad <- rowMeans(lad.divs[[3]]) ## Take the mean of the shannon diversity across horizontal slices
  
  ## Save voxel metrics to a list and rename the list elements with their associated
  ## metric names
  vox.metrics <- list(
    mlad <- mean(lad.voxels$LAD,na.rm=TRUE),
    cvmlad <- cvmlad,
    mcvlad <- sum(cvladij,na.rm=TRUE)/nrow(mlad.ij),
    mhorcvlad <- sum(horcvlad,na.rm=TRUE)/nrow(lad.h.bar),
    fhd <- FHD(lad.prof),
    slad,
    mverslad,
    mhorslad
  )
 names(vox.metrics)<-c("mlad","cvmlad","mcvlad","mhorcvlad","fhd","shanlad","mvershanlad","mhorshanlad") 
 return(vox.metrics)
}

## Function for the sensitivity analysis (to determine the number of repetitions
## needed to minimize error from random decimation). lidar = name of the point cloud
## file, pts = point density, reps = number of times to randomly decimate the lidar
las.sensitivity <- function(lidar,pts,reps){
  all.met <- data.frame() ## Empty data to store the results of the lidar quantification
  las <- readLAS(lidar,select="xyzrc") ## Read in the las data
  i <- 1 ## Initialize the value for i (superfluous)
  for(i in 1:reps){
    las.dec <- decimate_points(las,random(pts)) ## use the decimate points function from the lidar package to perform the random decimation
    chm <- as.data.frame(rasterize_canopy(las.dec,res=0.5)) ## create a canopy height model from the decimated lidar
    cld.met <- data.frame(cloud_metrics(las.dec,func=~m(z=Z,rn=ReturnNumber,th=4,
                                                        class=Classification))) ## calculate metrics from the lidar point cloud (from lidR package)
    vox.met <- data.frame(m.vox(las.dec)) ## calculate voxel metrics from the voxels
    cld.met$rugosity <-  sd(chm$Z) ## calculate rugosity by taking the standard deviation of the chm
    cld.met$gini <- GC(las.dec) ## calculate Gini (from leafR function)
    cld.met <- cbind(cld.met,vox.met) ## Collate the point cloud and voxel metrics 
    
    all.met <- rbind(all.met,cld.met) ## Collate the results across all repetitions
  }
  return(all.met)
}

## Function to take in the lidar point cloud. lidar = name of point cloud file
las.proc <- function(lidar){
  all.met <- data.frame() ## Initialize an empty data frame to store results
  las <- readLAS(lidar,select="xyzrc") ## Read in lidar
  chm <- as.data.frame(rasterize_canopy(las,res=0.5)) ## Create a canopy height model
  cld.met <- data.frame(cloud_metrics(las,func=~m(z=Z,rn=ReturnNumber,th=4,
                                                  class=Classification))) ## Calculate the cloud metrics (from lidR package)
  vox.met <- data.frame(m.vox(las)) ## Calculate the voxel metrics
  cld.met$rugosity <-  sd(chm$Z) ## Calculate rugosity by taking the standard deviation of the chm
  cld.met$gini <- GC(las) ## Calculate Gini (from leafR function)
  cld.met <- cbind(cld.met,vox.met) ## Collate the point cloud and voxel metrics
  
  all.met <- rbind(all.met,cld.met) ## Collate all of the metrics for each repetition
  
  return(all.met)
}

sem <- function(x) {sd(x)/sqrt(length(x))} ## Calculate standard error of mean for sensitivity analysis



