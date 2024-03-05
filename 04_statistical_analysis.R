## Code to run statistical analysis for lidar decimation study

## Import libraries
library(doBy)
library(tidyverse)
library(patchwork)
library(reshape2)
library(plotly)
library(Hmisc)
library(vroom)
library(cluster)
library(flextable)
library(car)
library(emmeans)

## This is the code for the ANOVA and Tukeys tests ##

rnames <- c("Density","Residuals")

## Define variables that are being tested
vars <- c("zmax","zmean","zmode","zskew","zkurt","zsd","zcov","cover",
          "cdens","lowveg","pground","zq25","zq50","zq75","zq90","zq95","zq99",
          "gini","zvd","mlad","cvmlad","mhorcvlad","shanlad","mvershanlad",
          "mhorshanlad","rugosity")

## These are the pairwise point densities
levs <- c("1-2","1-3", "1-4","1-5","1-10","1-15","1-20","1-25","1-40","1-50","1-75","1-100","1-125","1-150","1-175",
          "2-3","2-4","2-5","2-10","2-15","2-20","2-25","2-40","2-50","2-75","2-100","2-125","2-150","2-175",
          "3-4","3-5","3-10","3-15","3-20","3-25","3-40","3-50","3-75","3-100","3-125","3-150","3-175",
          "4-5","4-10","4-15","4-20","4-25","4-40","4-50","4-75","4-100","4-125","4-150","4-175",
          "5-10","5-15","5-20","5-25","5-40","5-50","5-75","5-100","5-125","5-150","5-175",
          "10-15","10-20","10-25","10-40","10-50","10-75","10-100","10-125","10-150","10-175",
          "15-20","15-25","15-40","15-50","15-75","15-100","15-125","15-150","15-175",
          "20-25","20-40","20-50","20-75","20-100","20-125","20-150","20-175",
          "25-40","25-50","25-75","25-100","25-125","25-150","25-175",
          "40-50","40-75","40-100","40-125","40-150","40-175",
          "50-75","50-100","50-125","50-150","50-175",
          "75-100","75-125","75-150","75-175",
          "100-125","100-150","100-175",
          "125-150","125-175",
          "150-175")

sizes <- c(10,15,20,25)

## Read in the full datasets for Assateague and KRC
assa <- "./Lidar/Lidar_metrics/Collated/assateague.csv"
krc <- "./Lidar/Lidar_metrics/Collated/krc.csv"

bart <- function(f,site,sizes){
  ## Sets up a function to perform Bartlett's test for determining homoscedasticity
  ## of variance
  
  ## Read in the file
  fi <- read.csv(f)
  
  ## Set up the data frame for the plots to get all of the relevant variables and
  ## only test densities of interest (<= 175); remove fhd and mcvlad from dataset
  ## because these variables were not kept in the final analysis
  if(site == "KRC"){
    ## If site is KRC, remove two bad plots
    fi.sm <- fi %>% dplyr::select(zmax:plot) %>% 
      filter(density <= 175 & plot != 51 & plot != 64) %>% 
      mutate(buffer = as.factor(buffer),
             density = as.factor(density)) %>% 
      select(-c(fhd,mcvlad))
  } else{
    
    fi.sm <- fi %>% dplyr::select(zmax:plot) %>%
      filter(density <= 175) %>% 
      mutate(buffer = as.factor(buffer),
             density = as.factor(density)) %>% 
      select(-c(fhd,mcvlad))
  }
  
  ## Get the column names
  cols <- names(fi.sm)[1:26]
  
  ## Initialize an empty data frame to store the results
  bart.all <- data.frame()
  
  for(k in 1:4){
    ## Select the plot size
    size <- sizes[k]
    
    ## Make a list to store the results for each variable
    barts <- vector(mode="list", length=length(cols))
    
    ## Filter out only the buffer of interst
    fi.filt <- fi.sm %>% filter(buffer == size)
    
    ## Perform Bartlett test and save results in the empty list
    bart.mod <- lapply(cols, function(x) bartlett.test(reformulate(termlabels = 
                                                                     "density",
                                                                   response = x),
                                                       data = fi.filt))
    
    ## Name the list elements the same as the column names (which are the names
    ## of the lidar-calculated variables of interest)
    names(bart.mod) <- cols
    
    ## Keep relevant information from Bartlett test to determine homoscedasticity 
    ## of variance
    bart.df <- bind_rows(lapply(seq_along(bart.mod), function(i) 
    {data.frame(metric = names(bart.mod)[[i]], p_value= round(bart.mod[[i]]$p.value,2), 
                buffer = size)}))
    
    ## Append the data frame to include the results from all buffer sizes
    bart.all <- rbind(bart.all,bart.df)
  }
  
  ## Filter out only those values that were significant (i.e., have heteroscedastic
  ## variance)
  bart.all <- bart.all %>% filter(p_value <= 0.05)
  
  ## Write out the results
  write.csv(bart.all, paste0("./Results/Bartletts/",site,"_bartletts.csv"))
}

## Use the bart function to perform Bartlett test on the data from Assateague
## and KRC
bart(f=assa,site="AINS",sizes=sizes)
bart(krc,"KRC",sizes)

anovs <- function(f,levs,sizes,site){
  ## Function to perform ANOVA on data. ANOVA type is dependent on results of
  ## Bartlett's test
  
  ## Read in the file
  fi <- read.csv(f)
  
  ## Select the variables of interest for point densities <= 175 pts/m2; turn
  ## buffer and density into factors instead of doubles; remove fhd and mcvlad
  ## because they were not part of the final variables we selected
  if(site == "KRC"){
    ## For KRC, two "bad" sites needed to be removed
    fi.sm <- fi %>% dplyr::select(zmax:plot) %>% 
      filter(density <= 175 & plot != 51 & plot != 64) %>% 
      mutate(buffer = as.factor(buffer),
             density = as.factor(density)) %>% 
      select(-c(fhd,mcvlad))
  } else{
    
    fi.sm <- fi %>% dplyr::select(zmax:plot) %>%
      filter(density <= 175) %>% 
      mutate(buffer = as.factor(buffer),
             density = as.factor(density)) %>% 
      select(-c(fhd,mcvlad))
  }
  
  for(k in 1:4){
    ## Select the buffer size from the list of sizes
    size <- sizes[k]
    
    ## Get the column (i.e., variable) names
    cols <- names(fi.sm)[1:26]
    
    ## These are the variables that had heteroscedastic variance across all sites
    white <- c("zmax","zmode","rugosity","mlad","cvmlad","mhorcvlad",
               "shanlad","mvershanlad","mhorshanlad")
    
    ## Add in variables that had heteroscedastic variance at specific buffers/sites
    if(size == 10 & site == "AINS"){
      white <- c(white, "zskew","zkurt","zq25","gini")
    } else if(size == 15 & site == "AINS"){
      white <- c(white, "zkurt","zq25","gini")
    } else if(size == 20 & site == "AINS" | size == 25 & site == "AINS"){
      white <- c(white, "zq25")
    } else if(size == 10 & site == "KRC"){
      white <- c(white, "zkurt","zcov","gini")
    } else if(size == 15 & site == "KRC"){
      white <- c(white,"zvd","gini")
    } else if(size == 20 & site == "KRC"){
      white <- c(white,"zkurt","zvd")
    } else {
      white <- white[!white %in% "zmax"]
    }
    
    ## Define variables as only those that were not in the white list (from Bartlett test)
    cols <- cols[which(!cols%in%white)]
    
    ## Filter out only the data from the specific buffer size
    fi.filt <- fi.sm %>% filter(buffer == size) 
    
    ## Select column names of variables with homoscedastic variance and the point
    ## density column
    f.cols <- fi.filt %>% select(c(cols,density))
    
    ## Create the linear model for the normal ANOVA
    fi.mod <- lapply(cols, function(x) lm(reformulate(termlabels = 
                                                        "density",
                                                      response = x),
                                          data = fi.filt))
    
    ## Name the list elements by the variable names
    names(fi.mod) <- cols
    
    ## Create the lineal models for the white corrected ANOVA
    white.mod <- lapply(white, function(x) lm(reformulate(termlabels = "density",
                                                          response = x),
                                              data = fi.filt))
    
    ## Name the list elements by the variable names
    names(white.mod)<-white
    
    ## Run the ANOVA function on the linear models of the homoscedastic variables
    fi.aov <- lapply(fi.mod, function(x) Anova(x))
    
    ## Run the ANOVA function on the linear models of the heteroscedastic variables
    white.aov <- lapply(white.mod, function(x) Anova(x, white.adjust=TRUE))
    
    ## Put the results of the ANOVA into a data fram
    fi.mod.df <- lapply(seq_along(fi.aov), function(i) data.frame(metric = c(names(fi.aov)[[i]],"residuals"),
                                                                  DF = fi.aov[[i]]$Df,
                                                                  F_value = fi.aov[[i]]$`F value`,
                                                                  p_value = round(fi.aov[[i]]$`Pr(>F)`,3)))
    
    white.mod.df <- lapply(seq_along(white.aov), function(i) data.frame(metric = c(names(white.aov)[[i]],"residuals"),
                                                                        DF = white.aov[[i]]$Df,
                                                                        F_value = white.aov[[i]]$F,
                                                                        p_value = round(white.aov[[i]]$`Pr(>F)`,3)))
    
    ## Put all the ANOVA results into a big data frame
    mods <- bind_rows(fi.mod.df,white.mod.df) %>% 
      filter(p_value <= 0.05) %>% 
      slice(order(factor(metric,levels = vars))) %>% 
      filter(metric != "residuals")
    
    mods.all <- bind_rows(fi.mod.df,white.mod.df) %>% 
      slice(order(factor(metric,levels = vars))) %>% 
      filter(metric != "residuals")
    
    ## Write out ANOVA results
    write.csv(mods,paste0("./Results/ANOVA/",site,size,"_anova.csv"))
    write.csv(mods.all,paste0("./Results/ANOVA/",site,size,"_anova_full.csv"))
    
    ## Filter out only metrics where the ANOVA was significant
    fi.sig <- fi.mod[mods$metric]
    fi.sig <- fi.sig[!sapply(fi.sig,is.null)]
    
    white.sig <- white.mod[mods$metric]
    white.sig <- white.sig[!sapply(white.sig,is.null)]
    
    ## Combine significant results across white-corrected and regular ANOVA
    all.sig <- c(fi.sig,white.sig)
    
    ## Set up the marginals to perform Tukeys
    marginal <- lapply(all.sig, function(x) emmeans(x, ~density))
    
    ## Perform Tukey's to do pairwise comparison between point densities
    tuk.pairs <- lapply(marginal, function(x) pairs(x, adjust = "Tukey", infer=TRUE))
    
    ## Turn the results of the Tukey's test into a tibble
    tk.df <- lapply(tuk.pairs, function(x) as_tibble(x))
    
    ## Filter out significant values from Tukey's pairings and split the point
    ## density values into separate columns
    tk.all <- bind_rows(tk.df,.id="metric") %>% filter(`p.value`<=0.05) %>% 
      mutate(`p.value`=round(`p.value`,3),
             contrast = map_chr(contrast, ~str_replace_all(.x, pattern="density","")),
             dens1 = str_split_fixed(contrast, " - ",2)[,1],
             dens2 = str_split_fixed(contrast, " - ",2)[,2]) 
    
    ## Save the full Tukey's results into a data frame
    tk.full <- bind_rows(tk.df,.id="metric") %>% 
      mutate(`p.value`=round(`p.value`,3),
             contrast = map_chr(contrast, ~str_replace_all(.x, pattern="density","")),
             dens1 = str_split_fixed(contrast, " - ",2)[,1],
             dens2 = str_split_fixed(contrast, " - ",2)[,2]) 
    
    ## Rename the data frames
    names(tk.all) <- c("metric","contrast","diff","SE","DF","lower","upper","t_ratio","p_value","dens1","dens2")
    names(tk.full) <- c("metric","contrast","estimate","SE","df","lower","upper","t_ratio","p_value","dens1","dens2")
    
    ## Create names for the save files
    out.name <- paste0("./Results/Tukey/",site,size,"_tukey.csv")
    full.out <- paste0("./Results/Tukey/",site,size,"_tukey_full.csv")
    
    ## Save the results from the Tukey's tests
    write.csv(tk.all,out.name)
    write.csv(tk.full,full.out)
  }
}

## Run the ANOVA and Tukey's tests for all variables of interest at both sites
anovs(f=krc,levs=levs,sizes=sizes,site="KRC")
anovs(f=assa,levs=levs,sizes=sizes,site="AINS")

##########################################################
##                                                      ##
## Below is the code for the reliability ratio analysis ##
##                                                      ##
##########################################################

## Pull in full results from Assateague and limit to point densities of interest
assa <- as.data.frame(read.csv("./Lidar/Lidar_metrics/Collated/assateague.csv")) %>% 
  dplyr::select(zmax:plot) %>% filter(density <= 175)


fun <- function(x){
  v = var(x)
}

# among-plot variance for Assateague
assa.apv <- summaryBy(zmax + zmean + zmode + zskew + zkurt + zsd + zcov + zvd +
                        cover + cdens + lowveg + pground + zq25 + zq50 + zq75 +
                        zq90 + zq95 + zq99 + rugosity + gini + mlad + cvmlad +
                        mcvlad + mhorcvlad + fhd + shanlad + mvershanlad +
                        mhorshanlad ~ density +buffer, data=assa, FUN=var) # among plot variance

# within-plot variance for Assateague
assa.wpv <- summaryBy(zmax + zmean + zmode + zskew + zkurt + zsd + zcov + zvd +
                        cover + cdens + lowveg + pground + zq25 + zq50 + zq75 +
                        zq90 + zq95 + zq99 + rugosity + gini + mlad + cvmlad +
                        mcvlad + mhorcvlad + fhd + shanlad + mvershanlad +
                        mhorshanlad ~ plot + buffer + density, data=assa, FUN = var) # within plot sample variance

assa.wpv <- summaryBy(zmax.var + zmean.var + zmode.var + zskew.var + zkurt.var + 
                      zsd.var + zcov.var + zvd.var + cover.var + cdens.var + 
                      lowveg.var + pground.var + zq25.var + zq50.var + zq75.var +
                      zq90.var + zq95.var + zq99.var + rugosity.var + gini.var + 
                      mlad.var + cvmlad.var + mcvlad.var + mhorcvlad.var +
                      fhd.var + shanlad.var + mvershanlad.var + mhorshanlad.var ~ 
                      buffer + density, data=assa.wpv, FUN=var) # variance of the within plot sample variance

## Calculate the reliability ratio
assa.rr <- round(assa.apv/(assa.apv+assa.wpv) ,2)

## Make columns for point density and buffer
assa.rr <- assa.rr %>% 
  mutate(density = assa.apv$density, buffer = assa.apv$buffer)


##Repeat for KRC
krc <- as.data.frame(read.csv("./Lidar/Lidar_metrics/Collated2/krc.csv") %>% 
                       filter(site != c(51,64)) %>% 
                       dplyr::select(zmax:plot) %>%
                       filter(density <= 175))


# among-plot variance of X 
krc.apv <- summaryBy(zmax + zmean + zmode + zskew + zkurt + zsd + zcov + zvd +
                       cover + cdens + lowveg + pground + zq25 + zq50 + zq75 +
                       zq90 + zq95 + zq99 + rugosity + gini + mlad + cvmlad +
                       mcvlad + mhorcvlad + fhd + shanlad + mvershanlad +
                       mhorshanlad ~ density + buffer, data=krc, FUN=var) # among plot variance

# within-plot variance of X
krc.wpv <- summaryBy(zmax + zmean + zmode + zskew + zkurt + zsd + zcov + zvd +
                       cover + cdens + lowveg + pground + zq25 + zq50 + zq75 +
                       zq90 + zq95 + zq99 + rugosity + gini + mlad + cvmlad +
                       mcvlad + mhorcvlad + fhd + shanlad + mvershanlad +
                       mhorshanlad ~ plot + buffer + density, data=krc, FUN = var) # within plot sample variance

krc.wpv <- summaryBy(zmax.var + zmean.var + zmode.var + zskew.var + zkurt.var + 
                       zsd.var + zcov.var + zvd.var + cover.var + cdens.var + 
                       lowveg.var + pground.var + zq25.var + zq50.var + zq75.var +
                       zq90.var + zq95.var + zq99.var + rugosity.var + gini.var + 
                       mlad.var + cvmlad.var + mcvlad.var + mhorcvlad.var +
                       fhd.var + shanlad.var + mvershanlad.var + mhorshanlad.var ~ 
                       buffer + density, data=krc.wpv, FUN=var) # variance of the within plot sample variance

krc.rr <- krc.rr %>% 
  mutate(density = krc.apv$density, buffer = krc.apv$buffer)


## Write out reliability ratio results for Assateague and KRC
write.csv(assa.rr,"./Results/Reliability_ratio/assateague_rr.csv")
write.csv(krc.rr,"./Results/Reliability_ratio/krc_rr.csv")

####################################################
##                                                ##
## Below is the code for the correlation analysis ##
##                                                ##
####################################################

## Define the density, buffers, and metrics
dens <- c(1,2,3,4,5,10,15,20,25,40,50,75,100,125,150,175)
buff <- c(10,15,20,25)
mets <- c("zmax","zmean","zmode","zskew","zkurt","zsd","zcov","cover",
          "cdens","lowveg","pground","zq25","zq50","zq75","zq90","zq95","zq99",
          "gini","zvd","mlad","cvmlad","mhorcvlad","shanlad","mvershanlad",
          "mhorshanlad","rugosity")


## Make empty data frames to store the results
assa.dats <- data.frame()
krc.dats <- data.frame()

for(i in 1:length(buff)){
  ## Loop through each buffer
  b <- buff[i]
  for(j in 1:length(dens)){
    ## Loop through each point density
    d <- dens[j]
    
    ## Filter out the buffer and point density of interest from the Assateague
    ## dataset
    assa.df <- assa %>% 
      filter(buffer==b & density==d) %>% 
      select(zmax:mhorshanlad)
    
    ## Perform a Pearson's r correlation between the variables for that buffer 
    ## size and point density
    assa.cor <- rcorr(as.matrix(assa.df),type="pearson")
    
    ## Keep the results from the Pearson's r correlation test
    assa.cor.df <- as.data.frame(assa.cor$r) %>% 
      mutate(var1=rownames(assa.cor$r)) %>% 
      pivot_longer(cols=zmax:mhorshanlad,names_to="var2",values_to="r_corr") %>% 
      mutate(r_corr=round(r_corr,2))
    
    ## Get the p-value for the correlation test
    assa.p <- as.data.frame(assa.cor$P) %>% 
      mutate(var1=rownames(assa.cor$P)) %>% 
      pivot_longer(cols=zmax:mhorshanlad,names_to="var2",values_to="p_value") %>% 
      mutate(p_value=round(p_value,2))
    
    ## Put the correlation results and signficance value together for plotting
    ## the data later
    assa.plot.data <- full_join(assa.cor.df,assa.p,join_by(var1,var2)) %>% 
      mutate(density=d,
             buffer=b,
             site="AINS")
    
    ## Put all of the Assateague data together
    assa.dats <- rbind(assa.dats,assa.plot.data)
    
    ## Perform the same steps for KRC
    krc.df <- krc %>% 
      filter(buffer==b & density==d) %>% 
      select(zmax:mhorshanlad)
    
    krc.cor <- rcorr(as.matrix(krc.df),type="pearson")
    
    krc.cor.df <- as.data.frame(krc.cor$r) %>% 
      mutate(var1=rownames(krc.cor$r)) %>% 
      pivot_longer(cols=zmax:mhorshanlad,names_to="var2",values_to="r_corr") %>% 
      mutate(r_corr=round(r_corr,2))
    
    krc.p <- as.data.frame(krc.cor$P) %>% 
      mutate(var1=rownames(krc.cor$P)) %>% 
      pivot_longer(cols=zmax:mhorshanlad,names_to="var2",values_to="p_value") %>% 
      mutate(p_value=round(p_value,2))
    
    krc.plot.data <- full_join(krc.cor.df,krc.p,join_by(var1,var2)) %>% 
      mutate(density=d,
             buffer=b,
             site="KRC")
    
    krc.dats <- rbind(krc.dats,krc.plot.data)
    
  }
}

## Put the Assateague and KRC datasets together
dats <- bind_rows(assa.dats,krc.dats)

## Reorder the metrics
dats$var1  <- fct_relevel(dats$var1,mets)
dats$var2 <- fct_relevel(dats$var2,mets)

## Write out the results of the correlation test
write.csv(dats,"./Results/Correlations/correlations.csv")
