## Code to generate images for paper
## March 21, 2023
## Christine Swanson

## Import necessary libraries
library(gtable)
library(cowplot)
library(tidyverse)
library(grid)
library(plotly)
library(lidR)
library(viridis)
library(ggplotify)
library(patchwork)
library(vroom)


## Function to move legend into empty grob
shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}
## Factor relevel setup (to get the variables in the desired order)
vars <- c("zmax","zmean","zmode","zskew","zkurt","zsd","zcov","cover",
          "cdens","lowveg","pground","zq25","zq50","zq75","zq90","zq95","zq99",
          "gini","zvd","mlad","cvmlad","mhorcvlad","shanlad","mvershanlad",
          "mhorshanlad","rugosity")

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

## Metrics derived from the full point cloud
cloud <- c("zmax","zmean","zmode","zskew","zkurt","zsd","zcov","cover","cdens",
           "lowveg","pground","zq25","zq50","zq75","zq80","zq90","zq95","zq99","gini")

## Metrics derived from voxelized point cloud
voxel <- c("zvd","mlad","cvmlad","mhorcvlad","shanlad","mvershanlad",
           "mhorshanlad")

var.name <- as_labeller(c("zmax" = "Max height","zmean" = "Mean height", 
                          "zmode" = "Mode height", "zskew" = "Height skewness",
                          "zkurt" = "Height kurtosis", "zsd" = "SD of height",
                          "zcov" = "COV of height", "zvd" = "Vertical diversity",
                          "cover" = "Canopy cover", "cdens" = "Canopy density",
                          "lowveg" = "% understory", "pground" = "% ground",
                          "zq25" = "25th quartile", "zq50" = "50th quartile",
                          "zq75" = "75th quartile", "zq90" = "90th quantile",
                          "zq95" = "95th quantile", "zq99" = "99th quantile",
                          "rugosity" = "Rugosity", "gini" = "Gini coefficient",
                          "mlad" = "Mean LAD", "mcvlad" = "COV of mean LAD",
                          "mhorcvlad" = "COV of mean horizontal LAD",
                          "fhd" = "FHD"))

buffs <- c(10,15,20,25)

titles <- c("(A) 10x10 m", "(B) 15x15 m", "(C) 20x20 m","(D) 25x25 m")

## SEM 
## Read in results from SEM analysis and filter out variables that were not kept
## in final analysis
sem.summ <- read.csv("./Lidar/SEM/Collated/sem_summary.csv") %>% 
  filter(variable != "fhd" & variable != "mcvlad")

## Reorder the variables
sem.summ$variable<-fct_relevel(sem.summ$variable,vars)

## Plot for results from Assateague
p1 <- sem.summ %>% filter(site=="assa" & buffer == 25) %>% 
  ggplot(aes(x=nreps,y=SEM,color=as.factor(density)))+geom_line(aes(group=as.factor(density)))+
  geom_vline(aes(xintercept=250))+
  facet_wrap(~variable,scales="free")+
  scale_color_viridis(discrete=TRUE,option="C","Point density")+
  labs(title="(A) AINS",x="Number of iterations")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

## Plot for results from KRC
p2 <- sem.summ %>% filter(site=="krc" & buffer == 25) %>% 
  ggplot(aes(x=nreps,y=SEM,color=as.factor(density)))+geom_line(aes(group=as.factor(density)))+
  geom_vline(aes(xintercept=250))+
  facet_wrap(~variable,scales="free")+
  scale_color_viridis(discrete=TRUE,option="C","Point density")+
  labs(title="(B) KRC",x="Number of iterations")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

## Plot layout
p1/p2 +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom',
                                          legend.direction = 'horizontal') &
  guides(color=guide_legend(nrow=1))

## Sav figures
ggsave("./Figures/SEM/sem_25.png",height=15,width=11,dpi=600)

## Tukey figures
## Read in all of the results from Tukey's test
tuks <- list.files("./Results/Tukey/",pattern="._full\\.csv")

## Create and empty data frame
tuks.df <- data.frame()

for(i in 1:length(tuks)){
  df <- read.csv(paste0("./Results/Tukey/",tuks[i])) %>% ## Read in the results
    mutate(Site = toupper(word(tuks[i],1,sep="[0-9]")), ## Site name from the file name
           Buffer = str_extract(tuks[i],"[0-9]{2}"), ## Extract the buffer value from the file name
           p_value = ifelse(p_value > 0.05,NA,p_value), ## Change non-significant p-values to NA
           estimate = ifelse(is.na(p_value),NA,estimate), ## Change estimate values to NA if non-significant
           lower = ifelse(is.na(p_value),NA,lower), ## Change lower limit of confidence interval to NA if non-significant
           upper = ifelse(is.na(p_value),NA,upper), ## Change upper limit of confidence interval to NA if non-significant
           contrast = str_replace_all(contrast," ","")) %>% ## Remove spaces in "contrast" column
    filter(metric != "fhd" & metric != "mcvlad") ## Remove metrics that were not used in study
  
  tuks.df <- rbind(tuks.df,df) ## Append to data frame of full results
}

tuks.df$metric <- fct_relevel(tuks.df$metric,vars) ## Reorder the metrics 

## Group p-values into categories
tuks.df <- tuks.df %>% 
  mutate(p_cat = case_when(p_value <= 0.05 & p_value > 0.01 ~ "<=0.05",
                           p_value <= 0.01 & p_value > 0.001 ~ "<= 0.01",
                           p_value <= 0.001 ~ "<= 0.001"))

## Filter out metrics derived from the full point cloud in a 25 m x 25 m buffer
tuk.cld <- tuks.df %>% 
  filter(Buffer==25&metric%in%cloud&!is.na(p_cat))

## Plot the results from Tukey's for the cloud metrics
a<-ggplot(tuk.cld,aes(x=as.factor(dens1),y=as.factor(dens2),
                        fill=p_cat,color="white"))+
  geom_tile(na.rm=TRUE)+
  facet_grid(rows=vars(Site),cols=vars(metric))+
  scale_fill_manual(values=c("#0d0887","#aa2395","#fdc527"),na.value="white",
                    name = "p-value",
                    labels=c(expression(""<=0.001),expression(""<=0.01),
                             expression(""<=0.05)))+
  scale_color_manual(values=c("white"),
                     name = "",
                     labels = "")+
  scale_x_discrete(limits=c("1","2","3","4","5","10","15","20","25","40",
                            "50","75","100","125","150"))+
  scale_y_discrete(limits=c("2","3","4","5","10","15","20","25","40",
                            "50","75","100","125","150","175"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")+
  labs(x="Point density",y="Point density") +
  guides(color="none")

a

## Save plot
ggsave("./Figures/Tukey/cloud_25.png",height=6,width=7,dpi=600)

## Filter metrics calculated from the voxels
tuk.vox <- tuks.df %>% 
  filter(Buffer==25&metric%in%voxel&!is.na(p_cat))

## Plot results from Tukey's for the voxel metrics
b<-ggplot(tuk.vox,aes(x=as.factor(dens1),y=as.factor(dens2),
                       fill=p_cat,color="white"))+
  geom_tile(na.rm=TRUE)+
  facet_grid(rows=vars(Site),cols = vars(metric))+
  scale_fill_manual(values=c("#0d0887","#aa2395","#fdc527"),na.value="white",
                    name = "p-value",
                    labels=c(expression(""<=0.001),expression(""<=0.01),
                             expression(""<=0.05)))+
  scale_color_manual(values=c("white"),
                     name = "",
                     labels = "")+
  scale_x_discrete(limits=c("1","2","3","4","5","10","15","20","25","40",
                            "50","75","100","125","150"))+
  scale_y_discrete(limits=c("2","3","4","5","10","15","20","25","40",
                            "50","75","100","125","150","175"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")+
  labs(x="Point density",y="Point density") +
  guides(color="none")

b

## Save plot
ggsave("./Figures/Tukey/voxel_25.png",height=6,width=13,dpi=600)

## Filter out metrics derived from the CHM
tuk.chm <- tuks.df %>% 
  filter(Buffer==25&metric=="rugosity"&!is.na(p_cat))

## Plot Tukey's for the CHM metrics
c<-ggplot(tuk.chm,aes(x=as.factor(dens1),y=as.factor(dens2),
                      fill=p_cat,color="white"))+
  geom_tile(na.rm=TRUE)+
  facet_grid(rows=vars(Site),cols = vars(metric))+
  scale_fill_manual(values=c("#0d0887","#aa2395","#fdc527"),na.value="white",
                    name = "p-value",
                    labels=c(expression(""<=0.001),expression(""<=0.01),
                             expression(""<=0.05)))+
  scale_color_manual(values=c("white"),
                     name = "",
                     labels = "")+
  scale_x_discrete(limits=c("1","2","3","4","5","10","15","20","25","40",
                            "50","75","100","125","150"))+
  scale_y_discrete(limits=c("2","3","4","5","10","15","20","25","40",
                            "50","75","100","125","150","175"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")+
  labs(x="Point density",y="Point density") +
  guides(color="none")

c

## Save plot
ggsave("./Figures/Tukey/chm_25.png",height=3,width=3.5,dpi=600)

## Boxplots for significantly different variables
## Read in results for Assateague and KRC
ains <- vroom("./Lidar/Lidar_metrics/Collated2/assateague.csv")
krc <- vroom("./Lidar/Lidar_metrics/Collated2/krc.csv") %>% 
  dplyr::select(-c(X))

## Group the results for both sites and keep only results from the 25 m x 25 m buffers
dat <- bind_rows(ains,krc) %>% filter(buffer==25)

## Put data into long format and filter out only results from decimated lidar
cld.dat <- dat %>% 
  pivot_longer(cols=zmax:mhorshanlad,names_to="metric",values_to="estimate") %>% 
  filter(density<=175) %>% 
  mutate(site=case_when(site=="assa"~"AINS",
                        site=="krc"~"KRC")) 

## Create boxplots for the metrics that were significantly different
ggplot(cld.dat %>% filter(metric=="zmode"),aes(x=as.factor(density),y=estimate))+
  geom_boxplot()+
  facet_wrap(~site,scales="free_y")+
  labs(y="Mode height (m)",x="Point density")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/zmode_box.png",height=4,width=6,dpi=600)

ggplot(cld.dat %>% filter(metric=="zmax"),aes(x=as.factor(density),y=estimate))+
  geom_boxplot()+
  facet_wrap(~site,scales="free_y")+
  labs(y="Maximum height (m)",x="Point density")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/zmax_box.png",height=4,width=6,dpi=600)

ggplot(cld.dat %>% filter(metric=="zvd"),aes(x=as.factor(density),y=estimate))+
  geom_boxplot()+
  facet_wrap(~site,scales="free_y")+
  labs(y="Vertical diversity",x="Point density")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/zvd_box.png",height=4,width=6,dpi=600)

ggplot(cld.dat %>% filter(metric=="mlad"),aes(x=as.factor(density),y=estimate))+
  geom_boxplot()+
  facet_wrap(~site,scales="free_y")+
  labs(y="Mean LAD",x="Point density")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/mlad_box.png",height=4,width=6,dpi=600)

ggplot(cld.dat %>% filter(metric=="cvmlad"),aes(x=as.factor(density),y=estimate))+
  geom_boxplot()+
  facet_wrap(~site,scales="free_y")+
  labs(y="CV of mean LAD",x="Point density")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/cvmlad_box.png",height=4,width=6,dpi=600)

ggplot(cld.dat %>% filter(metric=="mhorcvlad"),aes(x=as.factor(density),y=estimate))+
  geom_boxplot()+
  facet_wrap(~site,scales="free_y")+
  labs(y="MHORCVLAD",x="Point density")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/mhorcvlad_box.png",height=4,width=6,dpi=600)

ggplot(cld.dat %>% filter(metric=="shanlad"),aes(x=as.factor(density),y=estimate))+
  geom_boxplot()+
  facet_wrap(~site,scales="free_y")+
  labs(y="Shannon diversity of LAD",x="Point density")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/shanlad_box.png",height=4,width=6,dpi=600)

ggplot(cld.dat %>% filter(metric=="mvershanlad"),aes(x=as.factor(density),y=estimate))+
  geom_boxplot()+
  facet_wrap(~site,scales="free_y")+
  labs(y="Mean vertical Shannon diversity of LAD",x="Point density")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/mvershanlad_box.png",height=4,width=6,dpi=600)

ggplot(cld.dat %>% filter(metric=="mhorshanlad"),aes(x=as.factor(density),y=estimate))+
  geom_boxplot()+
  facet_wrap(~site,scales="free_y")+
  labs(y="Mean horizontal Shannon diversity of LAD",x="Point density")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/mhorshanlad_box.png",height=4,width=6,dpi=600)

ggplot(cld.dat %>% filter(metric=="rugosity"),aes(x=as.factor(density),y=estimate))+
  geom_boxplot()+
  facet_wrap(~site,scales="free_y")+
  labs(y="Rugosity",x="Point density")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/rugosity_box.png",height=4,width=6,dpi=600)


## Reliability ratio figures

## Read in results for the reliability ratio tests and collate results from both sites
assa.rr <- read.csv("./Results/Reliability_ratio/assateague_rr.csv") %>% 
  mutate(Site = "AINS")

krc.rr <- read.csv("./Results/Reliability_ratio/krc_rr.csv") %>% 
  mutate(Site = "KRC")

rr.all <- bind_rows(assa.rr,krc.rr) %>% 
  dplyr::select(-c(X,fhd.var,mcvlad.var))

## Remove ".var" from variable names
names(rr.all) <- str_replace(names(rr.all),"\\.var","")

## Change data to long format
rr.long <- rr.all %>% pivot_longer(cols = c(zmax:mhorshanlad), names_to="metric",
                                   values_to = "reliability_ratio")

## Select cloud metrics that had variability in reliability ratio
rr.cloud <- rr.long %>% 
  filter(metric%in%cloud&buffer==25) %>% 
  mutate(metric=fct_relevel(metric,vars)) %>% 
  filter(!metric%in%c("zmean","zsd","zcov","cover","cdens","lowveg","gini"))

## Select voxel metrics that had variability in reliability ratio
rr.vox <- rr.long %>% 
  filter(metric%in%voxel&buffer==25&!metric%in%c("zvd","mlad","mhorshanlad")) %>% 
  mutate(metric=fct_relevel(metric,vars))

## Plot reliability ratio results
cld.plot <- rr.cloud %>% ggplot(aes(x=as.factor(density),y=reliability_ratio,colour=Site))+
  geom_point()+
  geom_line(aes(group=Site))+
  geom_hline(aes(yintercept=0.90),linetype="dashed")+
  facet_wrap(~metric,scales="free")+
  scale_x_discrete(limits=c("1","2","3","4","5","10","15","20","25","40",
                            "50","75","100","125","150","175"))+
  labs(x="Point density",y="Reliability ratio") +
  scale_colour_manual(values=c("#0d0887","#fdc527"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/Reliability_Ratio/rr25_cloud.png",height=8,width=10,dpi=600)

vox.plot <- rr.vox %>% ggplot(aes(x=as.factor(density),y=reliability_ratio,colour=Site))+
  geom_point()+
  geom_line(aes(group=Site))+
  geom_hline(aes(yintercept=0.90),linetype="dashed")+
  facet_wrap(~metric,scales="free")+
  scale_x_discrete(limits=c("1","2","3","4","5","10","15","20","25","40",
                            "50","75","100","125","150","175"))+
  labs(x="Point density",y="Reliability ratio") +
  scale_colour_manual(values=c("#0d0887","#fdc527"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave("./Figures/Reliability_Ratio/rr25_voxel.png",height=5,width=7,dpi=600)


## Correlation plots

## Read in results from correlation analysis
cors <- read.csv("./Results/Correlations/correlations.csv") %>% 
  filter(buffer==25) %>% 
  mutate(strength = case_when(abs(r_corr) >= 0.7 ~ "strong",
                              abs(r_corr) >= 0.3 & abs(r_corr) < 0.7 ~ "weak",
                              TRUE ~ "NA"),
         strength=ifelse(strength=="NA",NA,strength),
         var2 = fct_relevel(var2,vars),
         var1 = fct_relevel(var1,vars))

## Create plots for correlations between each variable
for(i in 1:length(vars)){
  ggplot(cors %>% filter(var1==vars[i]&var2!=vars[i]),aes(x=as.factor(density),y=r_corr,color=site,shape=strength))+
    geom_line(aes(group=site))+
    geom_point(size=3)+
    geom_hline(aes(yintercept=0),linetype="dashed")+
    scale_colour_manual(values=c("#0d0887","#fdc527"),name="Site")+
    scale_shape_manual(values=c(16,1),
                       name="Correlation",
                       labels=c("|r| \u2265 0.7",
                                "0.3 \u2264 |r| < 0.7",
                                "|r| < 0.3"))+
    facet_wrap(~var2,scales="fixed")+
    labs(x="Point density",y="Pearson's r")+
    theme_bw()+
    theme(legend.position = "bottom")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste0("./Figures/Correlations/",vars[i],".png"),width=10,height=8)
}

