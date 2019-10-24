#### clear wk space
rm(list=ls(all=TRUE))

#### read stuff
library(HIEv)
setToken(tokenfile="~/Documents/Research/Projects/EucFACE_C_Balance/R_repo/tokenfile.txt") 

if(!require(pacman))install.packages("pacman")
pacman::p_load(dplyr, doBy, mgcv, stringr, lubridate, reshape, ggplot2, akima, imputeTS,lme4,
               cowplot)

#### Create output folder
if(!dir.exists("output")) {
    dir.create("output", showWarnings = FALSE)
}

################################ read in data ########################################
#### prepare necessary dataframes
source("make_litter_c_flux.R")
c_fraction <- 0.5
frass_basket_area <- 0.19
lit <- make_litter_c_flux(c_fraction,frass_basket_area)

# read 2017-`18 litter data`
# remove flower buds
l16 <- downloadCSV("FACE_P0017_RA_Litter_20161213_20170119_1.csv")
l17 <- downloadCSV("FACE_P0017_RA_Litter_2017")
l18 <-  downloadCSV("FACE_P0017_RA_Litter_2018")
names(l16) <- c("Ring","Date","Trap","Twig","Bark","Seed","Leaf",
                "Other","Insect","Comments","Source")
l16$Date <- parse_date_time(l16$Date,"d m y")
l17 <- l17[,1:11]
names(l17) <- c("Ring","Date","Trap","Twig","Bark","Seed","Leaf",
                "Other","Insect","Comments","Source")
l17$Date <- parse_date_time(l17$Date,"d m y")
l18 <- l18[,-10]
names(l18) <- c("Ring","Date","Trap","Twig","Bark","Seed","Leaf",
                "Other","Insect","Comments","Source")
l18$Date <- parse_date_time(l18$Date,"d m y")
lit2 <- rbind(l16,l17,l18)
lit2 <- subset(lit2,!is.na(Date))

lit2$Start_date <- parse_date_time(substr(lit2$Source,22,29),"y m d")
lit2$End_date <- parse_date_time(substr(lit2$Source,31,38),"y m d")
lit2$days.past <- with(lit2,as.numeric(End_date - Start_date))

### Conversion factor from g basket-1 to mg m-2
conv <- c_fraction * 1000 / frass_basket_area


### look at 2017 to 2018
litter1718 <- dplyr::mutate(lit2, 
                            Date = as.Date(lit2$Date, format = "%d/%m/%Y"),
                            Twig = as.numeric(Twig) * conv / days.past,
                            Bark = as.numeric(Bark) * conv / days.past,
                            Seed = as.numeric(Seed) * conv / days.past,
                            Leaf = as.numeric(Leaf) * conv / days.past)

# Averages by Ring
litter1718_a <- summaryBy(Twig + Bark + Seed + Leaf ~ Date + Ring, FUN=mean, na.rm=TRUE,
                          data=litter1718, id = ~Start_date + End_date, keep.names=TRUE)

litter1718_a <- as.data.frame(dplyr::rename(litter1718_a, 
                                            twig_flux = Twig,
                                            bark_flux = Bark,
                                            seed_flux = Seed,
                                            leaf_flux = Leaf))
litter1718_a$Days <- as.numeric(with(litter1718_a, End_date - Start_date)) + 1

#### combined all time series data together
litter_all <- rbind(lit,litter1718_a)
litter_all$Year <- year(litter_all$Date)
litter_all$CO2Treat <- "Amb"
litter_all$CO2Treat[litter_all$Ring %in% c(1,4,5)] <- "Elev"
litter_all$leaf_tot <- with(litter_all,leaf_flux*Days/1000)
litter_all$twig_tot <- with(litter_all,twig_flux*Days/1000)
litter_all$bark_tot <- with(litter_all,bark_flux*Days/1000)
litter_annual <- summaryBy(leaf_tot + twig_tot + bark_tot ~ Year + Ring + CO2Treat,
                           FUN=sum,data=litter_all,keep.names=T, na.rm=T)
with(subset(litter_annual,Year > 2012 & Year < 2019),
     plot(Year,leaf_tot,ylim=c(0,300),col=as.factor(CO2Treat),ylab="Leaf Litter (gC m-2 yr-1)"))
with(subset(litter_annual,Year > 2012 & Year < 2019),
     plot(Year,twig_tot,ylim=c(0,100),col=as.factor(CO2Treat),ylab="Twig Litter (gC m-2 yr-1)"))
with(subset(litter_annual,Year > 2012 & Year < 2019),
     plot(Year,bark_tot,ylim=c(0,100),col=as.factor(CO2Treat),ylab="Bark Litter (gC m-2 yr-1)"))
legend("topright",col=c("black","red"),legend=c("Amb","Elev"),pch=1)
with(subset(litter_annual,Year > 2012 & Year < 2019),
     plot(Year,leaf_tot+twig_tot+bark_tot,ylim=c(0,500),col=as.factor(CO2Treat),
          ylab="Total Litter (gC m-2 yr-1)"))

### stem increment
# ring diameter m
ring_diameter <- 25

# ring ground area m2
ring_area <- pi * (ring_diameter/2)^2

# calculate wood C pool
source("make_wood_pool.R")
wood_pool <- make_wood_pool(ring_area)

source("make_delta_wood_pool_function.R")
wood_production_flux <- make_delta_wood_pool_function(inDF=wood_pool, var.col=3)

wood_production_flux$year <- year(wood_production_flux$Date)
wood_production_flux$Trt[wood_production_flux$Ring%in%c(2,3,6)] <- "amb"
wood_production_flux$Trt[wood_production_flux$Ring%in%c(1,4,5)] <- "ele"

### plot
with(subset(wood_production_flux,year > 2012 & year < 2019),
     plot(year,delta,col=as.factor(Trt),ylab="Wood production (gC m-2 yr-1)"))
legend("topright",col=c("black","red"),legend=c("Amb","Elev"),pch=1)


### get annual precipitation
#rainDF <- read.csv("temp_files/rainfall_data_monthly.csv")
#rainDF$year <- year(rainDF$Month)
#rain <- summaryBy(Rain_mm_Tot~year, FUN=sum, data=rainDF, na.rm=T)
rain <- calculate_bom_data()

for (i in 2013:2018) {
    wood_production_flux$rain[wood_production_flux$year==i] <- rain$Rainfall[rain$Year==i]
}

with(wood_production_flux, plot(delta~rain,col=as.factor(Trt)))


### look at the difference
wood_diff <- summaryBy(delta+rain~year+Trt, data=wood_production_flux, FUN=c(mean,sd))

with(wood_diff, plot(delta.mean~rain.mean,col=as.factor(Trt)))

p1 <- ggplot() +
    geom_point(wood_diff, mapping=aes(x=rain.mean, y=delta.mean, col=Trt,
                                      shape=as.factor(year)), size = 4,
               position=position_dodge(width=5))+
    geom_errorbar(data=wood_diff, mapping=aes(x=rain.mean, ymin=delta.mean-delta.sd,
                                              ymax=delta.mean+delta.sd, color=Trt),
                  position=position_dodge(width=5))+
    xlab("Rain (mm)") +
    ylab(expression(Delta*C[stem]*" ( g C " * m^2 * " " * yr^1 * ")"))+
    scale_color_manual(values=c("blue3", "red2"))+
    scale_shape_manual(values=c(15, 16, 17, 18, 25, 23),
                       label=c("2013", "2014", "2015", "2016", "2017", "2018"))+
    theme(panel.grid.minor=element_blank(),
          axis.title.x = element_text(size=14), 
          axis.text.x = element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=12),
          panel.grid.major=element_blank(),
          legend.position="bottom",
          legend.text.align=0)


pdf("output/rainfall_vs_delta_Cstem.pdf")
plot(p1)
dev.off()





