##### Main program to look at EucFACE litter data over 2012-2018

################################ set up ########################################
#### system time
Sys.setenv(TZ="Australia/Sydney")

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
                           FUN=sum,data=litter_all,keep.names=T)
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

#### investigate change in lai
source("make_dLAI_litter.R")
source("canopy_C_pool/make_smooth_lai_variable.R")
source("canopy_C_pool/download_canopy_c_related_data.R")
source("canopy_C_pool/smoothplot.R")
source("canopy_C_pool/fitgam.R")
source("canopy_C_pool/splitbydate.R")
source("canopy_C_pool/make_sla_variable.R")
sla <- make_sla_variable()
dlai <- make_dLAI_litter(litter=litter_all,sla_variable=sla) 
with(subset(dlai,Ring==1),plot(Date,dLEAF,ylim=c(-50,100)))
with(subset(dlai,Ring==1),points(Date,leaflit,col="red"))
with(subset(dlai,Ring==1),points(Date,dLEAF+leaflit,col="green"))
abline(h=0)

dlai$Year <- year(dlai$Date)
dlai_yrly <- summaryBy(dLEAF + leaflit~Ring+Year, FUN=sum, dat=dlai,na.rm=T)
dlai_yrly$CO2Treat <- "Amb"
dlai_yrly$CO2Treat[dlai_yrly$Ring %in% c(1,4,5)] <- "Elev"
with(subset(dlai_yrly,Year > 2012 & Year < 2019),
     plot(Year,dLEAF.sum+leaflit.sum,col=as.factor(CO2Treat),ylim=c(0,300)))

lit_a <- litter_all
lit_a$Leafprod <- with(lit_a,leaf_flux*Days)
lit_a$Year <- year(lit_a$Date)
annual <- summaryBy(Leafprod~Year+Ring,FUN=sum,dat=lit_a)

head(annual)
with(annual,plot(Year,Leafprod.sum/1000,col=Ring))

with(subset(litter_all,year(Date)==2016),plot(yday(Date),leaf_flux,ylim=c(0,1500),xlim=c(0,365)))
with(subset(litter_all,year(Date)==2014),points(yday(Date),leaf_flux,col="blue"))
with(subset(litter_all,year(Date)==2017),points(yday(Date),leaf_flux,col="red",pch=19))
with(subset(litter_all,year(Date)==2018),points(yday(Date),leaf_flux,col="yellow",pch=19))

with(subset(litter_all,year(Date)==2016),plot(yday(Date),twig_flux,ylim=c(0,800),xlim=c(0,365)))
with(subset(litter_all,year(Date)==2014),points(yday(Date),twig_flux,col="blue"))
with(subset(litter_all,year(Date)==2017),points(yday(Date),twig_flux,col="red"))

with(subset(litter_all,year(Date)==2016),plot(yday(Date),bark_flux,ylim=c(0,800),xlim=c(0,365)))
with(subset(litter_all,year(Date)==2014),points(yday(Date),bark_flux,col="blue"))
with(subset(litter_all,year(Date)==2017),points(yday(Date),bark_flux,col="red"))


#### do statistics to 2017
dlai$Year <- year(dlai$Date)
dlai_yrly <- summaryBy(dLEAF + leaflit~Ring+Year, FUN=sum, dat=dlai,na.rm=T)
dlai_yrly$CO2Treat <- "Amb"
dlai_yrly$CO2Treat[dlai_yrly$Ring %in% c(1,4,5)] <- "Elev"
with(subset(dlai_yrly,Year > 2012 & Year < 2019),
     plot(Year,dLEAF.sum+leaflit.sum,col=as.factor(CO2Treat),ylim=c(0,300)))

### annual 2017
t2017 <- subset(dlai_yrly, Year == 2017)
t2017$tot <- t2017$dLEAF.sum + t2017$leaflit.sum

mod1 <- lme(tot~CO2Treat, random=~1|Ring,
            data=t2017, method="REML")
anova.lme(mod1, 
          type="sequential", 
          adjustSigma = FALSE)

### daily 2017
dlai2017 <- subset(dlai, Year == 2017)
dlai2017$tot <- with(dlai2017, dLEAF + leaflit)

dlai2017$CO2Treat <- "Amb"
dlai2017$CO2Treat[dlai2017$Ring %in% c(1,4,5)] <- "Elev"

mod2 <- lme(tot~CO2Treat+Date+CO2Treat*Date, random=~1|Ring,
            data=dlai2017, method="ML")

anova.lme(mod2, 
          type="sequential", 
          adjustSigma = FALSE)

### add cumulative leaf production
dlai2017$Date0 <- as.numeric(dlai2017$Date - min(dlai2017$Date))
for (i in 1:6) {
   dlai2017[dlai2017$Ring==i, "cumsum"] <- cumsum(dlai2017[dlai2017$Ring==i, "tot"])
}

mod3 <- lme(cumsum~CO2Treat+Date0+CO2Treat*Date0, random=~1|Ring,
            data=dlai2017, method="ML")

anova.lme(mod2, 
          type="sequential", 
          adjustSigma = FALSE)



#### Look at eCO2/aCO2 over the entire time series
trtDF1 <- summaryBy(Leafprod~Date+CO2Treat+Start_date+End_date+Days+Year, FUN=c(mean, sd), data=lit_a)
trtDF2 <- data.frame(unique(trtDF1$Date), unique(trtDF1$Start_date), unique(trtDF1$End_date))
colnames(trtDF2) <- c("Date", "Start_date", "End_date")

for (i in trtDF2$Date) {
    trtDF2[trtDF2$Date == i, "aCO2"] <- trtDF1$Leafprod.mean[trtDF1$Date == i & trtDF1$CO2Treat == "Amb"]
    trtDF2[trtDF2$Date == i, "eCO2"] <- trtDF1$Leafprod.mean[trtDF1$Date == i & trtDF1$CO2Treat == "Elev"]
    
}

trtDF2$Ratio <- with(trtDF2, eCO2/aCO2)

### test the slope
mod.lm <- lm(Ratio ~ Date, data=trtDF2)

summary(mod.lm)
a <- coefficients(mod.lm)[[2]]
b <- coefficients(mod.lm)[[1]]
rsq <- summary(mod.lm)$adj.r.squared

p1 <- ggplot(trtDF2) +
    geom_point(aes(Date, Ratio)) +
    geom_smooth(aes(Date, Ratio), method = lm)+
    labs(x="Year", y="eCO2/aCO2")+
    theme_linedraw()+
    theme(panel.grid.minor=element_blank(),
          axis.title.x = element_text(size=14), 
          axis.text.x = element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14),
          panel.grid.major=element_blank(),
          legend.position="none")+
    geom_hline(yintercept = 1.0)+
    annotate(geom="text", x=as.Date("2013-04-01"), y=1.6, 
             label=paste0("y = ", round(a,4), "x + ", round(b,4)), 
             color="black") +
    annotate(geom="text", x=as.Date("2013-04-01"), y=1.55, 
             label=paste0("r2 = ", round(rsq, 2), ", p < 0.0001"))

pdf("output/leaf_litter_ratio_over_time_all_time_points.pdf")
plot(p1)
dev.off()

### look at total new leaf production (i.e. litterfall + dLAI) at annual timestep
dlai_yrly$tot_prod <- with(dlai_yrly, (dLEAF.sum + leaflit.sum))
dlai_trt <- summaryBy(dLEAF.sum+leaflit.sum+tot_prod~Year+CO2Treat, data=dlai_yrly, FUN=c(mean, sd))

ann_prod <- data.frame(c(2012:2018), NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
colnames(ann_prod) <- c("Year", "dLEAF_aCO2_mean", "dLEAF_eCO2_mean",
                        "dLEAF_aCO2_sd", "dLEAF_eCO2_sd",
                        "lit_aCO2_mean", "lit_eCO2_mean",
                        "lit_aCO2_sd", "lit_eCO2_sd",
                        "tot_aCO2_mean", "tot_eCO2_mean",
                        "tot_aCO2_sd", "tot_eCO2_sd")

for (i in ann_prod$Year) {
    ann_prod[ann_prod$Year == i, "tot_aCO2_mean"] <- dlai_trt$tot_prod.mean[dlai_trt$Year == i & dlai_trt$CO2Treat == "Amb"]
    ann_prod[ann_prod$Year == i, "tot_eCO2_mean"] <- dlai_trt$tot_prod.mean[dlai_trt$Year == i & dlai_trt$CO2Treat == "Elev"]
    
    ann_prod[ann_prod$Year == i, "tot_aCO2_sd"] <- dlai_trt$tot_prod.sd[dlai_trt$Year == i & dlai_trt$CO2Treat == "Amb"]
    ann_prod[ann_prod$Year == i, "tot_eCO2_sd"] <- dlai_trt$tot_prod.sd[dlai_trt$Year == i & dlai_trt$CO2Treat == "Elev"]
    
    ann_prod[ann_prod$Year == i, "dLEAF_aCO2_mean"] <- dlai_trt$dLEAF.sum.mean[dlai_trt$Year == i & dlai_trt$CO2Treat == "Amb"]
    ann_prod[ann_prod$Year == i, "dLEAF_eCO2_mean"] <- dlai_trt$dLEAF.sum.mean[dlai_trt$Year == i & dlai_trt$CO2Treat == "Elev"]
    
    ann_prod[ann_prod$Year == i, "dLEAF_aCO2_sd"] <- dlai_trt$dLEAF.sum.sd[dlai_trt$Year == i & dlai_trt$CO2Treat == "Amb"]
    ann_prod[ann_prod$Year == i, "dLEAF_eCO2_sd"] <- dlai_trt$dLEAF.sum.sd[dlai_trt$Year == i & dlai_trt$CO2Treat == "Elev"]
    
    ann_prod[ann_prod$Year == i, "lit_aCO2_mean"] <- dlai_trt$leaflit.sum.mean[dlai_trt$Year == i & dlai_trt$CO2Treat == "Amb"]
    ann_prod[ann_prod$Year == i, "lit_eCO2_mean"] <- dlai_trt$leaflit.sum.mean[dlai_trt$Year == i & dlai_trt$CO2Treat == "Elev"]
    
    ann_prod[ann_prod$Year == i, "lit_aCO2_sd"] <- dlai_trt$leaflit.sum.sd[dlai_trt$Year == i & dlai_trt$CO2Treat == "Amb"]
    ann_prod[ann_prod$Year == i, "lit_eCO2_sd"] <- dlai_trt$leaflit.sum.sd[dlai_trt$Year == i & dlai_trt$CO2Treat == "Elev"]
    
}

ann_prod$tot_ratio <- with(ann_prod, tot_eCO2_mean/tot_aCO2_mean)
ann_prod$dLEAF_ratio <- with(ann_prod, dLEAF_eCO2_mean/dLEAF_aCO2_mean)
ann_prod$lit_ratio <- with(ann_prod, lit_eCO2_mean/lit_aCO2_mean)

### test the slope
mod2.lm <- lm(tot_ratio ~ Year, data=ann_prod)

summary(mod2.lm)
a <- coefficients(mod2.lm)[[2]]
b <- coefficients(mod2.lm)[[1]]
rsq <- summary(mod2.lm)$adj.r.squared

p2 <- ggplot(ann_prod) +
    geom_point(aes(Year, tot_ratio)) +
    geom_smooth(aes(Year, tot_ratio), method = lm)+
    labs(x="Year", y="eCO2/aCO2")+
    theme_linedraw()+
    theme(panel.grid.minor=element_blank(),
          axis.title.x = element_text(size=14), 
          axis.text.x = element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14),
          panel.grid.major=element_blank(),
          legend.position="none")+
    geom_hline(yintercept = 1.0)+
    annotate(geom="text", x=2013, y=1.6, 
             label=paste0("y = ", round(a,4), "x + ", round(b,4)), 
             color="black") +
    annotate(geom="text", x=2013, y=1.55, 
             label=paste0("r2 = ", round(rsq, 2), ", p = 0.3"))

plot(p2)

pdf("output/annual_leaf_production_ratio_over_time.pdf")
plot(p2)
dev.off()

### test the slope
mod3.lm <- lm(dLEAF_ratio ~ Year, data=ann_prod)

summary(mod3.lm)
a <- coefficients(mod3.lm)[[2]]
b <- coefficients(mod3.lm)[[1]]
rsq <- summary(mod3.lm)$adj.r.squared

p3 <- ggplot(ann_prod) +
    geom_point(aes(Year, dLEAF_ratio)) +
    geom_smooth(aes(Year, dLEAF_ratio), method = lm)+
    labs(x="Year", y="eCO2/aCO2")+
    theme_linedraw()+
    theme(panel.grid.minor=element_blank(),
          axis.title.x = element_text(size=14), 
          axis.text.x = element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14),
          panel.grid.major=element_blank(),
          legend.position="none")+
    geom_hline(yintercept = 1.0)+
    annotate(geom="text", x=2013, y=1.6, 
             label=paste0("y = ", round(a,4), "x + ", round(b,4)), 
             color="black") +
    annotate(geom="text", x=2013, y=1.55, 
             label=paste0("r2 = ", round(rsq, 2), ", p = 0.3"))

plot(p3)


### test the slope
mod4.lm <- lm(lit_ratio ~ Year, data=ann_prod)

summary(mod4.lm)
a <- coefficients(mod4.lm)[[2]]
b <- coefficients(mod4.lm)[[1]]
rsq <- summary(mod4.lm)$adj.r.squared

p4 <- ggplot(ann_prod) +
    geom_point(aes(Year, lit_ratio)) +
    geom_smooth(aes(Year, lit_ratio), method = lm)+
    labs(x="Year", y="eCO2/aCO2")+
    theme_linedraw()+
    theme(panel.grid.minor=element_blank(),
          axis.title.x = element_text(size=14), 
          axis.text.x = element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14),
          panel.grid.major=element_blank(),
          legend.position="none")+
    geom_hline(yintercept = 1.0)+
    annotate(geom="text", x=2013, y=1.6, 
             label=paste0("y = ", round(a,4), "x + ", round(b,4)), 
             color="black") +
    annotate(geom="text", x=2013, y=1.55, 
             label=paste0("r2 = ", round(rsq, 2), ", p = 0.3"))

plot(p4)




#### Look at eCO2/aCO2 over the entire time series
trtDF3 <- summaryBy(dLEAF+leaflit+leaf_pool~Date+Trt, FUN=c(mean, sd), data=dlai)
trtDF3 <- trtDF3[complete.cases(trtDF3$dLEAF.mea),]
trtDF4 <- data.frame(unique(trtDF3$Date), NA, NA, NA, NA, NA, NA, NA, NA)
colnames(trtDF4) <- c("Date", "dLEAF_aCO2_mean", "dLEAF_eCO2_mean",
                      "dLEAF_aCO2_sd", "dLEAF_eCO2_sd",
                      "lit_aCO2_mean", "lit_eCO2_mean",
                      "lit_aCO2_sd", "lit_eCO2_sd")

for (i in trtDF4$Date) {
    trtDF4[trtDF4$Date == i, "dLEAF_aCO2_mean"] <- trtDF3$dLEAF.mean[trtDF3$Date == i & trtDF3$Trt== "Amb"]
    trtDF4[trtDF4$Date == i, "dLEAF_eCO2_mean"] <- trtDF3$dLEAF.mean[trtDF3$Date == i & trtDF3$Trt == "Elev"]
    
    trtDF4[trtDF4$Date == i, "dLEAF_aCO2_sd"] <- trtDF3$dLEAF.sd[trtDF3$Date == i & trtDF3$Trt== "Amb"]
    trtDF4[trtDF4$Date == i, "dLEAF_eCO2_sd"] <- trtDF3$dLEAF.sd[trtDF3$Date == i & trtDF3$Trt == "Elev"]
    
    trtDF4[trtDF4$Date == i, "lit_aCO2_mean"] <- trtDF3$leaflit.mean[trtDF3$Date == i & trtDF3$Trt== "Amb"]
    trtDF4[trtDF4$Date == i, "lit_eCO2_mean"] <- trtDF3$leaflit.mean[trtDF3$Date == i & trtDF3$Trt == "Elev"]
    
    trtDF4[trtDF4$Date == i, "lit_aCO2_sd"] <- trtDF3$leaflit.sd[trtDF3$Date == i & trtDF3$Trt== "Amb"]
    trtDF4[trtDF4$Date == i, "lit_eCO2_sd"] <- trtDF3$leaflit.sd[trtDF3$Date == i & trtDF3$Trt == "Elev"]
}

trtDF4$dleaf_ratio <- with(trtDF4, dLEAF_eCO2_mean/dLEAF_aCO2_mean)
trtDF4$lit_ratio <- with(trtDF4, lit_eCO2_mean/lit_aCO2_mean)

### test the slope
mod5.lm <- lm(dleaf_ratio ~ Date, data=trtDF4)

summary(mod5.lm)
a5 <- coefficients(mod5.lm)[[2]]
b5 <- coefficients(mod5.lm)[[1]]
rsq5 <- summary(mod5.lm)$adj.r.squared

mod6.lm <- lm(lit_ratio ~ Date, data=trtDF4)

summary(mod6.lm)
a6 <- coefficients(mod6.lm)[[2]]
b6 <- coefficients(mod6.lm)[[1]]
rsq6 <- summary(mod6.lm)$adj.r.squared

p5 <- ggplot(trtDF4) +
    geom_point(aes(Date, dleaf_ratio)) +
    geom_smooth(aes(Date, dleaf_ratio), method = lm)+
    labs(x="Year", y="change in LAI eCO2/aCO2")+
    theme_linedraw()+
    theme(panel.grid.minor=element_blank(),
          axis.title.x = element_text(size=14), 
          axis.text.x = element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14),
          panel.grid.major=element_blank(),
          legend.position="none")+
    geom_hline(yintercept = 1.0)+
    annotate(geom="text", x=as.Date("2013-06-01"), y=20, 
             label=paste0("y = ", round(a5,4), "x + ", round(b5,4)), 
             color="black") +
    annotate(geom="text", x=as.Date("2013-06-01"), y=18, 
             label=paste0("r2 = ", round(rsq5, 2), ", p > 0.1"))


p6 <- ggplot(trtDF4) +
    geom_point(aes(Date, lit_ratio)) +
    geom_smooth(aes(Date, lit_ratio), method = lm)+
    labs(x="Year", y="leaflitter eCO2/aCO2")+
    theme_linedraw()+
    theme(panel.grid.minor=element_blank(),
          axis.title.x = element_text(size=14), 
          axis.text.x = element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14),
          panel.grid.major=element_blank(),
          legend.position="none")+
    geom_hline(yintercept = 1.0)+
    annotate(geom="text", x=as.Date("2013-06-01"), y=1.6, 
             label=paste0("y = ", round(a6,4), "x + ", round(b6,4)), 
             color="black") +
    annotate(geom="text", x=as.Date("2013-06-01"), y=1.55, 
             label=paste0("r2 = ", round(rsq6, 2), ", p < 0.0001"))


pdf("output/change_lai_and_litter_ratio_over_time_all_time_points.pdf")
plot_grid(p5, p6,
          labels="", ncol=1, align="v", axis = "l")

dev.off()


### repeat remko's GCB Figure 5
### i.e. total annual litter production / mean LAI
### and leaf lifespan

## all
lai1 <- summaryBy(LAI~Ring, data=dlai, FUN=mean, na.rm=T)

# Litter
lit <- subset(dlai, !is.na(Days))
lit$dLAIlitter <- with(lit, leaflit / c_fraction * (10^-4 * SLA))


lagg <- summaryBy(leaflit + dLEAF + dLAI + dLAIlitter + Days ~ Ring, FUN=sum, keep.names=TRUE,
                  data=lit)
trapArea <- 0.1979 # g m-2 year-1

lagg$Litter_annual <- (lagg$leaflit / trapArea) * 365.25 / lagg$Days
lagg$LAIlitter_annual <- lagg$dLAIlitter * 365.25 / lagg$Days 

lai1 <- merge(lai1, lagg[,c("Ring","LAIlitter_annual")])

# Leaf lifespan
lai1$LL <- with(lai1, LAI.mean / LAIlitter_annual)


### subset pre-2014
dlai_pre <- subset(dlai, Year < 2014)
lai2 <- summaryBy(LAI~Ring, data=dlai_pre, FUN=mean, na.rm=T)

# Litter
lit <- subset(dlai_pre, !is.na(Days))
lit$dLAIlitter <- with(lit, leaflit / c_fraction * (10^-4 * SLA))

lagg <- summaryBy(leaflit + dLEAF + dLAI + dLAIlitter + Days ~ Ring, FUN=sum, keep.names=TRUE,
                  data=lit)

lagg$Litter_annual <- (lagg$leaflit / trapArea) * 365.25 / lagg$Days
lagg$LAIlitter_annual <- lagg$dLAIlitter * 365.25 / lagg$Days 

lai2 <- merge(lai2, lagg[,c("Ring","LAIlitter_annual")])

# Leaf lifespan
lai2$LL <- with(lai2, LAI.mean / LAIlitter_annual)


### subset 2012-15
dlai_mid <- subset(dlai, Year < 2016)
lai4 <- summaryBy(LAI~Ring, data=dlai_mid, FUN=mean, na.rm=T)

# Litter
lit <- subset(dlai_mid, !is.na(Days))
lit$dLAIlitter <- with(lit, leaflit / c_fraction * (10^-4 * SLA))


lagg <- summaryBy(leaflit + dLEAF + dLAI + dLAIlitter + Days ~ Ring, FUN=sum, keep.names=TRUE,
                  data=lit)

lagg$Litter_annual <- (lagg$leaflit / trapArea) * 365.25 / lagg$Days
lagg$LAIlitter_annual <- lagg$dLAIlitter * 365.25 / lagg$Days 

lai4 <- merge(lai4, lagg[,c("Ring","LAIlitter_annual")])

# Leaf lifespan
lai4$LL <- with(lai4, LAI.mean / LAIlitter_annual)


### subset post 2016
dlai_pos <- subset(dlai, Year >= 2016)
lai3 <- summaryBy(LAI~Ring, data=dlai_pos, FUN=mean, na.rm=T)

# Litter
lit <- subset(dlai_pos, !is.na(Days))
lit$dLAIlitter <- with(lit, leaflit / c_fraction * (10^-4 * SLA))


lagg <- summaryBy(leaflit + dLEAF + dLAI + dLAIlitter + Days ~ Ring, FUN=sum, keep.names=TRUE,
                  data=lit)

lagg$Litter_annual <- (lagg$leaflit / trapArea) * 365.25 / lagg$Days
lagg$LAIlitter_annual <- lagg$dLAIlitter * 365.25 / lagg$Days 

lai3 <- merge(lai3, lagg[,c("Ring","LAIlitter_annual")])

# Leaf lifespan
lai3$LL <- with(lai3, LAI.mean / LAIlitter_annual)

lai4$Trt <- lai1$Trt <- lai2$Trt <- lai3$Trt <- "Amb"
lai4$Trt[lai4$Ring%in%c(1,4,5)] <- lai1$Trt[lai1$Ring%in%c(1,4,5)] <- lai2$Trt[lai2$Ring%in%c(1,4,5)] <- lai3$Trt[lai3$Ring%in%c(1,4,5)] <- "Ele"


p7 <- ggplot(lai1) +
    geom_point(aes(Trt, LL)) +
    ylab("leaf life span (yr)")+
    xlab("2012-2018")

p8 <- ggplot(lai2) +
    geom_point(aes(Trt, LL)) +
    ylab("leaf life span (yr)")+
    xlab("2012-2014")

p9 <- ggplot(lai4) +
    geom_point(aes(Trt, LL)) +
    ylab("leaf life span (yr)")+
    xlab("2012-2015")

p10 <- ggplot(lai3) +
    geom_point(aes(Trt, LL)) +
    ylab("leaf life span (yr)")+
    xlab("2016-2018")

pdf("output/leaf_life_span.pdf")
plot_grid(p8, p9, p10,
          labels="", ncol=1, align="v", axis = "l")

dev.off()

# stats
t.test(LL ~ Trt,data=lai3)
anova.lme(lme(LL ~ Trt, random=~1|Ring, 
    data=lai3, 
    method="REML"))

### only 2016-17
dlai_pos <- subset(dlai, Year >= 2016 & Year < 2018)
lai5 <- summaryBy(LAI~Ring, data=dlai_pos, FUN=mean, na.rm=T)

# Litter
lit <- subset(dlai_pos, !is.na(Days))
lit$dLAIlitter <- with(lit, leaflit / c_fraction * (10^-4 * SLA))


lagg <- summaryBy(leaflit + dLEAF + dLAI + dLAIlitter + Days ~ Ring, FUN=sum, keep.names=TRUE,
                  data=lit)

lagg$Litter_annual <- (lagg$leaflit / trapArea) * 365.25 / lagg$Days
lagg$LAIlitter_annual <- lagg$dLAIlitter * 365.25 / lagg$Days 

lai5 <- merge(lai5, lagg[,c("Ring","LAIlitter_annual")])

# Leaf lifespan
lai5$LL <- with(lai5, LAI.mean / LAIlitter_annual)

lai5$Trt <- "Amb"
lai5$Trt[lai5$Ring%in%c(1,4,5)] <- "Ele"


p11 <- ggplot(lai5) +
    geom_point(aes(Trt, LL)) +
    ylab("leaf life span (yr)")+
    xlab("2016-2017")

plot(p11)


### check the change in leaf life span over two periods
### compare 2016-18 against 2012-14
diffDF <- lai4
diffDF$LL_after <- lai5$LL
diffDF$diff <- with(diffDF, LL_after-LL)

# stats
t.test(diff ~ Trt,data=diffDF)
anova.lme(lme(diff ~ Trt, random=~1|Ring, 
              data=diffDF, 
              method="REML"))

p12 <- ggplot(diffDF) +
    geom_point(aes(Trt, diff)) +
    ylab(expression(paste(Delta, " leaf life span (yr)")))+
    xlab("2016 to 2018 over 2012 to 2015")

plot(p12)
