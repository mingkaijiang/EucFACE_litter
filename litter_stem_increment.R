#### clear wk space
rm(list=ls(all=TRUE))

### prepare
source("prepare.R")


### stem increment

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





