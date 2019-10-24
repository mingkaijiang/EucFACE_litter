#### clear wk space
rm(list=ls(all=TRUE))

### prepare
source("prepare.R")


### get annual precipitation
rain <- calculate_ros_data()


# calculate wood C pool
wood_pool <- make_wood_pool(ring_area)

wood_production_flux <- make_delta_wood_pool_function(inDF=wood_pool, var.col=3)

wood_production_flux$year <- year(wood_production_flux$Date)
wood_production_flux$Trt[wood_production_flux$Ring%in%c(2,3,6)] <- "amb"
wood_production_flux$Trt[wood_production_flux$Ring%in%c(1,4,5)] <- "ele"

### assign rainfall
delta_wood <- assign_rainfall(wood=wood_production_flux,
                              rainfall=rain)

### look at the difference
wood_diff <- summaryBy(delta+Rainfall~year+Trt, data=delta_wood, FUN=c(mean,sd))


p1 <- ggplot() +
    geom_point(wood_diff, mapping=aes(x=Rainfall.mean, y=delta.mean, col=Trt,
                                      shape=as.factor(year)), size = 4,
               position=position_dodge(width=1))+
    geom_errorbar(data=wood_diff, mapping=aes(x=Rainfall.mean, ymin=delta.mean-delta.sd,
                                              ymax=delta.mean+delta.sd, color=Trt),
                  position=position_dodge(width=1))+
    xlab("Rain (mm)") +
    ylab(expression(Delta*C[stem]*" ( g C " * m^2 * " " * yr^1 * ")"))+
    scale_color_manual(values=c("blue3", "red2"))+
    scale_shape_manual(name="Year", values=c(15, 16, 17, 18, 25, 23, 21),
                       label=c("2012", "2013", "2014", "2015", "2016", "2017", "2018"))+
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





