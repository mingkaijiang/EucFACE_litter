wood_increment_CO2_drought_interaction(rain) {
    
    
    ### calculate wood C pool
    wood_pool <- make_wood_pool(ring_area)
    
    wood_production_flux <- make_delta_wood_pool_function(inDF=wood_pool, var.col=3)
    
    wood_production_flux$year <- year(wood_production_flux$Date)
    wood_production_flux$Trt[wood_production_flux$Ring%in%c(2,3,6)] <- "amb"
    wood_production_flux$Trt[wood_production_flux$Ring%in%c(1,4,5)] <- "ele"
    
    ### assign rainfall
    delta_wood <- assign_rainfall(wood=wood_production_flux,
                                  rainfall=rain)
    
    ### treatment summary
    wood_trt <- summaryBy(delta+Rainfall~year+Trt, data=delta_wood, FUN=c(mean,sd))
    
    ### calculate CO2 response
    wood_diff <- subset(wood_trt, Trt=="amb")
    tmp <- subset(wood_trt, Trt=="ele")
    wood_diff$ele <- tmp$delta.mean
    wood_diff$CO2_effect <- wood_diff$ele - wood_diff$delta.mean
    
    ################## statistics
    #### assign dry wet factor
    delta_wood$dw <- ifelse(delta_wood$Rainfall > 850, "Wet", "Dry")
    wood_diff$dw <- ifelse(wood_diff$Rainfall.mean > 850, "Wet", "Dry")
    
    ### check on delta Cstem
    modelt1 <- lm(delta~Trt * dw, data=delta_wood)
    summary(modelt1)
    
    ### check on CO2 effect of delta Cstem
    modelt2 <- lm(CO2_effect~ dw, data=wood_diff)
    summary(modelt2)
    
    ### summary table for dry wet comparison
    dwDF <- summaryBy(CO2_effect~dw, data=wood_diff, FUN=c(mean,sd))
    
    
    #### Plot treatment comparison
    p1 <- ggplot() +
        geom_point(wood_trt, mapping=aes(x=Rainfall.mean, y=delta.mean, col=Trt,
                                         shape=as.factor(year)), size = 4,
                   position=position_dodge(width=1))+
        geom_errorbar(data=wood_trt, mapping=aes(x=Rainfall.mean, ymin=delta.mean-delta.sd,
                                                 ymax=delta.mean+delta.sd, color=Trt),
                      position=position_dodge(width=1))+
        xlab("Rain (mm)") +
        ylab(expression(Delta*C[stem]*" ( g C " * m^2 * " " * yr^1 * ")"))+
        scale_color_manual(values=c("blue3", "red2"))+
        scale_shape_manual(name="Year", values=c(15, 16, 17, 18, 25, 23, 21),
                           label=c("2012", "2013", "2014", "2015", "2016", "2017", "2018"))+
        theme(panel.grid.minor=element_blank(),
              axis.title.x = element_text(size=12), 
              axis.text.x = element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=12),
              legend.text=element_text(size=12),
              legend.title=element_text(size=12),
              panel.grid.major=element_blank(),
              legend.position="right",
              legend.text.align=0,
              legend.direction="vertical")+
        xlim(500, 1000)+
        guides(size = "none",
               shape = "none",
               #fill = viridis,
               color = guide_legend(ncol=1, override.aes = 
                                        list(size = 4, shape= 15, 
                                             colour=c("blue3", "red2"))))
    
    
    ### plot CO2 effect
    p2 <- ggplot() +
        geom_point(wood_diff, mapping=aes(x=Rainfall.mean, y=CO2_effect, 
                                          shape=as.factor(year)), size = 4,
                   position=position_dodge(width=1))+
        xlab("Rain (mm)") +
        ylab(expression(CO[2] * " response"))+
        scale_shape_manual(name="Year", values=c(15, 16, 17, 18, 25, 23, 21),
                           label=c("2012", "2013", "2014", "2015", "2016", "2017", "2018"))+
        theme(panel.grid.minor=element_blank(),
              axis.title.x = element_text(size=12), 
              axis.text.x = element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=12),
              legend.text=element_text(size=12),
              legend.title=element_text(size=12),
              panel.grid.major=element_blank(),
              legend.position="right",
              legend.text.align=0,
              legend.direction="vertical")+
        xlim(500, 1000)+
        geom_vline(xintercept=850, lty=2)
    
    
    ### plot CO2 effect on dry wet comparison
    p3 <- ggplot()+
        geom_bar(data=dwDF, stat = "identity", aes(dw, CO2_effect.mean, color=dw, fill=dw), 
                 position="dodge") +
        geom_errorbar(data=dwDF, stat = "identity", aes(dw, ymin=CO2_effect.mean-CO2_effect.sd,
                                                        ymax = CO2_effect.mean + CO2_effect.sd), 
                 position="dodge", width=0.5) +
        geom_point(data=wood_diff, mapping=aes(x=dw, y=CO2_effect, shape=as.factor(year), fill=dw), 
                   size=4, position = position_dodge(0.9), color="black")+
        xlab("") + ylab(expression(CO[2] * " effect"))+
        theme_linedraw() +
        theme(panel.grid.minor=element_blank(),
              axis.title.x = element_text(size=12), 
              axis.text.x = element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=12),
              legend.text=element_text(size=12),
              legend.title=element_text(size=12),
              panel.grid.major=element_blank(),
              legend.position="right",
              legend.text.align=0,
              legend.direction="vertical")+
        scale_fill_manual(name="Rainfall", values = c("Dry" = "orange", "Wet" = "lightblue"),
                           labels=c("Dry", "Wet"))+
        scale_color_manual(name="Rainfall", values = c("Dry" = "orange", "Wet" = "lightblue"),
                           labels=c("Dry", "Wet"))+
        scale_shape_manual(name="Year", values=c(15, 16, 17, 18, 25, 23, 21),
                           label=c("2012", "2013", "2014", "2015", "2016", "2017", "2018"))+
        scale_x_discrete("",  
                         labels=c("Dry",
                                  "Wet"))+
        guides(size = "none",
               shape = "none",
               fill = guide_legend(ncol=1, override.aes = 
                                       list(colour=c("orange", "lightblue"),
                                            shape=15)))
        
    
    ### plot file structure
    wood.increment.plot <-
        ggdraw() +
        draw_plot(p1, x = 0, y = 0.7, width =1.0, height =0.3)+
        draw_plot(p2, x = 0, y = 0.35, width =1.0, height =0.3)+
        draw_plot(p3, x = 0, y = 0.0, width =1.0, height =0.3)
    
    
    
    ### save output
    ggsave(filename = "output/wood_increment_plot.pdf", 
           plot = wood.increment.plot,
           width = 12, 
           height = 18,
           units = "cm",
           dpi = 300)
    
    
}