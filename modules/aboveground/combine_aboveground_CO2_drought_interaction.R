combine_aboveground_CO2_drought_interaction <- function(wood, litter) {
    
    ### combine the plots
    wood$Year <- year(wood$Date)
    
    for (i in unique(litter$Year)) {
        for (j in 1:6) {
            litter$wood_ann[litter$Year == i & litter$Ring == j] <- wood$delta[wood$Year == i & wood$Ring == j]
        }
    }
    
    ### calculate aboveground total
    litter$abg_tot <- with(litter, tot_ann + wood_ann)
    
    ### exclude 2012 because it is partial data
    myDF <- subset(litter, Year > 2012)

    ### sumamarize by Trt and year for total
    abgDF <- summaryBy(abg_tot+leaf_ann+twig_ann+bark_ann+wood_ann~Year+CO2Treat, FUN=c(mean, sd), data=myDF)
    
    ### change to long format
    plotDF1 <- abgDF[,c("Year", "CO2Treat", "abg_tot.mean", "abg_tot.sd")]
    colnames(plotDF1) <- c("Year", "Trt", "mean", "sd")
    
    subDF1 <- abgDF[,c("Year", "CO2Treat", "leaf_ann.mean", "leaf_ann.sd")]
    subDF2 <- abgDF[,c("Year", "CO2Treat", "twig_ann.mean", "twig_ann.sd")]
    subDF3 <- abgDF[,c("Year", "CO2Treat", "bark_ann.mean", "bark_ann.sd")]
    subDF4 <- abgDF[,c("Year", "CO2Treat", "wood_ann.mean", "wood_ann.sd")]
    
    subDF1$Component <- "Leaf"
    subDF2$Component <- "Twig"
    subDF3$Component <- "Bark"
    subDF4$Component <- "Wood"
    
    colnames(subDF1) <- colnames(subDF2) <- colnames(subDF3) <- colnames(subDF4) <- c(
        "Year", "Trt", "mean", "sd", "Component"
    )
    
    plotDF2 <- do.call("rbind", list(subDF1, subDF2, subDF3, subDF4))
    
    ### assign plot labels
    plotDF2$x.brk <- ifelse(plotDF2$Trt == "Amb", plotDF2$Year - 0.15, plotDF2$Year + 0.15)
    plotDF1$x.brk <- ifelse(plotDF1$Trt == "Amb", plotDF1$Year - 0.15, plotDF1$Year + 0.15)
    
    
    #### Plot treatment comparison
    p1 <- ggplot(plotDF2,
                 aes(x.brk, mean, group=Trt)) +  
        geom_hline(yintercept=0)+
        geom_bar(stat = "identity", aes(fill=Component, alpha=Trt),
                 width=0.2, col="black") +
        geom_errorbar(data=plotDF1, mapping=aes(x=x.brk, ymin=mean-sd, ymax=mean+sd, group=Trt), 
                      width=0.1, size=0.8, color="black", stat = "identity") + 
        geom_point(data=plotDF1, mapping=aes(x=x.brk, y=mean), size=2, shape=21, fill="white")+
        xlab("") + ylab(expression(Delta * C[abg] * " (g C " * m^-2 * " " * yr^-1 * ")")) +
        scale_x_continuous(breaks = c(2013, 2014, 2015, 2016, 2017, 2018))+
        scale_alpha_manual(name="Treatment", values=c(1.0, 0.4))+
        theme_linedraw() +
        theme(panel.grid.minor=element_blank(),
              axis.title.x = element_text(size=16), 
              axis.text.x = element_text(size=20),
              axis.text.y=element_text(size=16),
              axis.title.y=element_text(size=18),
              legend.text=element_text(size=14),
              legend.title=element_text(size=14),
              panel.grid.major=element_blank(),
              legend.position="bottom",
              legend.text.align=0)+
        guides(fill=guide_legend(ncol=2),
               alpha=guide_legend(ncol=1))

    #plot(p1)
    

    ### save output
    ggsave(filename = "output/aboveground_increment_plot.pdf", 
           plot = p1,
           width = 16, 
           height = 12,
           units = "cm",
           dpi = 300)
        
    
}