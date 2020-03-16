make_litterfall_2019_plot <- function() {
    myDF <- read.csv("temp_files/EucFACE_litter_2019_data.csv")
    
    # date
    myDF$Date <- parse_date_time(myDF$Date,"d m y")
    
    # summaryBy
    myDF$Trt <- "aCO2"
    myDF$Trt[myDF$Ring%in%c(1,4,5)] <- "eCO2"
    
    plotDF <- summaryBy(Total~Date+Trt, FUN=c(mean,sd), data=myDF, keep.names=T, na.rm=T)
    
    p1 <- ggplot(plotDF,
                       aes(Date, Total.mean, group=Trt)) +  
        geom_bar(stat = "identity", aes(fill=Trt), position="dodge2") +
        geom_errorbar(data=plotDF, mapping=aes(x=Date, ymin=Total.mean-Total.sd, 
                                               ymax=Total.mean+Total.sd, group=Trt), 
                      stat = "identity", position="dodge2") + 
        xlab("") + ylab(expression("Litter flux (g C " * m^-2 * " " * yr^-1 * ")")) +
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
               alpha=guide_legend(ncol=1))+
        scale_fill_manual(values=c("blue3", "red2"))
    
    #plot(p1)
    
    
    ### save output
    ggsave(filename = "output/litterfall_2019_plot.pdf", 
           plot = p1,
           width = 16, 
           height = 12,
           units = "cm",
           dpi = 300)
    
}