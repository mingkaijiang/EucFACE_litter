calculate_ros_data <- function() {
    
    downloadHIEv(hiev=searchHIEv("ROS_rain_gauge_2012.dat"))
    downloadHIEv(hiev=searchHIEv("ROS_rain_gauge_2013.dat"))
    downloadHIEv(hiev=searchHIEv("ROS_rain_gauge_2014.dat"))
    downloadHIEv(hiev=searchHIEv("ROS_rain_gauge_2015.dat"))
    downloadHIEv(hiev=searchHIEv("ROS_rain_gauge_2016.dat"))
    downloadHIEv(hiev=searchHIEv("ROS_rain_gauge_2017.dat"))
    downloadHIEv(hiev=searchHIEv("ROS_rain_gauge_2018.dat"))
    downloadHIEv(hiev=searchHIEv("ROS_rain_gauge_2019.dat"))
    
    
    ### read in
    myDF12 <- read.csv("download/ROS_rain_gauge_2012.dat", skip=3)
    myDF13 <- read.csv("download/ROS_rain_gauge_2013.dat", skip=3)
    myDF14 <- read.csv("download/ROS_rain_gauge_2014.dat", skip=3)
    myDF15 <- read.csv("download/ROS_rain_gauge_2015.dat", skip=3)
    myDF16 <- read.csv("download/ROS_rain_gauge_2016.dat", skip=3)
    myDF17 <- read.csv("download/ROS_rain_gauge_2017.dat", skip=3)
    myDF18 <- read.csv("download/ROS_rain_gauge_2018.dat", skip=3)
    myDF19 <- read.csv("download/ROS_rain_gauge_2019.dat", skip=3)
    
    colnames(myDF12) <- colnames(myDF13) <- colnames(myDF14) <- colnames(myDF15) <- colnames(myDF16) <- colnames(myDF17) <- colnames(myDF18) <- colnames(myDF19) <- c("DateTime", "Rainfall", "Other")
    
    test <- do.call("rbind", list(myDF12, myDF13, myDF14, myDF15, myDF16, myDF17, myDF18, myDF19))

    ### add year information
    test$DateTime <- as.character(test$DateTime)
    test$Year <- as.numeric(substr(test$DateTime, 1, 4))
    test$Date <- as.Date(test$DateTime,'%Y-%m-%d %H:%M:%S') 
    
    out <- summaryBy(Rainfall~Date+Year, FUN=sum, data=test, keep.names=T, na.rm=T)
    out <- out[out$Year >= 2012, ]
    
    return(out)
}