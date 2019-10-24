assign_rainfall <- function(wood, rainfall) {
    
    s.list <- unique(wood$Start_date)
    e.list <- unique(wood$End_date)
    
    for (i in 1:length(s.list)) {
        
        ### calculate sum of rainfall
        subDF <- subset(rainfall, Date <= e.list[i] & Date > s.list[i])
        sum.rainfall <- sum(subDF$Rainfall)
        wood$Rainfall[wood$Start_date == s.list[i]] <- sum.rainfall
    }
    
    return(wood)
    
}