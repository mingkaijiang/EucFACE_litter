calculate_bom_data <- function() {
    ### read in
    myDF <- read.csv("temp_files/BOM_rainfall_RAAF_Richmond.csv")
    
    out <- summaryBy(Rainfall~Year, FUN=sum, data=myDF, keep.names=T)
    
    return(out)
}