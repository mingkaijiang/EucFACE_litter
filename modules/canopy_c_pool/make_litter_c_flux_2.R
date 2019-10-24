make_litter_c_flux <- function(c_fraction, frass_basket_area){
    
    # download_leaflitter()
    
    #### read in data
    f17 <- downloadCSV("FACE_P0017_RA_Litter_20170101-20171130-L1.csv")
    f17$Date <- parse_date_time(f17$Date,"d m y")
    f16 <- downloadCSV("FACE_P0017_RA_Litter_20160101-20161212-L1-V3.csv")
    f16$Date <- parse_date_time(f16$Date,"d m y")
    f15 <- downloadCSV("FACE_P0017_RA_Litter_20150101-20151217-L1-V2.csv")
    f15$Date <- parse_date_time(f15$Date,"d m y")
    f14 <- downloadCSV("FACE_P0017_RA_Litter_20140101-20141216-L1-V2.csv")
    f14$Date <- parse_date_time(f14$Date,"d m y")
    f13 <- downloadCSV("FACE_P0017_RA_Litter_20121001-20131231-L1-V2.csv")
    f13$Date <- parse_date_time(f13$Date,"d m y")
    
    #### set up col names
    colnames(f13) <- c("Ring", "Date", "Trap", "Twig", "Bark", "Seed", "Leaf", "Other", "Insect", "Comments", "days.past","Source")
    colnames(f14) <- c("Ring", "Date", "Trap", "Twig", "Bark", "Seed", "Leaf", "Other", "Insect", "Comments", "days.past","Source")
    colnames(f15) <- c("Ring", "Date", "Trap", "Twig", "Bark", "Seed", "Leaf", "Other", "Insect", "Comments", "days.past","Source")
    colnames(f16) <- c("Ring", "Date", "Trap", "Twig", "Bark", "Seed", "Leaf", "Other", "Insect", "Comments", "days.past","Source")
    
    # extract the flower buds from 2017
    f17buds <- f17[,c("Ring","Date","Trap","Flower Buds")]
    names(f17buds)[4] <- "Flower_Buds"
    f17 <- f17[,-10]
    
    #### Merge the files
    litter_raw <- rbind(f13, f14, f15, f16, f17) 
    litter_raw <- merge(litter_raw,f17buds,by=c("Ring","Date","Trap"),all.x=T)
    
    # glitch fix
    litter_raw$Ring <- as.character(litter_raw$Ring)
    litter_raw$Trap <- as.character(litter_raw$Trap)
    #litter_raw$Ring[is.na(litter_raw$Ring)] <- litter_raw$RING[is.na(litter_raw$Ring)]
    #litter_raw$TRAP[is.na(litter_raw$Ring)] <- litter_raw$RING[is.na(litter_raw$Ring)]
    litter_raw$Twig <- as.numeric(litter_raw$Twig)
    litter_raw$Bark <- as.numeric(litter_raw$Bark)
    litter_raw$Seed <- as.numeric(litter_raw$Seed)
    litter_raw$Leaf <- as.numeric(litter_raw$Leaf)
    litter_raw$Other <- as.numeric(litter_raw$Other)
    litter_raw$Insect <- as.numeric(litter_raw$Insect)
    
    # remove three data points where big branches fall into litter bascket
    line.num <- which.max(litter_raw$Twig)
    litter_raw <- litter_raw[-line.num,]
    line.num <- which.max(litter_raw$Twig)
    litter_raw <- litter_raw[-line.num,]
    line.num <- which.max(litter_raw$Twig)
    litter_raw <- litter_raw[-line.num,]
    
    # Conversion factor from g basket-1 to mg m-2
    conv <- c_fraction * 1000 / frass_basket_area
    
    litter <- dplyr::mutate(litter_raw, 
                            Date = as.Date(litter_raw$Date, format = "%d/%m/%Y"),
                            Start_date = Date - days.past,
                            End_date = Date,
                            Twig = as.numeric(Twig) * conv / days.past,
                            Bark = as.numeric(Bark) * conv / days.past,
                            Seed = as.numeric(Seed) * conv / days.past,
                            Leaf = as.numeric(Leaf) * conv / days.past)
    
    # Averages by Ring
    litter_a <- summaryBy(Twig + Bark + Seed + Leaf ~ Date + Ring, FUN=mean, na.rm=TRUE,
                          data=litter, id = ~Start_date + End_date, keep.names=TRUE)
    
    litter_a <- as.data.frame(dplyr::rename(litter_a, 
                                            twig_flux = Twig,
                                            bark_flux = Bark,
                                            seed_flux = Seed,
                                            leaf_flux = Leaf))
    
    litter_a$Days <- as.numeric(with(litter_a, End_date - Start_date)) + 1
    
    return(litter_a)
}
