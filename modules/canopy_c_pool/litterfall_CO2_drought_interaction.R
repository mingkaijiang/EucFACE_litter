##### Main program to look at EucFACE litter data over 2012-2018

litterfall_CO2_drought_interaction <- function (rain) {
    
    ### get litter flux
    lit <- make_litter_c_flux(c_fraction,frass_basket_area)
    
    l17 <- downloadCSV("FACE_P0017_RA_Litter_20171116_20180116.csv")
    l18 <-  downloadCSV("FACE_P0017_RA_Litter_2018")
    names(l17) <- c("Ring","Date","Trap","Twig","Bark","Seed","Leaf",
                    "Other","Insect","Flower_Buds","Comments","Source")
    l17$Date <- parse_date_time(l17$Date,"d m y")
    names(l18) <- c("Ring","Date","Trap","Twig","Bark","Seed","Leaf",
                    "Other","Insect","Flower_Buds","Comments","Source")
    l18$Date <- parse_date_time(l18$Date,"d m y")
    lit2 <- rbind(l17,l18)
    lit2 <- subset(lit2,!is.na(Date))
    
    # remove one data point where big branch fell into litter bascket
    line.num <- which.max(lit2$Twig)
    lit2 <- lit2[-line.num,]
    
    lit2$Start_date <- parse_date_time(substr(lit2$Source,22,29),"y m d")
    lit2$End_date <- parse_date_time(substr(lit2$Source,31,38),"y m d")
    lit2$days.past <- with(lit2,as.numeric(End_date - Start_date))
    
    # remove one data point where big branch fell into litter basket ?
    # This is a bit debatable
    lit2$Twig[lit2$Twig > 60] <- NA
    
    # Conversion factor from g basket-1 to mg m-2
    conv <- c_fraction * 1000 / frass_basket_area
    
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
    
    litter_all <- rbind(lit,litter1718_a)
    litter_all$Year <- year(litter_all$Date)
    litter_all$CO2Treat <- "Amb"
    litter_all$CO2Treat[litter_all$Ring %in% c(1,4,5)] <- "Elev"
    
    
    ### calculate sum of flux over a period
    litter_all$leaf_tot <- with(litter_all,leaf_flux*Days/1000)
    litter_all$twig_tot <- with(litter_all,twig_flux*Days/1000)
    litter_all$bark_tot <- with(litter_all,bark_flux*Days/1000)
    litter_all$tot <- with(litter_all,leaf_tot+twig_tot+bark_tot)

    
    ### for output
    litter_ann <- summaryBy(tot+leaf_tot+twig_tot+bark_tot+Days~Year+Ring+CO2Treat, 
                            FUN=sum, data=litter_all, keep.names=T, na.rm=T)
    
    litter_ann$leaf_ann <- with(litter_ann, leaf_tot)#/Days*365)
    litter_ann$twig_ann <- with(litter_ann, twig_tot)#/Days*365)
    litter_ann$bark_ann <- with(litter_ann, bark_tot)#/Days*365)
    litter_ann$tot_ann <- with(litter_ann, leaf_ann+twig_ann+bark_ann)
    
    
    return(litter_ann)
    
    ### for plotting
    #litter_av <- summaryBy(tot+leaf_tot+twig_tot+bark_tot~CO2Treat+Date,data=litter_all,
    #                       FUN=mean,keep.names=T)
    #litter_avA <- subset(litter_av,CO2Treat=="Amb")
    #litter_avE <- subset(litter_av,CO2Treat=="Elev")
    #cf_litter <- merge(litter_avA,litter_avE,by="Date")
    #with(cf_litter,plot(Date,log(tot.y/tot.x)))
    #trend <- lm(log(tot.y/tot.x)~Date,data=cf_litter)
    #abline(trend)
    #abline(h=0,lty=2)
    #with(cf_litter,plot(Date,log(leaf_tot.y/leaf_tot.x)))
    #trend_leaf <- lm(log(leaf_tot.y/leaf_tot.x)~Date,data=cf_litter)
    #abline(trend_leaf)
    #summary(trend_leaf)
    #
    #with(cf_litter,plot(Date,log(bark_tot.y/bark_tot.x)))
    #with(cf_litter,plot(Date,log(twig_tot.y/twig_tot.x)))
    #trend_twig <- lm(log(twig_tot.y/twig_tot.x)~Date,data=cf_litter)
    #abline(trend_twig)
    #summary(trend_twig)
    #
    #
    #litter_annual <- summaryBy(leaf_tot + twig_tot + bark_tot ~ Year + Ring + CO2Treat,
    #                           FUN=sum,data=litter_all,keep.names=T)
    #with(subset(litter_annual,Year > 2012 & Year < 2019),
    #     plot(Year,leaf_tot,ylim=c(0,300),col=as.factor(CO2Treat),ylab="Leaf Litter (gC m-2 yr-1)"))
    #with(subset(litter_annual,Year > 2012 & Year < 2019),
    #     plot(Year,twig_tot,ylim=c(0,100),col=as.factor(CO2Treat),ylab="Twig Litter (gC m-2 yr-1)"))
    #with(subset(litter_annual,Year > 2012 & Year < 2019),
    #     plot(Year,bark_tot,ylim=c(0,100),col=as.factor(CO2Treat),ylab="Bark Litter (gC m-2 yr-1)"))
    #legend("topright",col=c("black","red"),legend=c("Amb","Elev"),pch=1)
    #with(subset(litter_annual,Year > 2012 & Year < 2019),
    #     plot(Year,leaf_tot+twig_tot+bark_tot,ylim=c(0,500),col=as.factor(CO2Treat),
    #          ylab="Total Litter (gC m-2 yr-1)"))
    #annual_av <- summaryBy(leaf_tot+twig_tot+bark_tot~CO2Treat+Year, FUN=mean,
    #                       keep.names=T,data=litter_annual)
    #annual_avA <- subset(annual_av,CO2Treat=="Amb")
    #annual_avE <- subset(annual_av,CO2Treat=="Elev")
    #cf_annual <- merge(annual_avA,annual_avE,by="Year")
    #cf_annual$leaf.resp <- with(cf_annual,leaf_tot.y/leaf_tot.x)
    #cf_annual$tot.resp <- with(cf_annual,(leaf_tot.y+twig_tot.y+bark_tot.y) /
    #                               (leaf_tot.x+twig_tot.x+bark_tot.x))
    
    
    
    
    
    
    #sla <- make_sla_variable()
    #dlai <- make_dLAI_litter(litter=litter_all,sla_variable=sla) 
    #with(subset(dlai,Ring==1),plot(Date,dLEAF,ylim=c(-50,100)))
    #with(subset(dlai,Ring==1),points(Date,leaflit,col="red"))
    #with(subset(dlai,Ring==1),points(Date,dLEAF+leaflit,col="green"))
    #abline(h=0)
    #
    #dlai$Year <- year(dlai$Date)
    #dlai$CO2Treat <- "Amb"
    #dlai$CO2Treat[dlai$Ring %in% c(1,4,5)] <- "Elev"
    #with(dlai,plot(Date,leaflit,col=as.factor(CO2Treat)))
    #dlai_yrly <- summaryBy(dLEAF + leaflit~Ring+Year+CO2Treat, FUN=sum, dat=dlai,na.rm=T)
    #
    #with(subset(dlai_yrly,Year > 2012 & Year < 2019),
    #     plot(Year,dLEAF.sum+leaflit.sum,col=as.factor(CO2Treat),ylim=c(0,300)))
    #
    #lit_a <- litter_all
    #lit_a$Leafprod <- with(lit_a,leaf_flux*Days)
    #lit_a$Year <- year(lit_a$Date)
    #annual <- summaryBy(Leafprod~Year+Ring,FUN=sum,dat=lit_a)
    #
    #head(annual)
    #with(annual,plot(Year,Leafprod.sum/1000,col=Ring))
    #
    #with(subset(litter_all,year(Date)==2016),plot(yday(Date),leaf_flux,ylim=c(0,1500),xlim=c(0,365)))
    #with(subset(litter_all,year(Date)==2014),points(yday(Date),leaf_flux,col="blue"))
    #with(subset(litter_all,year(Date)==2017),points(yday(Date),leaf_flux,col="red",pch=19))
    #with(subset(litter_all,year(Date)==2018),points(yday(Date),leaf_flux,col="yellow",pch=19))
    #
    #with(subset(litter_all,year(Date)==2016),plot(yday(Date),twig_flux,ylim=c(0,800),xlim=c(0,365)))
    #with(subset(litter_all,year(Date)==2014),points(yday(Date),twig_flux,col="blue"))
    #with(subset(litter_all,year(Date)==2017),points(yday(Date),twig_flux,col="red"))
    #
    #with(subset(litter_all,year(Date)==2016),plot(yday(Date),bark_flux,ylim=c(0,800),xlim=c(0,365)))
    #with(subset(litter_all,year(Date)==2014),points(yday(Date),bark_flux,col="blue"))
    #with(subset(litter_all,year(Date)==2017),points(yday(Date),bark_flux,col="red"))
    #
    #
    #
    #### do statistics to 2017
    #dlai$Year <- year(dlai$Date)
    #dlai_yrly <- summaryBy(dLEAF + leaflit~Ring+Year, FUN=sum, dat=dlai,na.rm=T)
    #dlai_yrly$CO2Treat <- "Amb"
    #dlai_yrly$CO2Treat[dlai_yrly$Ring %in% c(1,4,5)] <- "Elev"
    #with(subset(dlai_yrly,Year > 2012 & Year < 2019),
    #     plot(Year,dLEAF.sum+leaflit.sum,col=as.factor(CO2Treat),ylim=c(0,300)))
    #
    ### annual 2017
    #t2017 <- subset(dlai_yrly, Year == 2017)
    #t2017$tot <- t2017$dLEAF.sum + t2017$leaflit.sum
    #
    #simple <- lm(tot~CO2Treat,data=t2017)
    #confint(simple)
    #
    #mod1 <- lme(tot~CO2Treat, random=~1|Ring,
    #            data=t2017, method="REML")
    #anova.lme(mod1, 
    #          type="sequential", 
    #          adjustSigma = FALSE)
    #
    ### daily 2017
    #dlai2017 <- subset(dlai, Year == 2017)
    #dlai2017$tot <- with(dlai2017, dLEAF + leaflit)
    #
    #dlai2017$CO2Treat <- "Amb"
    #dlai2017$CO2Treat[dlai2017$Ring %in% c(1,4,5)] <- "Elev"
    #
    #mod2 <- lme(tot~CO2Treat+Date+CO2Treat*Date, random=~1|Ring,
    #            data=dlai2017, method="ML")
    #
    #anova.lme(mod2, 
    #          type="sequential", 
    #          adjustSigma = FALSE)
    #
    ### add cumulative leaf production
    #dlai2017$Date0 <- as.numeric(dlai2017$Date - min(dlai2017$Date))
    #for (i in 1:6) {
    #    dlai2017[dlai2017$Ring==i, "cumsum"] <- cumsum(dlai2017[dlai2017$Ring==i, "tot"])
    #}
    #
    #mod3 <- lme(cumsum~CO2Treat+Date0+CO2Treat*Date0, random=~1|Ring,
    #            data=dlai2017, method="ML")
    #
    #anova.lme(mod2, 
    #          type="sequential", 
    #          adjustSigma = FALSE)
    #
    
    
    
    
    
    
    
    
    
    
    
}








