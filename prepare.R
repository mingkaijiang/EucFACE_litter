#### read stuff
library(HIEv)
setToken(tokenfile="~/Documents/Research/Projects/EucFACE_C_Balance/R_repo/tokenfile.txt") 

if(!require(pacman))install.packages("pacman")
pacman::p_load(dplyr, 
               doBy, 
               mgcv, 
               stringr, 
               lubridate, 
               reshape, 
               ggplot2, 
               akima, 
               imputeTS,
               lme4,
               cowplot)

#### Create output folder
if(!dir.exists("output")) {
    dir.create("output", showWarnings = FALSE)
}

if(!dir.exists("download")) {
    dir.create("download", showWarnings = FALSE)
}

setToPath("download")

##### Define constants
### ring diameter m
ring_diameter <- 25

### ring ground area m2
ring_area <- pi * (ring_diameter/2)^2

### c fraction in pools
c_fraction <- 0.5

### litter basket size
frass_basket_area <- 0.19



##### Sourcing all R files in the modules subdirectory
sourcefiles <- dir("modules", pattern="[.]R$", recursive = TRUE, full.names = TRUE)
for(z in sourcefiles)source(z)


