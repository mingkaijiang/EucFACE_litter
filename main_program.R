#### Code scripts to invesigate drought x CO2 interaction at EucFACE

################################################################################################
#### clear wk space
rm(list=ls(all=TRUE))

#### prepare
source("prepare.R")

################################################################################################
#### get annual precipitation
rain <- calculate_ros_data()


################################################################################################
#### investigate CO2 x drought interaction on wood increment
wood_increment_CO2_drought_interaction(rain=rain)

#### investigate CO2 x drought interaction on leaflitter and total litterfall
litterfall_CO2_drought_interaction(rain=rain)