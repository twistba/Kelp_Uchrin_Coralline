## --- Kelp & Urchin Coralline Settlement ---
##Relative coralline abundance- based on Hind et al. 2019
## Author: Brenton Twist 
## Date Created: 2022-10-05
#---------------------------------------------------------

# This script contains code for loading data from Hind et al. 2019 paper on coralline abundance in kelp and urchin barrens
# relative abundances in the two habitats are calculated based off this data


######################
#   Load packages   ##
######################

library(tidyverse)
library(magrittr)


##################
#   Load data   ##
##################

# Read in slightly modified Hind et al. 2019 data set
# Found on github https://github.com/martonelab/SubtidalCorallines2018
quads<-read_csv("./coralline_relative-abund/Hind_2019_quad_data.csv")
str(quads)


##################
#   Tidy data   ##
##################

quads %<>%  mutate_at( c("Site_name","Type","Species_name_paper","Species_name_martone"), as.factor)

#Need to add zeros for each time the species wasn't found
quads_wide <- quads %>% dplyr::select(Site,Site_name,Quad,Type,Species_name_martone,Cover) %>%
  group_by(Site,Site_name,Quad,Type,Species_name_martone)%>% #need only one entry of species per site quad
  dplyr::summarise(Cover= mean(Cover))%>%
  spread(Species_name_martone, Cover)
quads_wide %<>% replace(is.na(.), 0) # Replace NA with zero

#back to long format with zeros
quads2<-quads_wide %>% gather( Species_name_martone, Cover, 'Bossiella exarticulata':'Neopolyporolithon sp' )

#####################
#   Analyze data   ##
#####################

#Summaries data
cover_site<-quads2 %>% group_by(Site_name,Type, Species_name_martone)%>%
  dplyr::summarise(Cover= mean(Cover))

cover_type <-cover_site %>% group_by(Type, Species_name_martone)%>%
  dplyr::summarise(Cover= mean(Cover))

cover_type_wide<-as.data.frame(cover_type) %>% spread(Type, Cover)  

cover_type_wide %<>% mutate(Kelp_relative_abun= Kelp/(Kelp+Urchin),Urchin_relative_abun= Urchin/(Kelp+Urchin))%>%
  dplyr::rename(Kelp_cover_total= Kelp, Urchin_cover_total= Urchin) %>%
  arrange(Kelp_relative_abun)

cover_type_wide

write.csv(cover_type_wide,"./coralline_relative-abund/relative_abundance_hind2019.csv")

(cor_rel_abun<-cover_type_wide %>% dplyr::select(-Kelp_cover_total,-Urchin_cover_total))
save(cor_rel_abun,cover_type_wide,file="./coralline_relative-abund/relative_coralline.R")



