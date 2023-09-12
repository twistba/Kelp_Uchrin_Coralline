## --- Kelp & Urchin Coralline Settlement ---
## Kelp Settlement
## Author:       Brenton Twist ##
## Date Created: 2020-10-07   
#-------------------------------------------------

# This script contains code for loading, tidying and analysis kelp spore settlement data onto coralline algal species


######################
#   Load packages   ##
######################

library(tidyverse)
library(lubridate)
library(ggthemes)
library(Rmisc)
library(car)
library(cowplot)
library(magrittr)
library(lme4)
library(lmerTest)
library(ggrepel)
library(rstatix)
library(ggpubr)
library(multcomp)
library(emmeans)
library(FSA)
library(sandwich)
library(xtable)


##################
#   Load data   ##
##################
#   The following code loads the kelp spore settlement data

#Read the metadata file with row numbers for each coralline species
spore<-read_csv("./kelp/Kelp_spore_counts_raw.csv")

spore %<>% mutate(UID=paste(Kelp_species,Plate,sep="-"))  %>% 
  mutate_at( c("Kelp_species","Period"), as.factor)
str(spore)


##################
#   Tidy data   ##
##################
##Select for the highest density out of the two count periods

spore_wide<- spore %>% dplyr::select(UID,Plate,Plate_numb,Kelp_species,Morpho_species,Species_updated,
                                     Type,Unique_site,Rep,Period,Density_cm) %>%
  gather(key, value, Density_cm  ) %>%  
  unite(new.col, c(key, Period)) %>%   
  spread(new.col, value) 

str(spore_wide)

spore_wide %<>% mutate(Density_cm_max=pmax(Density_cm_1,Density_cm_2,na.rm = TRUE), 
                       Density_cm_diff = Density_cm_2-Density_cm_1)


spore_wide$max_density_time<-colnames(spore_wide[,c("Density_cm_1","Density_cm_2")])[
  max.col(replace(spore_wide[,c("Density_cm_1","Density_cm_2")], is.na(spore_wide[,c("Density_cm_1","Density_cm_2")]), 0),
          ties.method="first")] # where the max density was recorded - If the measure was zero for both times always selects time 1


spore_wide %<>%  mutate(Density_cm_log = log(Density_cm_max + 1),Density_cm_p1= (Density_cm_max + 1)) # add log of density


## Add relative abundance of coralline species to dataset
load(file="./coralline_relative-abund/relative_coralline.Rdata") # load in relative abundance
cor_rel_abun<-cor_rel_abun %>% dplyr::rename(Species_updated=Species_name_martone)

(relative_abund<-left_join(data.frame(Species_updated=unique(spore_wide$Species_updated)),
                           cor_rel_abun, by="Species_updated"))

relative_abund %<>% mutate(Kelp_relative_abun = round(Kelp_relative_abun, 2),
                           Urchin_relative_abun = round(Urchin_relative_abun, 2))  #round to 2 dp

fleshy<-cor_rel_abun%>% dplyr::filter(Species_updated=="Fleshy red crust")
kelp_fleshy<-round(fleshy$Kelp_relative_abun,2)
urchin_fleshy<-round(fleshy$Urchin_relative_abun,2)

#lets add a relative abundance for two ACA and fleshy reds
str(relative_abund)
unique(relative_abund$Species_updated)

relative_abund %<>% 
  mutate(Kelp_relative_abun=replace(Kelp_relative_abun, 
                                    Species_updated=="Bossiella orbigniana"| 
                                      Species_updated=="Calliarthron tuberculosum", 1.00),
         Urchin_relative_abun=replace(Urchin_relative_abun, 
                                      Species_updated=="Bossiella orbigniana"| 
                                        Species_updated=="Calliarthron tuberculosum", 0.00),
         Kelp_relative_abun=replace(Kelp_relative_abun, 
                                    Species_updated=="Hildenbrandia spp"| 
                                      Species_updated=="Peyssonnelia spp", kelp_fleshy),
         Urchin_relative_abun=replace(Urchin_relative_abun, 
                                      Species_updated=="Hildenbrandia spp"| 
                                        Species_updated=="Peyssonnelia spp", urchin_fleshy)) 


spore_wide<-left_join(spore_wide,relative_abund, by="Species_updated")


##Reorder levels
spore_wide$Species_updated %<>% as.factor 
spore_wide$Species_updated <- factor(spore_wide$Species_updated, 
                                     levels=c("Cover slip",  "Rock", "Hildenbrandia spp","Peyssonnelia spp",  "Dead Crusticorallina" ,
                                              "Lithophyllum Corsp27BC", "Chiharaea americana f. americana", "Crusticorallina muricata",
                                              "Crusticorallina painei","Chiharaea americana f. bodegensis" ,"Lithophyllum Corsp4BC",
                                              "Lithothamnion glaciale", "Bossiella schmittii", "Calliarthron tuberculosum", 
                                              "Bossiella orbigniana" )) 
levels(spore_wide$Species_updated)

#add a count of each genetically ID coralline species for each kelp-  so can filter
species_counts<- spore_wide%>%                  
  filter(!is.na(Species_updated)) %>%    
  group_by(Kelp_species,Species_updated) %>%        
  dplyr::summarise(Species_updated_count = n()) 

spore_wide<-left_join(spore_wide, species_counts, by=c("Kelp_species","Species_updated"))


#subset into kelp species
nereo<- as.data.frame(spore_wide) %>% dplyr::filter(Kelp_species == "Nereocystis" )
macro<- as.data.frame(spore_wide) %>% dplyr::filter(Kelp_species == "Macrocystis" )
cost<- as.data.frame(spore_wide) %>% dplyr::filter(Kelp_species == "Costaria Costata" )


#############################
#   Statistical analysis   ##
#############################


#### --- Costaria ---- ####
#--------------------------

cost3<-cost %>%   filter(!is.na(Density_cm_max)) %>% 
  filter(Species_updated!="Cover slip" ) %>% # Remove cover slip- not intended as control based on reveiwer comments
  filter(Species_updated_count >= 3) # removes species with less than 3 replicates & those without densities

cost_lm <- cost3 %>% lm(Density_cm_max ~ Species_updated , data = .)
# check assumptions - use plots and not test because of small sample size 
# https://stats.stackexchange.com/questions/2492/is-normality-testing-essentially-useless
#shapiro.test(residuals(cost_lm))#P less than 0.05- not normally distributed
#leveneTest(Density_cm_max ~ Species_updated, data = cost3) #p less than 0.05 significant difference between variances across groups
hist(cost_lm$residuals,breaks=20)
ggqqplot(residuals(cost_lm)) #look at qqplot to see normality 
plot(cost_lm, 1) # look at residuals to check the homogeneity of variances.
#doesn't look normally distributed- lets log transform


#Log transformed model
cost_lm <- cost3 %>%   lm(Density_cm_log ~ Species_updated , data = .)
# check assumptions- as above use plots and not test because of small sample size 
#shapiro.test(residuals(cost_lm))#P less than 0.05- not normally distributed- looks better below though
#leveneTest(Density_cm_max ~ Species_updated, data = cost3) #p less than 0.05 significant difference between variances across groups
hist(cost_lm$residuals,breaks=20)
ggqqplot(residuals(cost_lm))
plot(cost_lm, 1) # look at residuals to check the homogeneity of variances.
# Looks good in terms of homogeneity of variances and normality- anova robust to slight deviations
# Unbalanced design- Only a one way anova so don't need to worry about changing to typeII or III with car::Anova(., type = 2) 

cost_aov<-aov(cost_lm) #Type I anova
summary(cost_aov)

#write results
cost_out1<-xtable(cost_aov)
write.csv(cost_out1,file="./kelp/results/cost_aov_results.csv")

#Posthoc Tukeys procedure- We have unbalanced design so we can use ones tailored towards unbalanced designs
# start with using Glht as good for unbalanced designs
#adding heteroskedasticity-consistent covariance matrix estimation (vcov = vcovHC) following Herberich et al. 2010

cost_tukey<-glht(model=cost_lm, linfct = mcp(Species_updated = "Tukey"),  vcov = vcovHC)
summary(cost_tukey)
cld(cost_tukey,alpha = .05, Letters = letters) # returns the letters
# warning message returned
#https://stats.stackexchange.com/questions/328463/tukey-with-interaction-and-extracting-letters-in-multcomp


#lets try Emmeans- which is comparable and similar method for unbalanced designs
#https://stackoverflow.com/questions/62963945/glht-and-emmeans-returning-crazy-compact-letters-for-unbalanced-dataset-in-r
cost_emmeans <- emmeans(cost_lm, "Species_updated",  vcov = vcovHC) # removed C. americana f. american weird results with it in
summary(cost_emmeans)
(cost_out2<-pairs(cost_emmeans)) #default Tukey method
(cost_cld<-cld(cost_emmeans, Letter = letters) )
#same results as GLHT posthoc tukey test

#write results
write.csv(cost_out2,file="./kelp/results/cost_tukey_results.csv")

#add letters to data set
(cost_tukey_letters<-cost_cld %>% mutate(Kelp_species="Costaria Costata") %>% 
    dplyr::rename("letter"=".group")%>%
    dplyr::select(Kelp_species,Species_updated,letter))


(tukey_letters<-cost_tukey_letters)



#### --- Nereo --- ####
#--------------------


nereo3<-nereo %>%   filter(!is.na(Density_cm_max)) %>% 
  filter(Species_updated!="Cover slip" ) %>% # Remove cover slip- not intended as control based on reveiwer comments
  filter(Species_updated_count >= 3) # removes species with less than 3 replicates & those without densities

#linear model
nereo_lm <- nereo3 %>% lm(Density_cm_max ~ Species_updated , data = .)
#check assumptions - as above use plots and not test because of small sample size 
#shapiro.test(residuals(nereo_lm))#P less than 0.05- not normally distributed
#leveneTest(Density_cm_max ~ Species_updated, data = nereo3) #p less than 0.05 significant difference between variances across groups
hist(nereo_lm$residuals,breaks=20)
ggqqplot(residuals(nereo_lm))
plot(nereo_lm, 1) # look at residuals to check the homogeneity of variances.
#Not looking good lets try log transformation

#Use log values
nereo_lm <- nereo3 %>%   lm(Density_cm_log ~ Species_updated , data = .)
#check assumptions - as above use plots and not test because of small sample size
#shapiro.test(residuals(nereo_lm))#P greater than 0.05-  normally distributed
#leveneTest(Density_cm_max ~ Species_updated, data = nereo3) #p less than 0.05 significant difference between variances across groups
hist(nereo_lm$residuals,breaks=20)
ggqqplot(residuals(nereo_lm))
plot(nereo_lm, 1) # look at residuals to check the homogeneity of variances.
# Looks good in terms of homogeneity of variances and normality anova robust to slight deviations


#Anova results
nereo_aov<-aov(nereo_lm)
summary(nereo_aov)

#write results
nereo_out1<-xtable(nereo_aov)
write.csv(nereo_out1,file="./kelp/results/nereo_aov_results.csv")

#Posthoc Tukeys procedure- We have unbalanced design so we can use ones tailored towards unbalanced designs
# start with using Glht as good for unbalanced designs
nereo_tukey<-glht(model=nereo_lm, linfct = mcp(Species_updated = "Tukey"), vcov = vcovHC)
summary(nereo_tukey)
cld(nereo_tukey,alpha = .05, Letters = letters) # returns the letters
# warning message returned

#Try Emmeans- which is comparable and similar method for unbalanced designs
nereo_emmeans <- emmeans(nereo_lm, "Species_updated", vcov = vcovHC)
summary(nereo_emmeans)
(nereo_out2<-pairs(nereo_emmeans))
(nereo_cld<-cld(nereo_emmeans, Letter = letters))
#same results as GLHT posthoc tukey test 
#note slightly different lettering for C. muricata & B. orbigninna between GLHT & Emmeans but is same result

#write results
write.csv(nereo_out2,file="./kelp/results/nereo_tukey_results.csv")


#add letters to data set
(nereo_tukey_letters<-nereo_cld %>% mutate(Kelp_species="Nereocystis") %>% 
    dplyr::rename("letter"=".group")%>%
    dplyr::select(Kelp_species,Species_updated,letter))

(tukey_letters<-rbind(tukey_letters,nereo_tukey_letters))



#### ---Macro--- ####
#--------------------

macro3<-macro %>%   filter(!is.na(Density_cm_max)) %>% 
  filter(Species_updated!="Cover slip" ) %>% # Remove cover slip- not intended as control based on reveiwer comments
  filter(Species_updated_count >= 3) # removes species with less than 3 replicates & those without densities

#linear model
macro_lm <- macro3 %>% lm(Density_cm_max ~ Species_updated , data = .)
#check assumptions - as above use plots and not test because of small sample size
#shapiro.test(residuals(macro_lm))#P greater than 0.05-  normally distributed
#leveneTest(Density_cm_max ~ Species_updated, data = macro3) #p greater than 0.05 no significant difference between variances across groups
hist(nereo_lm$residuals,breaks=20)
ggqqplot(residuals(macro_lm))
plot(macro_lm, 1) # look at residuals to check the homogeneity of variances. looks bad even though test levene tests good
#Not looking good lets try log transformation

#Use log values
macro_lm <- macro3 %>%   lm(Density_cm_log ~ Species_updated , data = .)
#Check assumptions - as above use plots and not test because of small sample size
#shapiro.test(residuals(macro_lm))#P greater than 0.05-  normally distributed- looks better
#leveneTest(Density_cm_max ~ Species_updated, data = macro3) #p less than 0.05 significant difference between variances across groups
hist(nereo_lm$residuals,breaks=20)
ggqqplot(residuals(macro_lm))
plot(macro_lm, 1) # look at residuals to check the homogeneity of variances.
#Looks better in terms of normalitiy but pattern to residuals with two bell shapes towards center

#Use WHITE.Adjust ANOVA to account for unequal variances
#https://stats.stackexchange.com/questions/91872/alternatives-to-one-way-anova-for-heteroskedastic-data
macro_aov<-Anova(macro_lm,  white.adjust=TRUE) #Type I anova
macro_aov


#write results
write.csv(macro_aov,file="./kelp/results/macro_anova_results.csv")


#Posthoc Tukeys procedure- We have unbalanced design and unequal variance - can still use this procedure of Turkeys
# start with using Glht as good for unbalanced designs
macro_tukey<-glht(model=macro_lm, linfct = mcp(Species_updated = "Tukey"), vcov = vcovHC)
summary(macro_tukey)
cld(macro_tukey,alpha = .05, Letters = letters) # returns the letters
# warning message returned

#lets try Emmeans- which is comparable and similar method for unbalanced designs
#https://stackoverflow.com/questions/62963945/glht-and-emmeans-returning-crazy-compact-letters-for-unbalanced-dataset-in-r
macro_emmeans <- emmeans(macro_lm, "Species_updated", vcov = vcovHC)
summary(macro_emmeans)
(macro_out2<-pairs(macro_emmeans)) #default Tukey method
(macro_cld<-cld(macro_emmeans, Letter = letters, ))
#same results as GLHT posthoc tukey test

#write results
write.csv(macro_out2,file="./kelp/results/macro_tukey_results.csv")


#add letters to data set
(macro_tukey_letters<-macro_cld %>% mutate(Kelp_species="Macrocystis") %>% 
    dplyr::rename("letter"=".group")%>%
    dplyr::select(Kelp_species,Species_updated,letter))

(tukey_letters<-rbind(tukey_letters,macro_tukey_letters))

tukey_letters<-tukey_letters %>% mutate(letter=gsub(" ","",letter)) # remove spaces around letters


#################
#   Graphing   ##
#################

Species_xlabs<-c("Bossiella orbigniana" = "Bossiella\norbigniana", "Bossiella schmittii" = "Bossiella\nschmittii", 
                 "Chiharaea americana f. americana" = "C. americana\nf. americana","Chiharaea americana f. bodegensis" = "C. americana\nf. bodegensis",
                 "Crusticorallina 1" = "Crusticorallina\nmorph 1","Crusticorallina 2" = "Crusticorallina\nmorph 2",
                 "Lithophyllum Corsp27BC"="Lithophyllum\nCorsp27BC" ,"Crusticorallina painei"= "Crusticorallina\npainei","Lithophyllum Corsp4BC"= "Lithophyllum\nCorsp4BC",
                 "Crusticorallina muricata"= "Crusticorallina\nmuricata","Lithothamnion glaciale"= "Lithothamnion\nglaciale",
                 "Calliarthron tuberculosum" = "Calliarthron\ntuberculosum","Cover slip" = "Control\nCover slip", "Hildinbrandia" = "Non-corallline\nred crust",
                 "unknown crust 1" = "unknown\ncrust", "Rock" = "Bare Rock","Dead Crusticorallina" = "Dead\nCrusticorallina",
                 "Hildenbrandia spp"="Hildenbrandia spp\n(non-coralline)","Peyssonnelia spp"="Peyssonnelia spp\n(non-coralline)")


all_d_log<-spore_wide%>% 
  filter(Species_updated_count >= 3)%>% 
  filter(Species_updated!="Cover slip" ) %>% # Remove cover slip- not intended as control based on reveiwer comments
  summarySE(., measurevar="Density_cm_log", groupvars=c("Kelp_species","Species_updated"), na.rm=TRUE)
all_d_log<-left_join(all_d_log,tukey_letters,by=c("Kelp_species","Species_updated")) #add letters to dataset
all_d_log

p_log<-ggplot(data=all_d_log, aes(x=Species_updated, y= Density_cm_log, 
                                  ymin=Density_cm_log -se, ymax= Density_cm_log +se))+
  geom_bar(stat="identity",colour="black",fill="grey60") +
  geom_errorbar( width=0.25) +
  scale_x_discrete(labels=Species_xlabs)+
  geom_text(aes(label=letter, y=Density_cm_log +se), vjust=-1) +
  facet_wrap(~Kelp_species,ncol=1,scales = "free_y", 
             labeller = labeller(Kelp_species = 
                                   c("Costaria Costata" = "A) Costaria costata", 
                                     "Macrocystis" = "B) Macrocystis pyrifera",
                                     "Nereocystis" = "C) Nereocystis luetkeana")),
             strip.position="top")+
  geom_blank(aes(y = (Density_cm_log +se)+0.3*(Density_cm_log +se)))+ # used to push the y axis up
  theme_few()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text=element_text(size=12), axis.title.y=element_text(size=14))+
  theme(strip.text.x = element_text(size = 16))+
  ylab(expression(paste("Kelp Sporophyte density (log   ", cm^{-2}, ")"))) + xlab(NULL)
p_log
ggsave(p_log, width = 21, height = 24, units = "cm", dpi = 2400,
       filename = "./figures/kelp_settlement_pannel_plot_letters.png")


#Easier to add data relative abundances to add to figure later - Update names and italicize later
#write relative abundances to add after
abund<-relative_abund %>% arrange(Kelp_relative_abun)
write.csv(abund,"./kelp/results/relative_abund_coralline_kelp.csv")

str(spore_wide)

#Lets write all species IDS with mean densities to a Excel file
mean_values_all<-spore_wide %>% 
  filter(!is.na( Species_updated)) %>%
  summarySE(., measurevar="Density_cm_max", groupvars=c("Kelp_species","Species_updated"), na.rm=TRUE)
mean_values_all

write.csv(mean_values_all,"./kelp/results/kelp_species_means_counts.csv")
