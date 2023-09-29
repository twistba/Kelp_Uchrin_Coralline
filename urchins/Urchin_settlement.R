## --- Kelp & Urchin Coralline Settlement ---
## Urchin Metamorphosis
## Author: Brenton Twist
## Date Created: 2021-03-05
#-------------------------------------------------

# This script contains code for loading and analysis urchin larvae settlment/metamorphosis data

######################
#   Load packages   ##
######################

library(tidyverse)
library(ggthemes)
library(Rmisc)
library(car)
library(cowplot)
library(grid)
library(lme4)
library(lmerTest)
library(lubridate)
library(multcomp) # removes select function from dplyr
library(emmeans)
library(ggpubr)
library(magrittr)
library(sandwich)
library(rstatix)
library(ggsignif)



##################
#   Load data   ##
##################

#   The following code loads the kelp spore settlement data

#read the metadata file with row numbers for each coralline species
urchin<-read_csv("./urchins/Urchin_metamorphosis.csv")

urchin<- urchin %>% mutate_at(vars(Block,Type,Species_genetic_only,Cor_type,
                                   Morpho_species,Genus, Species_code, Time_period,Recorder, Treated_species),
                              factor)
str(urchin)


##################
#   Tidy data   ##
##################
##time differences in hours between start and finish
# Not used- but using time adjusted metamorphosis does not affect results # metomorphis rates might not be linear with time
urchin<- urchin %>% mutate(#Start_date = dmy(Start_date), Date= dmy(Date), 
                           Start_date_time= ymd_hms(paste(Start_date, Start_time),tz="Canada/Pacific"),
                           End_date_time= ymd_hms(paste(Date, Time),tz="Canada/Pacific"),
                           Time_difference= paste(End_date_time-Start_date_time)
                           )
urchin$Time_difference<-as.numeric(urchin$Time_difference)
str(urchin)



###Creating merged categories based on counts
# Lets just look at those that were completely metmorphosed vs unmeatmorphosed
urchin<- urchin %>% mutate(Meta_complete=(Meta_coralline + Meta_jar),  Not_meta=(Meta_undergoing + Free_swimming) )

urchin_raw<-urchin



#Convert from long to wide to look at errors of reduced metamorphose after 48hrs compared to 24hrs
urchin_wide<- urchin %>% dplyr::select(Jar,Time_period, Meta_complete,Meta_undergoing,
                                       Free_swimming, Not_meta, Total_urchins) %>%
  gather( condition, measurement, Meta_complete:Total_urchins, factor_key=TRUE)%>%
  mutate(condition= paste(condition,Time_period,sep="_"))%>% dplyr::select(Jar, condition, measurement)%>%
  spread(condition, measurement)


urchin_wide_errors<- urchin_wide %>% mutate( Meta_complete_diff= Meta_complete_48hrs- Meta_complete_24hrs,
                                     Max_count=pmax(Total_urchins_24hrs,Total_urchins_48hrs),
                                     Missing_48hrs= Max_count- Total_urchins_48hrs ,
                                     Meta_complete_errors= Meta_complete_diff + Missing_48hrs
                                       #a decrease in metamorphose from24 to 48 that can't be accounted for differences in total coun
                                     )
(large_errors<- urchin_wide_errors %>% filter(Meta_complete_errors<=-2) %>% 
  dplyr::select( Jar,Meta_complete_24hrs, Meta_complete_48hrs,Meta_complete_diff, Max_count, 
          Missing_48hrs, Meta_complete_errors))


##Big issues with B03, C04, I05 and N03
#Lets remove these corallines with big issues and update the 48hrs metamorphosed to max value out of 24 and 48hrs

urchin_wide<-urchin_wide %>% 
  filter(!Jar %in% large_errors$Jar) %>% #Remove those with big issues
  mutate(Meta_complete_adj_24hrs= Meta_complete_24hrs,
         Meta_complete_adj_48hrs= pmax(Meta_complete_48hrs, Meta_complete_24hrs),
         Total_urchins_adj_24hrs= Meta_complete_24hrs +Not_meta_24hrs,
          Total_urchins_adj_48hrs= Meta_complete_48hrs +Not_meta_48hrs 
         )

str(urchin_wide)



urchin_long<- urchin_wide  %>%
  gather( condition, measurement,Free_swimming_24hrs:Total_urchins_adj_48hrs, 
                                      factor_key=TRUE)%>% 
  tidyr::extract(condition, into = c("condition", "Time_period"), "(.*)_([^_]+)$") %>%
  spread(condition, measurement)  %>% 
  dplyr::select(Jar, Time_period, Meta_complete_adj,Total_urchins_adj)

str(urchin_long)

urchin <-  left_join(urchin_long, urchin_raw, by=c("Jar","Time_period"))
str(urchin)


##Calculate percentage metamorphosed
# It doesn't make sense to calculate the number metamorphosed per hour as this could not be linear
# and plus accounted for time difference by the block designs
urchin <- urchin %>% mutate(Meta_per= (Meta_complete/Total_urchins)*100,
                            Meta_adj_per=(Meta_complete_adj/Total_urchins_adj )*100,
                            Meta_per_time= (Meta_per/Time_difference),
                            Meta_adj_per_time=(Meta_adj_per/Time_difference)
                            )

str(urchin)


#Add a species count of each genetically ID coralline species-  so can filter

urchin_species_counts_gen<- urchin%>%                  
  filter(!is.na(Species_genetic_only)) %>%    
  group_by(Time_period,Type,Species_genetic_only) %>%        
  dplyr::summarise(Species_genetic_only_count = n()) 


urchin<-left_join( urchin, urchin_species_counts_gen, by=c("Time_period","Type","Species_genetic_only"))


## add relative abundance numbers
load(file="./coralline_relative-abund/relative_coralline.Rdata") # load in relative abundance
cor_rel_abun<-cor_rel_abun %>% dplyr::rename(Species_genetic_only=Species_name_martone)

(relative_abund<-left_join(data.frame(Species_genetic_only=unique(urchin$Species_genetic_only)),
                           cor_rel_abun, by="Species_genetic_only"))

relative_abund %<>% mutate(Kelp_relative_abun = round(Kelp_relative_abun, 2),
                           Urchin_relative_abun = round(Urchin_relative_abun, 2))  #round to 2 dp

fleshy<-cor_rel_abun %>% dplyr::filter(Species_genetic_only=="Fleshy red crust")
kelp_fleshy<-round(fleshy$Kelp_relative_abun,2)
urchin_fleshy<-round(fleshy$Urchin_relative_abun,2)

#lets add a relative abundance for two ACA and fleshy reds
str(relative_abund)

relative_abund %<>% 
  mutate(Kelp_relative_abun=replace(Kelp_relative_abun, 
                                      Species_genetic_only=="Calliarthron tuberculosum", 1.00),
         Urchin_relative_abun=replace(Urchin_relative_abun, 
                                        Species_genetic_only=="Calliarthron tuberculosum", 0.00),
         Kelp_relative_abun=replace(Kelp_relative_abun, 
                                    Species_genetic_only=="Hildenbrandia spp"| 
                                      Species_genetic_only=="Peyssonnelia spp", kelp_fleshy),
         Urchin_relative_abun=replace(Urchin_relative_abun, 
                                      Species_genetic_only=="Hildenbrandia spp"| 
                                        Species_genetic_only=="Peyssonnelia spp", urchin_fleshy),
         Kelp_relative_abun=replace(Kelp_relative_abun, 
                                    Species_genetic_only=="Corsp23BCcrust", 0.24),
         Urchin_relative_abun=replace(Urchin_relative_abun, 
                                      Species_genetic_only=="Corsp23BCcrust", 0.76)
         ) 

(abund<-relative_abund %>% arrange(Kelp_relative_abun))
write.csv(abund,"./urchins/results/relative_abund_coralline_urchin.csv")


urchin<-left_join(urchin,relative_abund, by="Species_genetic_only")


##Make species a factor and put in order of relative abundance
urchin$Species_genetic_only <- factor(urchin$Species_genetic_only, 
                                      levels=c("Control",  "Hildenbrandia spp", "Peyssonnelia spp","Lithothamnion sp Rhodolith", 
                                               "Lithophyllum Corsp27BC" ,"Lithophyllum sp4","Corsp23BCcrust","Crusticorallina muricata",
                                               "Crusticorallina painei", "Chiharaea americana f. bodegensis", "Lithothamnion glaciale",
                                               "Lithophyllum Corsp4BC", "Corsp7BCcrsut",  "Bossiella schmittii", 
                                               "Calliarthron tuberculosum"  ))
levels(urchin$Species_genetic_only)


##############################
#  Untreated species only   ##
##############################

#### --- Statistical analysis --- ####
#-------------------------------------

#### -- 24hrs -- ####

urchin3_24hrs<-urchin %>% drop_na(Species_genetic_only)%>%
  filter(Type =="Untreated") %>%
  filter(Time_period =="24hrs") %>%
  filter(Species_genetic_only_count >= 3) # removes species with less than 3 replicates & those without densities
str(urchin3_24hrs)

urchin_lm_24hrs <- urchin3_24hrs %>% lm(Meta_adj_per ~ Species_genetic_only, data = .)
# check assumptions - use plots and not test because of small sample size 
# https://stats.stackexchange.com/questions/2492/is-normality-testing-essentially-useless
  #shapiro.test(residuals(urchin_lm_24hrs))#P greater than 0.05- normally distributed
  #leveneTest(Meta_adj_per ~ Species_genetic_only, data = urchin3_24hrs) #p greater 0.05 no significant difference between variances across groups
hist(urchin_lm_24hrs$residuals,breaks=20)
ggqqplot(residuals(urchin_lm_24hrs)) #look at qqplot to see normality 
plot(urchin_lm_24hrs, 1) # look at residuals to check the homogeneity of variances.
#doesn't look like homogeneity of variances try log transformation

#Log transformed model
urchin_log_lm_24hrs <- urchin3_24hrs %>% lm(log(Meta_adj_per+1) ~ Species_genetic_only , data = .)
# check assumptions
  #shapiro.test(residuals(urchin_log_lm_24hrs))#P less than 0.05-  not normally distributed
  #leveneTest(log(Meta_adj_per+1) ~ Species_genetic_only, data = urchin3_24hrs) #p less than 0.05 unequal variance
hist(urchin_log_lm_24hrs$residuals,breaks=20)
ggqqplot(residuals(urchin_log_lm_24hrs)) #look at qqplot to see normality 
plot(urchin_log_lm_24hrs, 1) # look at residuals to check the homogeneity of variances.
# Makes it worse than pre-transformed. Switch to back to non-transformed as it is normally distributed

###Use WHITE.Adjust ANOVA to account for unequal variances
#https://stats.stackexchange.com/questions/91872/alternatives-to-one-way-anova-for-heteroskedastic-data

urchin_aov_24hrs<-Anova(urchin_lm_24hrs,  white.adjust=TRUE) #Type I anova
urchin_aov_24hrs


#write results
write.csv(urchin_aov_24hrs,file="./urchins/results/urchin_24hrs_anova_results.csv")

#Posthoc Tukeys procedure- We have unbalanced design and unequal variance - can still use this procedure of Turkeys
#adding heteroskedasticity-consistent covariance matrix estimation (vcov = vcovHC) following Herberich et al. 2010
#lets start with Emmeans due to error with GLHT- which is comparable and similar method for unbalanced designs
urchin_emmeans_24hrs <- emmeans(urchin_lm_24hrs, "Species_genetic_only",  vcov = vcovHC) 
summary(urchin_emmeans_24hrs)
(urchin_out1<-pairs(urchin_emmeans_24hrs)) #default Tukey method
(urchin_cld_24hrs<-cld(urchin_emmeans_24hrs, Letter = letters) )

#lets just compare to GLHT to see if same results despite probable warning
urchin_glht_24hrs<-glht(model=urchin_lm_24hrs, linfct = mcp(Species_genetic_only = "Tukey"),  vcov = vcovHC)
summary(urchin_glht_24hrs)
cld(urchin_glht_24hrs,alpha = .05, Letters = letters) # returns the letters
# warning message returned
# same results as emmeans posthoc tukey test

#write results
write.csv(urchin_out1,file="./urchins/results/urchin_24hrs_tukey_results.csv")


#add letters to data set
(urchin_24hrs_tukey_letters<-urchin_cld_24hrs %>% mutate(Type ="Untreated",Time_period ="24hrs") %>% 
    dplyr::rename("letter"=".group")%>%
    dplyr::select(Type,Time_period,Species_genetic_only,letter))

(tukey_letters<-urchin_24hrs_tukey_letters)


####  -- 48hrs -- ####

urchin3_48hrs<-urchin %>% drop_na(Species_genetic_only)%>%
  filter(Type =="Untreated") %>%
  filter(Time_period =="48hrs") %>%
  filter(Species_genetic_only_count >= 3) # removes species with less than 3 replicates & those without densities
str(urchin3_48hrs)

urchin_lm_48hrs <- urchin3_48hrs %>% lm(Meta_adj_per ~ Species_genetic_only, data = .)
# check assumptions - use plots and not test because of small sample size 
# https://stats.stackexchange.com/questions/2492/is-normality-testing-essentially-useless
  #shapiro.test(residuals(urchin_lm_48hrs))#P greater than 0.05- normally distributed
  #leveneTest(Meta_adj_per ~ Species_genetic_only, data = urchin3_48hrs) #p less than 0.05 significant difference between variances across groups
hist(urchin_lm_48hrs$residuals,breaks=20)
ggqqplot(residuals(urchin_lm_48hrs)) #look at qqplot to see normality 
plot(urchin_lm_48hrs, 1) # look at residuals to check the homogeneity of variances.
#doesn't look like homogeneity of variances try log transformation

#Log transformed model
urchin_log_lm_48hrs <- urchin3_48hrs %>% lm(log(Meta_adj_per+1) ~ Species_genetic_only , data = .)
# check assumptions
  #shapiro.test(residuals(urchin_log_lm_48hrs))#P less than 0.05-  not normally distributed
  #leveneTest(log(Meta_adj_per+1) ~ Species_genetic_only, data = urchin3_48hrs) #p greater than 0.05 equal variance
hist(urchin_log_lm_48hrs$residuals,breaks=20)
ggqqplot(residuals(urchin_log_lm_48hrs)) #look at qqplot to see normality 
plot(urchin_log_lm_48hrs, 1) # look at residuals to check the homogeneity of variances.
# Makes it worse than pre-transformed. Switch to back to non-transformed 


###Use WHITE.Adjust ANOVA to account for unequal variances
#https://stats.stackexchange.com/questions/91872/alternatives-to-one-way-anova-for-heteroskedastic-data

urchin_aov_48hrs<-Anova(urchin_lm_48hrs,  white.adjust=TRUE) #Type I anova
urchin_aov_48hrs

#write results
write.csv(urchin_aov_48hrs,file="./urchins/results/urchin_48hrs_anova_results.csv")

#Posthoc Tukeys procedure- We have unbalanced design and unequal variance - can still use this procedure of Turkeys
#adding heteroskedasticity-consistent covariance matrix estimation (vcov = vcovHC) following Herberich et al. 2010
#lets start with Emmeansdue to error with GLHT- which is comparable and similar method for unbalanced designs
urchin_emmeans_48hrs <- emmeans(urchin_lm_48hrs, "Species_genetic_only",  vcov = vcovHC) 
summary(urchin_emmeans_48hrs)
(urchin_out2<-pairs(urchin_emmeans_48hrs)) #default Tukey method
(urchin_cld_48hrs<-cld(urchin_emmeans_48hrs, Letter = letters) )

#lets just compare to GLHT to see if same results despite probable warning
urchin_glht_48hrs<-glht(model=urchin_lm_48hrs, linfct = mcp(Species_genetic_only = "Tukey"),  vcov = vcovHC)
summary(urchin_glht_48hrs)
cld(urchin_glht_48hrs,alpha = .05, Letters = letters) # returns the letters
# warning message returned
# same results as emmeans posthoc tukey test

#write results
write.csv(urchin_out2,file="./urchins/results/urchin_48hrs_tukey_results.csv")


#add letters to data set
(urchin_48hrs_tukey_letters<-urchin_cld_48hrs %>% mutate(Type ="Untreated",Time_period ="48hrs") %>% 
    dplyr::rename("letter"=".group")%>%
    dplyr::select(Type,Time_period,Species_genetic_only,letter))

(tukey_letters<-rbind(tukey_letters,urchin_48hrs_tukey_letters))


tukey_letters<-tukey_letters %>% mutate(letter=gsub(" ","",letter)) # remove spaces around letters





#### --- Graphing untreated --- ####
#-----------------------------------

Species_xlabs<-c("Bossiella schmittii" = "Bossiella\nschmittii", "Chiharaea americana f. bodegensis" = "C. americana\nf. bodegensis",
                 "Lithophyllum Corsp27BC"="Lithophyllum\nCorsp27BC" ,"Crusticorallina painei"= "Crusticorallina\npainei",
                 "Lithothamnion sp Rhodolith" = "Lithothamnion sp\nRhodolith" ,"Lithothamnion glaciale"= "Lithothamnion\nglaciale",
                 "Calliarthron tuberculosum" = "Calliarthron\ntuberculosum","Lithophyllum Corsp4BC" = "Lithophyllum\nCorsp4BC",
                 "Hildenbrandia spp"="Hildenbrandia spp\n(non-coralline)","Peyssonnelia spp"= "Peyssonnelia spp\n(non-coralline)")



urchin_plot<-urchin%>% drop_na(Species_genetic_only)%>%
  filter(Species_genetic_only_count >= 3)%>% # selects species 3 or more occurances
  summarySE(., measurevar="Meta_adj_per", groupvars=c("Type","Time_period","Species_genetic_only"), na.rm=TRUE)

urchin_plot<-left_join(urchin_plot,tukey_letters,
                         by=c("Type","Time_period","Species_genetic_only")) #add letters to dataset
urchin_plot

(urchin_plot_untreat<- urchin_plot %>% filter(Type =="Untreated"))



p1<-ggplot(data=urchin_plot_untreat, aes(x=Species_genetic_only, y= Meta_adj_per, 
                                  ymin=Meta_adj_per -se, ymax= Meta_adj_per +se))+
  geom_bar(stat="identity",colour="black",fill="grey60") +
  geom_errorbar( width=0.25) +
  scale_x_discrete(labels=Species_xlabs)+
  geom_text(aes(label=letter, y=Meta_adj_per +se), vjust=-1) +
  facet_wrap(~Time_period,ncol=1,
#             scales = "free_y", 
             labeller = labeller(Time_period = 
                                   c("24hrs" = "A) 24 Hours", 
                                     "48hrs" = "B) 48 Hours")),
             strip.position="top") +
  geom_blank(aes(y = (Meta_adj_per +se)+0.15*(Meta_adj_per +se)))+ # used to push the y axis up
  theme_few()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text=element_text(size=12), axis.title.y=element_text(size=14))+
  theme(strip.text.x = element_text(size = 16))+
  ylab("Metamorphosed Urchin Larvae (%)") + xlab(NULL)
p1
ggsave(p1, width = 21, height = 16, units = "cm", dpi = 2400,
       filename = "./figures/Urchin_settlement_pannel_plot_letters.png")

##Add relative abundances and fix up names later


############################
#  Antibiotic treatment   ##
############################

#### --- cleaning data --- ####
#------------------------------

##difference between treated and Untreated
#Find species in both untreated and antibiotic to compare
sp_untreated<-urchin %>% filter(Type =="Untreated") %>% drop_na(Species_genetic_only)%>%
  filter(Species_genetic_only_count >= 3)%>% 
  distinct(Species_genetic_only) %>% dplyr::select(Species_genetic_only)

sp_antibiotic<-urchin %>% filter(Type =="Antibiotic") %>%  drop_na(Species_genetic_only)%>%
  filter(Species_genetic_only_count >= 3)%>% 
  distinct(Species_genetic_only) %>% dplyr::select(Species_genetic_only)

#select species only in both datasets 
(sp_both<-inner_join(sp_untreated,sp_antibiotic, by="Species_genetic_only")) #species in both data sets

urchin_antiboitic<-left_join(sp_both,urchin, by="Species_genetic_only")


#### --- Statistical analysis --- ####
#-------------------------------------

#24hours t-test between antibiotic and untreated for given species
urchin_antiboitic_24hrs<-  urchin_antiboitic%>% filter(Time_period =="24hrs")
  
t_urchin_24 <- urchin_antiboitic_24hrs %>%
  group_by(Species_genetic_only) %>%
  t_test(Meta_adj_per ~ Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
t_urchin_24 #all non-significant

write.csv(t_urchin_24,file="./urchins/results/urchin_antibiotics_24hrs_t-test.csv")

#24hours t-test between antibiotic and untreated for given species
urchin_antiboitic_48hrs<-  urchin_antiboitic%>% filter(Time_period =="48hrs")

t_urchin_48 <- urchin_antiboitic_48hrs %>%
  group_by(Species_genetic_only) %>%
  t_test(Meta_adj_per ~ Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
t_urchin_48 #all non-significant

write.csv(t_urchin_48,file="./urchins/results/urchin_antibiotics_48hrs_t-test.csv")





#### --- Graphing antibiotic --- ####
#-------------------------------------

#use the dataframe created in the graphing before with mean and SE
(urchin_plot_anti<-left_join(sp_both,urchin_plot, by="Species_genetic_only"))
 

p2<-ggplot(data=urchin_plot_anti, aes(x=Species_genetic_only, y= Meta_adj_per,fill=Type, 
                                      ymin= Meta_adj_per - se, ymax= Meta_adj_per + se))+
  geom_bar(position="dodge", stat="identity",colour="black")+
  scale_fill_manual(values=c("white","gray")) +
  geom_errorbar( position = position_dodge(0.9),width=0.25)+
  scale_x_discrete(labels=Species_xlabs)+
  facet_wrap(~Time_period,ncol=1,
             #             scales = "free_y", 
             labeller = labeller(Time_period = 
                                   c("24hrs" = "A) 24 Hours", 
                                     "48hrs" = "B) 48 Hours")),
             strip.position="top") +
  theme_few()+ 
  ylab("Metamorphosed Urchin Larvae (%)") + xlab(NULL) + 
  theme(legend.title=element_blank(), legend.justification=c(0,1), legend.position=c(0.04,0.98), 
        legend.text = element_text(size = 14))+
  theme(axis.text=element_text(size=12), axis.title.y=element_text(size=14))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(strip.text.x = element_text(size = 16))
p2

ggsave(p2, width = 21, height = 16, units = "cm", dpi = 2400,
       filename = "./figures/Urchin_antibiotic_settlement_pannel_plot_BW.png")

#lest add significance bars
#need to get heights of max + se
(sig_height_long<-urchin_plot_anti %>% 
  dplyr::mutate(y_height=Meta_adj_per+se, Type_time=paste(Type,Time_period,sep="_")) %>%
  dplyr::select(Species_genetic_only,Type_time,y_height))
  
sig_height_wide <- spread(sig_height_long, Type_time, y_height)
sig_height_wide 

(sig_height_wide %<>% rowwise() %>%
  mutate(max_val = max(c_across(where(is.numeric))))) #gets max value
#Coralline species seem to be in right order


p2a<- p2 + geom_signif(annotations = rep("ns",7),
                 y_position = sig_height_wide$max_val+3, xmin=seq(0.8, 6.8, 1), xmax=seq(1.2, 7.2, 1))+
  geom_blank(aes(y = (Meta_adj_per +se)+0.1*(Meta_adj_per +se)))
p2a

ggsave(p2a, width = 21, height = 16, units = "cm", dpi = 2400,
       filename = "./figures/Urchin_antibiotic_settlement_pannel_plot_BW_sig.png")



#write mean percentage metamorphosis
urchin_per_all<-urchin %>% drop_na(Species_genetic_only)%>%
  summarySE(., measurevar="Meta_adj_per", groupvars=c("Type","Time_period","Species_genetic_only"), na.rm=TRUE)
urchin_per_all

write.csv(urchin_per_all,"./urchins/results/urchin_species_means_percentage.csv")
write.csv(urchin_plot_anti,"./urchins/results/urchin_antibiotics_species_means_percentage.csv")

