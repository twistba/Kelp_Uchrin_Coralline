## --- Kelp & Urchin Coralline Settlement ---
## Microbiome - Beta diversity analysis
## Author: Florent Mazel 
## Edited by: Brenton Twist
## Date Created: 2021-03-05   
#-------------------------------------------------


# This script contains code for analyzing  the micriobiome community structure


######################
#   Load packages   ##
######################

library(tidyverse)
library(phyloseq)
library(ggConvexHull)
library(MicEco) # for omega2
library(RColorBrewer)
library(vegan)
library(ggthemes)
library(usedist)
library(ggpubr)
library(ape)

#Install microbiome package
# library(BiocManager)
# BiocManager::install("microbiome")


##################
#   Load data   ##
##################

phyloseq_obj_rare = readRDS("./microbiome/Phyloseq_object_rarefied.1000.RDS")
phyloseq_obj_nonR = readRDS("./microbiome/Phyloseq_object_non_rarefied.RDS")

sample_info_rare = as_tibble(sample_data(phyloseq_obj_rare))
sample_info_nonR = as_tibble(sample_data(phyloseq_obj_nonR))


##########################
#  Supp. Mat. figure    ##
##########################


## Rarefied approach
# Bray–Curtis
bray_div = phyloseq_obj_rare %>% 
  phyloseq::distance("bray")
hist(c(bray_div))


## Non Rarefied 
# Compositionally -aware method; we can use the « Aitchison » distance 
aitc_div_nonR = phyloseq_obj_nonR %>% 
  microbiome::transform("clr") %>% 
  phyloseq::distance("euclidean")
hist(c(aitc_div_nonR))



# NMDS for rarefied Bray-Curtis data
set.seed(258)
MDS_bray <- bray_div %>% 
  metaMDS() 

NMDS_coord = MDS_bray$points %>% as_tibble(rownames = "swab.ID") %>% 
  dplyr::left_join(sample_info_rare)

table(sample_info_rare$Substrate)

# Non-Rarefied CoDa approach
set.seed(258)
PCOA_aitc_nonR <-aitc_div_nonR %>% 
  pcoa() 
PCOA_coord_nonR = PCOA_aitc_nonR$vectors[,1:3]%>%
  as_tibble(rownames = "swab.ID") %>% 
  dplyr::left_join(sample_info_nonR)


#PERMANOVA bray 
set.seed(245)
model_beta = adonis2(bray_div~Substrate, data = as_tibble(sample_data(phyloseq_obj_rare)),
                     permutations = 9999)
model_beta2 = adonis_OmegaSq(model_beta)
model_beta2
# Df SumOfSqs      F parOmegaSq Pr(>F)    
# Substrate  5   6.2411 4.5595    0.25499  0.001 ***
# Residual  46  12.5933                             
# Total     51  18.8344   

#PERMANOVA aitchison 
set.seed(245)
model_betaAi = adonis2(aitc_div_nonR~Substrate, data = as_tibble(sample_data(phyloseq_obj_rare)),
                     permutations = 9999)
model_betaAi2 = adonis_OmegaSq(model_betaAi)
model_betaAi2
#Df SumOfSqs     F parOmegaSq Pr(>F)    
#Substrate  5    58757 4.364    0.24441  1e-04 ***
#Residual  46   123867                            
#Total     51   182624  


# FIGURE 

#Coralline spelt wrong

NMDS_coord<-NMDS_coord%>%
  mutate(Substrate= str_replace(Substrate,"Corraline", "Coralline"))


Fig1a <- NMDS_coord %>% 
  ggplot(aes(y=MDS2,x=MDS1,shape=Substrate)) + #ylim(-.3,.3)+
  geom_point(size=2)+
  guides(shape=guide_legend(ncol=2),fill=guide_legend(ncol=2))+
  theme_few()+ 
  annotate("text",x=-0.18, y=.8, label=paste0("Stress value  = ",round(MDS_bray$stress,2))) +
theme(#legend.spacing.y = unit(0, "mm"), 
      panel.border = element_rect(colour = "black", fill=NA),
      legend.background = element_blank(),
      #legend.box.background = element_rect(colour = "black"),
      legend.position = c(0.83, 0.87),legend.title=element_text(size=14),
      legend.text=element_text(size=12)
      )

Fig1a

Fig1a_col<- Fig1a +  geom_convexhull(aes(fill=Substrate),alpha=0.15,size=.2)+
  #scale_fill_manual(values= c("#67ae62","white","white", "#7f63b8", "#b87e39","#b94971"))
  scale_fill_manual(values= c("red","white","white", "#b0923b", "#006400","blue"))
Fig1a_col
#Fig1_bw<- Fig1 +  geom_convexhull(alpha=0.05,size=.2)
#Fig1_bw

Fig1a_col2<- Fig1a_col +theme(legend.box.background = element_rect(colour = "black"))
Fig1a_col2


ggsave(Fig1a_col2, width = 21, height = 16, units = "cm", dpi = 600,
      filename = "./figures/Microbiome_multivariate_corallines_controls_Bray_NMDS_rarefied.png")

PCOA_coord_nonR<-PCOA_coord_nonR%>%
  mutate(Substrate= str_replace(Substrate,"Corraline", "Coralline"))

Fig1b <- PCOA_coord_nonR %>% 
  ggplot(aes(y=Axis.2,x= Axis.1,shape=Substrate)) + #ylim(-.3,.3)+
  geom_point(size=2)+
  xlab("Axis I")+ ylab("Axis II")+
  guides(shape=guide_legend(ncol=2),fill=guide_legend(ncol=2))+
  theme_few()+ 
  #annotate("text",x=-0.3, y=.8, label=paste0("Stress value  = ",round(MDS_bray$stress,2))) +
  theme(#legend.spacing.y = unit(0, "mm"), 
    panel.border = element_rect(colour = "black", fill=NA),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    legend.position = c(0.83, 0.83),legend.title=element_text(size=14),
    legend.text=element_text(size=12)
  )

Fig1b

Fig1b_col<- Fig1b +  geom_convexhull(aes(fill=Substrate),alpha=0.15,size=.2)+
  #scale_fill_manual(values= c("#67ae62","white","white", "#7f63b8", "#b87e39","#b94971"))
  scale_fill_manual(values= c("red","white","white", "#b0923b", "#006400","blue"))
Fig1b_col
#Fig1b_bw<- Fig1b +  geom_convexhull(alpha=0.05,size=.2)
#Fig1b_bw



ggsave(Fig1b_col, width = 21, height = 16, units = "cm", dpi = 600,
       filename = "./figures/Microbiome_multivariate_corallines_controls_PCoA_Aitchison_Non_rarefied.png")


#### -- Combined figure -- ####

Fig1a_comb<-Fig1a_col +  theme(legend.position="top", legend.title = element_blank())+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
Fig1a_comb

Fig1b_comb<-Fig1b_col +  theme(legend.position="none")
Fig1b_comb



Fig1c <- ggarrange(Fig1a_comb, Fig1b_comb, 
                  labels = c("A)", "B)"),
                  ncol = 1, nrow = 2,
                  heights = c(1,0.8))
Fig1c

ggsave(Fig1c, width = 21, height = 24, units = "cm", dpi = 600,
       filename = "./figures/Microbiome_multivariate_corallines_controls_nMDS_and_PCoA.png")





##########################
#  Main      figure    ##
##########################

#### --- Statistical Analysis ---####
#-------------------------------------

## Rarefied approach
phyloseq_obj_rare_subset = phyloseq_obj_rare %>% subset_samples(Substrate %in% c("Corraline","Red_crust"))
sample_info_rare_sub = as_tibble(sample_data(phyloseq_obj_rare_subset ))

# Bray-Curtis method
bray_div_sub = phyloseq_obj_rare_subset  %>% 
  distance("bray")
hist(c(bray_div_sub))

## Non Rarefied 
# Compositionally -aware method « Aitchison » distance 

phyloseq_obj_nonR_subset = phyloseq_obj_nonR %>% subset_samples(Substrate %in% c("Corraline","Red_crust"))
sample_info_nonR_sub = as_tibble(sample_data(phyloseq_obj_nonR_subset))

aitc_div_nonR_sub = phyloseq_obj_nonR_subset %>% 
  microbiome::transform("clr") %>% 
  phyloseq::distance("euclidean")
hist(c(aitc_div_nonR_sub))


# NMDS for rarefied bray curtis data
set.seed(946)
MDS_bray2 <- bray_div_sub %>% 
  metaMDS() 

NMDS_coord2 = MDS_bray2$points %>% as_tibble(rownames = "swab.ID") %>% 
  left_join(sample_info_rare_sub)

# PCoA for non rarefied compositionally -aware method
set.seed(946)
PCOA_aitc_nonR2 <- aitc_div_nonR_sub %>% 
  pcoa() 

PCOA_coord_nonR2 = PCOA_aitc_nonR2$vectors[,1:3]%>%
  as_tibble(rownames = "swab.ID") %>% 
  dplyr::left_join(sample_info_nonR_sub)


#PERMANOVA rarefied data 
set.seed(253)
model_beta = adonis2(bray_div_sub~Boulder+DNA.ID.updated, 
                     data = as_tibble(sample_data(phyloseq_obj_rare_subset)),
                     by = "margin", permutations = 9999)
model_beta2 = adonis_OmegaSq(model_beta)
model_beta2

#Df SumOfSqs     F parOmegaSq Pr(>F)    
#Boulder        22   4.4434 1.546    0.23095  0.001 ***
#  DNA.ID.updated  8   3.4374 3.289    0.31403  0.001 ***
#  Residual        9   1.1758                            
#Total          39  12.0945  

#PERMANOVA redcrust rarefied data
dataRedcrust = as_tibble(sample_data(phyloseq_obj_rare_subset)) %>% 
  mutate(redcrust=Substrate=="Red_crust") %>% 
  group_by(redcrust) %>% sample_n(5)


set.seed(123)
model_beta = adonis2(dist_subset(bray_div_sub,dataRedcrust$swab.ID)~redcrust, 
                     data = dataRedcrust,by = "margin", permutations = 9999)
model_beta2 = adonis_OmegaSq(model_beta)
model_beta2
#Df SumOfSqs     F parOmegaSq Pr(>F)    
#Boulder        22   4.4434 1.546    0.23095  0.001 ***
#  DNA.ID.updated  8   3.4374 3.289    0.31403  0.001 ***
#  Residual        9   1.1758                            
#Total          39  12.0945  


#PERMANOVA unrarefied data CoDa approach  
set.seed(253)
model_beta_CoDA = adonis2(aitc_div_nonR_sub~Boulder+DNA.ID.updated, 
                     data = as_tibble(sample_data(phyloseq_obj_rare_subset)),
                     by = "margin", permutations = 9999)
model_beta_CoDA2 = adonis_OmegaSq(model_beta_CoDA)
model_beta_CoDA2

#Df SumOfSqs      F parOmegaSq Pr(>F)    
#Boulder        22    42390 1.2792    0.13312 0.0193 *  
#  DNA.ID.updated  8    26183 2.1728    0.19000 0.0001 ***
#  Residual        9    13556                             
#Total          39   102283     


#### --- Family level diversity --- #####
#----------------------------------------

Taxo_counts_Family <-  phyloseq_obj_rare %>%
  tax_glom(taxrank = 'family' ,NArm=F) %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()

Taxonomy = as.tibble(tax_table(phyloseq_obj_rare))
class(Taxonomy)
dim(Taxonomy)



tableC <- table(c(Taxonomy$family))
ASV_c <- tibble(family=names(tableC),
                ASV_counts = tableC)


Summary_Family = Taxo_counts_Family %>% 
  group_by(family) %>% 
  summarise(Sum_RelAb_bySample=sum(Abundance),
            Mean_RelAb_bySample=mean(Abundance),
            Sd_RelAb_bySample=sd(Abundance)) %>% 
  mutate(Perc_tot_relatoveReadCounts=Sum_RelAb_bySample/sum(Sum_RelAb_bySample)) %>% 
  left_join(ASV_c)

Family_Prev = Taxo_counts_Family %>% 
  group_by(family) %>% 
  summarise(Prev=sum(Abundance>0)) %>%
  subset(Prev<70) %>% 
  pull(family)


Taxo_counts_sample = Taxo_counts_Family %>% 
  mutate(Family_plot = ifelse((Abundance<.005)|(family%in%Family_Prev),"Rare Family",family))


#### --- Graphing --- ####
#-------------------------
 
#### -- Multivariate plot --####
Species_labs<-c("Bossiella exarticulata" = "*Bossiella exarticulata*" , "Bossiella mayae" = "*Bossiella mayae*",
                "Chiharaea americana"= "*Chiharaea americana*" , "Corsp23BC"="Hapalidiales sp. 1",
                "Crusticorallina muricata"= "*Crusticorallina muricata*","Crusticorallina painei"= "*Crusticorallina painei*",
                "Lithophyllum impressum" = "*Lithophyllum impressum*", "Lithophyllum sp4" = "*Lithophyllum* sp. 3",
                "Peyssonnelia sp10BC"  =  "*Peyssonnelia* sp")
ordered(unique(NMDS_coord2$DNA.ID.updated))

str(NMDS_coord2)
NMDS_coord2<- NMDS_coord2 %>% dplyr::rename(Species=DNA.ID.updated)


Fig2a_col <- NMDS_coord2 %>% 
  ggplot(aes(y=MDS2,x=MDS1,shape=Species)) + 
  #xlim(-.45,.85)+
  #ylim(-0.2,.5)+
  geom_point(alpha=1,size=2.5,stroke=1.1) + 
  geom_convexhull(aes(fill=Species), alpha=0.15,size=.2)+
  scale_shape_manual(values=c(1,7,3,4,15,6,17,8,19), labels=Species_labs)+
  scale_fill_manual(values= c("#F8766D", "#D39200","white", "white" , "#00C19F", "#00B9E3", "#00BA38", "#DB72FB",
                              "#FF61C3"),labels=Species_labs )+
  theme_few() + 
  guides(shape=guide_legend(ncol=2),fill=guide_legend(ncol=2))+
  annotate("text",x=-.35, y=.4, label=paste0("Stress value  = ",round(MDS_bray2$stress,2))) +
  theme(#legend.spacing.y = unit(0, "mm"), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.position = c(0.72, 0.83),legend.title=element_text(size=13),
        legend.text=ggtext::element_markdown(size=11)
  )
Fig2a_col


Fig2a_bw <- NMDS_coord2 %>% 
  ggplot(aes(y=MDS2,x=MDS1,shape=Species)) + 
  #xlim(-.45,.85)+
  #ylim(-0.2,.5)+
  geom_point(alpha=1,size=2.5,stroke=1.1) + 
  geom_convexhull( alpha=0.07,size=.2)+
  scale_shape_manual(values=c(1,7,3,4,15,6,17,8,19), labels=Species_labs) +
  theme_few() + 
  guides(shape=guide_legend(ncol=2),fill=guide_legend(ncol=2))+
  annotate("text",x=-.35, y=.4, label=paste0("Stress value  = ",round(MDS_bray2$stress,2))) +
  theme(#legend.spacing.y = unit(0, "mm"), 
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    legend.position = c(0.72, 0.83),legend.title=element_text(size=13),
    legend.text=ggtext::element_markdown(size=11)
  )
Fig2a_bw



ggsave(Fig2a_col, width = 21, height = 16, units = "cm", dpi = 2400,
       filename = "./figures/Microbiome_multivariate_corallines_Bray_NMDS_Rarefied.png")

ggsave(Fig2a_bw, width = 21, height = 16, units = "cm", dpi = 2400,
       filename = "./figures/Microbiome_multivariate_corallines_Bray_NMDS_Rarefied_bw.png")



#PCOA CoDA unrarefied data
str(PCOA_coord_nonR2)
PCOA_coord_nonR2<- PCOA_coord_nonR2 %>% dplyr::rename(Species=DNA.ID.updated)

#remove random #N/A in dataset. not sure how this entered dataset
#PCOA_coord2<- PCOA_coord2 %>% filter(Species!="#N/A")
ordered(unique(PCOA_coord_nonR2$Species))



Fig2b_col <- PCOA_coord_nonR2 %>% 
  ggplot(aes(y=Axis.2,x=Axis.1,shape=Species)) + 
   ylab("Axis II")+ xlab("Axis I")+
  geom_point(alpha=1,size=2.5,stroke=1.1) + 
  geom_convexhull(aes(fill=Species), alpha=0.15,size=.2)+
  scale_shape_manual(values=c(1,7,3,4,15,6,17,8,19,12), labels=Species_labs)+
  scale_fill_manual(values= c("#F8766D", "#D39200","white", "white" , "#00C19F", "#00B9E3", "#00BA38", "#DB72FB",
                              "#FF61C3"),labels=Species_labs )+
  theme_few() + 
  theme(legend.position="top",legend.title = element_blank(),
        legend.text=ggtext::element_markdown(size=11))+
  #guides(fill=guide_legend(nrow=4,byrow=TRUE))
  guides(shape=guide_legend(nrow=3),fill=guide_legend(nrow=3))





Fig2b_col



Fig2b2_col <- PCOA_coord_nonR2 %>% 
  ggplot(aes(y=Axis.3,x=Axis.1,shape=Species)) + 
  ylab("Axis III")+ xlab("Axis I")+ 
  geom_point(alpha=1,size=2.5,stroke=1.1) + 
  geom_convexhull(aes(fill=Species), alpha=0.15,size=.2)+
  scale_shape_manual(values=c(1,7,3,4,15,6,17,8,19,12), labels=Species_labs)+
  scale_fill_manual(values= c("#F8766D", "#D39200","white", "white" , "#00C19F", "#00B9E3", "#00BA38", "#DB72FB",
                              "#FF61C3"),labels=Species_labs )+
  theme_few() + 
  theme(legend.position = "none")

Fig2b2_col

Fig2bFINAL = ggarrange(Fig2b_col,Fig2b2_col,ncol = 1)
Fig2bFINAL

Fig2bFINAL <- ggarrange(Fig2b_col,Fig2b2_col, 
                   ncol = 1, nrow = 2,
                   heights = c(1,0.68))
Fig2bFINAL

ggsave(Fig2bFINAL, width = 21, height = 28, units = "cm", dpi = 900,
       filename = "./figures/Microbiome_multivariate_corallines_Aitchison_PCOA_Non_Rarefied.png")



#### --Family level plot -- ####

colo = c("#f47f37",
         "#4ae486",
         "#a048ba",
         "#e0cc3d",
         "#42c4af",
         "#007a35",
         "#ba5392",
         "#5a7000",
         "#b40036",
         "#d3db7a",
         "#631700",
         "#ff8376",
         "#9f6a00",
         "#b05853")

Fig3 <-  Taxo_counts_sample %>% 
  subset(Substrate %in% c("Corraline","Red_crust")) %>% 
  ggplot(aes(y=Abundance,x=swab.ID,fill=Family_plot)) + 
  facet_wrap(vars(DNA.ID.updated),scale="free") + #
  geom_bar(stat="identity",position = "stack") + 
  scale_fill_manual(values=colo) +
  ylab("Relative read counts") + xlab("Samples") +
  theme_few() + 
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  guides(fill=guide_legend(title="Families"))
Fig3

ggsave(Fig3, width = 21, height = 16, units = "cm", dpi = 2400,
       filename = "./Figures/Coralline_bacteria_family_plot_updated.png")

Fig3_alt<-Fig3 +  theme(legend.position="top")+
  guides(fill=guide_legend(nrow=3,byrow=TRUE,title=NULL))

Fig3_alt

ggsave(Fig3_alt, width = 21, height = 23, units = "cm", dpi = 3000,
       filename = "./Figures/Coralline_bacteria_family_plot_updated_new.png")

#Suggested for better quality plots
#doesn't look any better than ggplot
library(Cairo)
CairoPNG("./Figures/Coralline_bacteria_family_plot_updated_new_HQ.png",
         width = 21.3, height = 23,units = "cm", dpi = 2400,)
Fig3_alt
dev.off()


#### -- Combined figure -- ####

Fig4 <- ggarrange(Fig2_col, Fig3, 
                    labels = c("A)", "B)"),
                    ncol = 1, nrow = 2,
                    heights = c(1,0.7))
Fig4

ggsave(Fig4, width = 21, height = 24, units = "cm", dpi = 600,
       filename = "./Figures/Coralline_bacteria_nMDS_and_family_plot_updated.png")



#Lets look at  mean cover for each group
sum_abund<-Taxo_counts_sample %>% 
  subset(Substrate %in% c("Corraline","Red_crust"))%>%
  group_by(DNA.ID.updated,Family_plot)%>%
  summarise(sum_abd= sum(Abundance))%>%
  ungroup()

rel_abund<-sum_abund %>% 
  group_by(DNA.ID.updated)%>%
  mutate(rel_abd = (sum_abd / sum(sum_abd)*100))%>%
  ungroup()

str(rel_abund)

# lets go to wide format
rel_abund_wide <-rel_abund%>% select(-sum_abd)  %>%
  spread(., DNA.ID.updated, rel_abd)

write.csv(rel_abund_wide,file="./microbiome/results/relative_abundance_bacteria_families_per_species.csv")
