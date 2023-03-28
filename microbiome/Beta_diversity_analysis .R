## --- Kelp & Urchin Coralline Settlement ---
## Microbiome - Beta diversity analysis
## Author: Florent Mazel 
## Edited by: Brenton Twist
## Date Created: 2021-03-05   
#-------------------------------------------------


# This script contains code for analyzing  the micriobiome community composition


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


##################
#   Load data   ##
##################

phyloseq_obj = readRDS("./microbiome/Phyloseq_object_rarefied.1000.RDS")


betaM="bray"
sample_info = as_tibble(sample_data(phyloseq_obj))


##########################
#  Supp. Mat. figure    ##
##########################


# Beta-div 
beta_div = phyloseq_obj %>% 
  phyloseq::distance(betaM)
hist(c(beta_div))


# NMDS
set.seed(258)
MDS_bray <- beta_div %>% 
  metaMDS() 

NMDS_coord = MDS_bray$points %>% as_tibble(rownames = "swab.ID") %>% 
  dplyr::left_join(sample_info)

table(sample_info$Substrate)

#PERMANOVA
model_beta = adonis2(beta_div~Substrate, data = as_tibble(sample_data(phyloseq_obj)))
model_beta2 = adonis_OmegaSq(model_beta)
model_beta2
# Df SumOfSqs      F parOmegaSq Pr(>F)    
# Substrate  5   6.2411 4.5595    0.25499  0.001 ***
# Residual  46  12.5933                             
# Total     51  18.8344   


# FIGURE 
Fig1 <- NMDS_coord %>% 
  ggplot(aes(y=MDS1,x=MDS2,shape=Substrate)) + #ylim(-.3,.3)+
  geom_point(size=2)+
  guides(shape=guide_legend(ncol=2),fill=guide_legend(ncol=2))+
  theme_few()+ 
  annotate("text",x=-0.3, y=.8, label=paste0("Stress value  = ",round(MDS_bray$stress,2))) +
theme(#legend.spacing.y = unit(0, "mm"), 
      panel.border = element_rect(colour = "black", fill=NA),
      legend.background = element_blank(),
      legend.box.background = element_rect(colour = "black"),
      legend.position = c(0.80, 0.82),legend.title=element_text(size=14),
      legend.text=element_text(size=12)
      )

Fig1

Fig1_col<- Fig1 +  geom_convexhull(aes(fill=Substrate),alpha=0.15,size=.2)+
  #scale_fill_manual(values= c("#67ae62","white","white", "#7f63b8", "#b87e39","#b94971"))
  scale_fill_manual(values= c("red","white","white", "#b0923b", "#006400","blue"))
Fig1_col
Fig1_bw<- Fig1 +  geom_convexhull(alpha=0.05,size=.2)
Fig1_bw

ggsave(Fig1_col, width = 21, height = 16, units = "cm", dpi = 600,
       filename = "./figures/Microbiome_multivariate_corallines_controls.png")

ggsave(Fig1_bw, width = 21, height = 16, units = "cm", dpi = 2400,
       filename = "./figures/Microbiome_multivariate_corallines_controls_bw.png")



##########################
#  Main      figure    ##
##########################

#### --- Statistical Analysis ---####
#-------------------------------------

# Beta-div 
phyloseq_obj_subset = phyloseq_obj %>% subset_samples(Substrate %in% c("Corraline","Red_crust"))

beta_div_sub = phyloseq_obj_subset  %>% 
  distance(betaM)
hist(c(beta_div))

# NMDS
set.seed(946)
MDS_bray2 <- beta_div_sub %>% 
  metaMDS() 

NMDS_coord2 = MDS_bray2$points %>% as_tibble(rownames = "swab.ID") %>% 
  left_join(sample_info)

#PERMANOVA
model_beta = adonis2(beta_div~Boulder+DNA.ID.updated, data = as_tibble(sample_data(phyloseq_obj_subset)),by = "margin")
model_beta2 = adonis_OmegaSq(model_beta)
model_beta2
model_beta
#Df SumOfSqs     F parOmegaSq Pr(>F)    
#Boulder        22   4.4434 1.546    0.23095  0.001 ***
#  DNA.ID.updated  8   3.4374 3.289    0.31403  0.001 ***
#  Residual        9   1.1758                            
#Total          39  12.0945  

#PERMANOVA redcrust

dataRedcrust = as_tibble(sample_data(phyloseq_obj_subset)) %>% 
  mutate(redcrust=Substrate=="Red_crust") %>% 
  group_by(redcrust) %>% sample_n(5)



model_beta = adonis2(dist_subset(beta_div,dataRedcrust$swab.ID)~redcrust, data = dataRedcrust,by = "margin")
model_beta2 = adonis_OmegaSq(model_beta)
model_beta2
model_beta
#Df SumOfSqs     F parOmegaSq Pr(>F)    
#Boulder        22   4.4434 1.546    0.23095  0.001 ***
#  DNA.ID.updated  8   3.4374 3.289    0.31403  0.001 ***
#  Residual        9   1.1758                            
#Total          39  12.0945  


#### --- Family level diversity --- #####
#----------------------------------------

Taxo_counts_Family <-  phyloseq_obj %>%
  tax_glom(taxrank = 'family' ,NArm=F) %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()

Taxonomy = as.tibble(tax_table(phyloseq_obj))
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

str(NMDS_coord2)
NMDS_coord2<- NMDS_coord2 %>% dplyr::rename(Species=DNA.ID.updated)

Fig2 <- NMDS_coord2 %>% 
  ggplot(aes(y=MDS2,x=MDS1,shape=Species)) + 
  #xlim(-.45,.85)+
  #ylim(-0.2,.5)+
  geom_point(alpha=1,size=2.5,stroke=1.1) + 
  scale_shape_manual(values=c(1,7,3,4,15,6,17,8,19))+
  theme_few() + 
  guides(shape=guide_legend(ncol=2),fill=guide_legend(ncol=2))+
  annotate("text",x=-.35, y=.4, label=paste0("Stress value  = ",round(MDS_bray2$stress,2))) +
  theme(#legend.spacing.y = unit(0, "mm"), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.position = c(0.72, 0.83),legend.title=element_text(size=13),
        legend.text=element_text(size=11)
  )
Fig2



Fig2_col <- Fig2 + geom_convexhull(aes(fill=Species), alpha=0.15,size=.2)+
  #scale_fill_manual(values= c("#F8766D", "#D39200","white", "white" , "#00C19F", "#00B9E3", "#619CFF", "#DB72FB",
                             #"#FF61C3"))
  scale_fill_manual(values= c("#F8766D", "#D39200","white", "white" , "#00C19F", "#00B9E3", "#00BA38", "#DB72FB",
                            "#FF61C3"))
Fig2_col


Fig2_bw <- Fig2 + geom_convexhull( alpha=0.07,size=.2)
Fig2_bw


ggsave(Fig2_col, width = 21, height = 16, units = "cm", dpi = 2400,
       filename = "./figures/Microbiome_multivariate_corallines.png")

ggsave(Fig2_bw, width = 21, height = 16, units = "cm", dpi = 2400,
       filename = "./figures/Microbiome_multivariate_corallines_bw.png")


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

#### -- Combined figure -- ####

Fig4 <- ggarrange(Fig2_col, Fig3, 
                    labels = c("A)", "B)"),
                    ncol = 1, nrow = 2,
                    heights = c(1,0.7))
Fig4

ggsave(Fig4, width = 21, height = 24, units = "cm", dpi = 600,
       filename = "./Figures/Coralline_bacteria_nMDS_and_family_plot_updated.png")


