################################################################################
### R code for plotting enrichment in the AB compartments of zebra finch micro 
### and macro chromosomes 
### written by Linnéa Smeds 2026

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(patchwork) 
library(ggpmisc)
library(ggpubr)
library(dplyr)
plotdir="plots/"


compfile1="compart/macro.summary.200kb.txt"
compfile2="compart/micro.summary.200kb.txt"
compfile3="compart/dot_summary.200kb.txt"

# Reading in the data 
DATA1<-compfile1 %>% read.table(header=TRUE) %>% 
  as_tibble() %>% mutate(NonB=if_else(NonB=="Any", "All", NonB))
DATA1$NonB <- factor(DATA1$NonB, levels=c("All","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))

DATA2<-compfile2 %>% read.table(header=TRUE) %>% 
  as_tibble() %>% mutate(NonB=if_else(NonB=="Any", "All", NonB))
DATA2$NonB <- factor(DATA2$NonB, levels=c("All","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))

DATA3<-compfile3 %>% read.table(header=TRUE) %>% 
  as_tibble() %>% mutate(NonB=if_else(NonB=="Any", "All", NonB))
DATA3$NonB <- factor(DATA3$NonB, levels=c("All","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))


# Set colors
comp_col=c("A"="#E38AAA","B"="#6F9DD0", "."="lightgray")
################################################################################
# PLOT FUNCTION FOR NON-B MOTIF ENRICHMENT IN A vs B COMPARTMENTS 

make_ab_violinplot<-function(data, group)  {
  # Test differences between boxes with wilcoxon 
  pvals <- data %>%
    group_by(NonB) %>%
    summarise(
      p = wilcox.test(Density ~ Compartment)$p.value
    ) %>%
    mutate(p_adj = p.adjust(p, method = "fdr"))
  
  SIGN <- pvals %>%
    mutate(y.position =  max(data$Density)*102,  # adjust height
           group1 = "A",
           group2 = "B",
           label = case_when(
             p_adj <= 0.001 ~ "***",
             p_adj <= 0.01  ~ "**",
             p_adj <= 0.05  ~ "*",
             TRUE           ~ "ns"))
  
  pA<-ggplot(data, aes(x=Compartment, y=Density*100, fill=Compartment, color=Compartment))+
 #   geom_hline(yintercept=1, linetype="dashed", color = "red3")+
    geom_line(aes(group = Chr), color="lightgray", linewidth=0.2) + 
    geom_violin(show.legend=FALSE, alpha=0.7)+
    geom_point(size = 0.5, show.legend=FALSE) + 
    facet_wrap(NonB~., nrow=1)+
    scale_fill_manual(values=comp_col)+
    scale_color_manual(values=comp_col)+
    labs(x="Compartment", y="%Non-B DNA motif coverage") +
    stat_pvalue_manual(SIGN, label = "label",
                       tip.length = 0.01, bracket.size = 0.5,
                       inherit.aes = FALSE)+
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = 'white', colour="black"),
          strip.background = element_rect(fill = 'white', colour="white"),
          strip.text = element_text(size=12),
          plot.title = element_text(size=18, face="bold"),
          axis.title = element_text(size=12),
          axis.text = element_text(size=11))+
    ggtitle(group)
  return(pA)
}

### Run function 
p1<-make_ab_violinplot(DATA1, "A")
p2<-make_ab_violinplot(DATA2, "B")
p3<-make_ab_violinplot(DATA3, "C")


pcomb <- p1 / p2 /p3
pcomb

outfile=paste(plotdir,"Fig23_nonB_in_AB_chromCategory.png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 300,limitsize = TRUE,width=9,height=12)


