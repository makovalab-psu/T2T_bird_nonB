################################################################################
### R code for plotting enrichment in the AB compartments of zebra finch dot 
### chromosomes and long-read sequenceing coverage in relation to non-B motifs. 
### written by Linnéa Smeds 17-July-2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(patchwork) 
library(ggpmisc)
library(ggpubr)
library(dplyr)
plotdir="plots/"

# Input files 
compfile="compart/dot_summary.10kb.txt"
hifiABfile="coverage/bTaeGut7v0.4_MT_rDNA.hifi.nonB_and_seqCov.dotComp.tsv"   
ontABfile="coverage/bTaeGut7v0.4_MT_rDNA.ont.nonB_and_seqCov.dotComp.tsv"   

# Reading in the data 
DATAA<-compfile %>% read.table(header=TRUE) %>% as_tibble()
DATAA$NonB <- factor(DATAA$NonB, levels=c("Any","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))

DATAB <- hifiABfile %>% read.table(header=TRUE) %>% as_tibble()
DATAB$Comp <- factor(DATAB$Comp, levels=c("A","B","."))
DATAB$NonB <- factor(DATAB$NonB, levels=c("Any","APR", "DR", "STR", "IR","TRI", "G4", "Z"))

DATAC <- ontABfile %>% read.table(header=TRUE) %>% as_tibble()
DATAC$Comp <- factor(DATAC$Comp, levels=c("A","B","."))
DATAC$NonB <- factor(DATAC$NonB, levels=c("Any","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))


# Set colors
comp_col=c("A"="#E38AAA","B"="#6F9DD0", "."="lightgray")
################################################################################
# PLOT A) NON-B MOTIF ENRICHMENT IN A vs B COMPARTMENTS ON DOT CHROMOSOMES 

# Test differences between boxes with wilcoxon 
pvals <- DATAA %>%
  group_by(NonB) %>%
  summarise(
    p = wilcox.test(Enrichment ~ Compartment)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p, method = "fdr"))

SIGN <- pvals %>%
  mutate(y.position =  max(DATAA$Enrichment) + 0.5,  # adjust height
         group1 = "A",
         group2 = "B",
         label = case_when(
           p_adj <= 0.001 ~ "***",
           p_adj <= 0.01  ~ "**",
           p_adj <= 0.05  ~ "*",
           TRUE           ~ "ns"))

pA<-ggplot(DATAA, aes(x=Compartment, y=Enrichment, fill=Compartment, color=Compartment))+
  geom_hline(yintercept=1, linetype="dashed", color = "red3")+
  geom_line(aes(group = Chr), color="lightgray", linewidth=0.2) + 
  geom_violin(show.legend=FALSE, alpha=0.7)+
  geom_point(size = 0.5, show.legend=FALSE) + 
  facet_wrap(NonB~., nrow=1)+
  scale_fill_manual(values=comp_col)+
  scale_color_manual(values=comp_col)+
  labs(x="Compartment", y="Fold enrichment") +
  stat_pvalue_manual(SIGN, label = "label",
                     tip.length = 0.01, bracket.size = 0.5,
                     inherit.aes = FALSE)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour="black"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        legend = element_blank(),
        plot.title = element_text(size=18, face="bold"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11))+
  ggtitle("A")
pA


################################################################################
# PLOT B) NON-B MOTIF DENSITY VS HIFI COVERAGE 

pB <- ggplot(DATAB, aes(x = SeqCov, y = NonBCov*100/1024)) +
  geom_point(aes(fill = Comp, color=Comp), shape=16, alpha = 0.1, size = 1, show.legend = FALSE) +
  stat_smooth(method = "lm", se = TRUE, level=0.95, linewidth=0.5, linetype = "dashed", color="red3", fill="red", show.legend=FALSE) + # Customize method, standard error, color, and linetype
  facet_wrap(NonB~., nrow=1) +  # single-column facet
  scale_fill_manual(values = comp_col) +
  scale_color_manual(values = comp_col) +
  stat_cor(aes(label=after_stat(rr.label)), color="red3", # Display only the R-squared label
           label.x = 95, label.y = 97, hjust=1, show.legend=FALSE)+ # Adjust position as needed
  stat_cor(method="spearman", aes(label = after_stat(sub("R", "ρ", r.label))), color="blue", # Display only the rho label
           label.x = 95, label.y = 87, hjust=1, show.legend=FALSE)+ # Adjust position as needed
  labs(x="HiFi coverage (X)", y="Non-B motif coverage (%)")+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    strip.text.x = element_text(size = 12),
    strip.background = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size=18, face="bold"),
  )+
  scale_x_continuous(limits=c(0,100), expand=c(0.03,0.03))+
  scale_y_continuous(limits=c(0,101), expand=c(0.01,0.01))+
  ggtitle("B")
pB

################################################################################
# PLOT C) NON-B MOTIF DENSITY VS ONT COVERAGE 

pC <- ggplot(DATAC, aes(x = SeqCov, y = NonBCov*100/1024)) +
  geom_point(aes(fill = Comp, color=Comp), shape=16, alpha = 0.1, size = 1, show.legend = FALSE) +
  stat_smooth(method = "lm", se = TRUE, level=0.95, linewidth=0.5, linetype = "dashed", color="red3", fill="red", show.legend=FALSE) + # Customize method, standard error, color, and linetype
  facet_wrap(NonB~., nrow=1) +  # single-column facet
  scale_fill_manual(values = comp_col) +
  scale_color_manual(values = comp_col) +
  stat_cor(aes(label=after_stat(rr.label)), color="red3", # Display only the R-squared label
           label.x = 95, label.y = 97, hjust=1, show.legend=FALSE)+ # Adjust position as needed
  stat_cor(method="spearman", aes(label = after_stat(sub("R", "ρ", r.label))), color="blue", # Display only the rho label
           label.x = 95, label.y = 87, hjust=1, show.legend=FALSE)+ 
  labs(x="ONT coverage (X)", y="Non-B motif coverage (%)")+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    strip.text.x = element_text(size = 12),
    strip.background = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size=18, face="bold"),
  )+
  scale_x_continuous(limits=c(0,100), expand=c(0.03,0.03))+
  scale_y_continuous(limits=c(0,101), expand=c(0.01,0.01))+
  ggtitle("C")
pC

pcomb<- pA + pB + pC + plot_layout(heights=c(1,1,1))
pcomb

outfile=paste(plotdir,"Fig6_nonB_in_ABdots.png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=9)
outfile=paste(plotdir,"Fig6_nonB_in_ABdots.svg", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 600,limitsize = TRUE,width=10,height=9)







################################################################################
# SAME AS ABOVE BUT ALL DOTS IN GRAY

DATAB<-DATAB %>% filter(NonB=="ALL")
DATAB$NonB <- factor(DATAB$NonB, levels=c("ALL"))

DATAC<-DATAC %>% filter(NonB=="ALL")
DATAC$NonB <- factor(DATAC$NonB, levels=c("ALL"))

pB <- ggplot(DATAB, aes(x = SeqCov, y = NonBCov*100/1024)) +
  geom_point(fill = "lightgrey", color="lightgrey", shape=16, alpha = 0.5, size = 1, show.legend = FALSE) +
  stat_smooth(method = "lm", se = TRUE, level=0.95, linewidth=0.5, linetype = "dashed", color="red3", fill="red", show.legend=FALSE) + # Customize method, standard error, color, and linetype
  facet_wrap(NonB~., nrow=1) +  # single-column facet
  stat_cor(aes(label=after_stat(rr.label)), color="red3", # Display only the R-squared label
           label.x = 95, label.y = 97, hjust=1, show.legend=FALSE)+ # Adjust position as needed
  stat_cor(method="spearman", aes(label = after_stat(sub("R", "ρ", r.label))), color="blue", # Display only the rho label
           label.x = 95, label.y = 87, hjust=1, show.legend=FALSE)+ # Adjust position as needed
  labs(x="HiFi coverage (X)", y="Non-B motif coverage (%)")+
  theme(
    panel.grid = element_blank(),
    strip.text = element_blank(),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    strip.background = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size=18, face="bold"),
  )+
  scale_x_continuous(limits=c(0,100), expand=c(0.03,0.03))+
  scale_y_continuous(limits=c(0,101), expand=c(0.01,0.01))+
  ggtitle("PacBio")
pB

pC <- ggplot(DATAC, aes(x = SeqCov, y = NonBCov*100/1024)) +
  geom_point(fill = "lightgrey", color="lightgrey", shape=16, alpha = 0.5, size = 1, show.legend = FALSE) +
  stat_smooth(method = "lm", se = TRUE, level=0.95, linewidth=0.5, linetype = "dashed", color="red3", fill="red", show.legend=FALSE) + # Customize method, standard error, color, and linetype
  facet_wrap(NonB~., nrow=1) +  # single-column facet
  stat_cor(aes(label=after_stat(rr.label)), color="red3", # Display only the R-squared label
           label.x = 95, label.y = 97, hjust=1, show.legend=FALSE)+ # Adjust position as needed
  stat_cor(method="spearman", aes(label = after_stat(sub("R", "ρ", r.label))), color="blue", # Display only the rho label
           label.x = 95, label.y = 87, hjust=1, show.legend=FALSE)+ 
  labs(x="ONT coverage (X)", y="Non-B motif coverage (%)")+
  theme(
    panel.grid = element_blank(),
    strip.text = element_blank(),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    strip.background = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size=18, face="bold"),
  )+
  scale_x_continuous(limits=c(0,100), expand=c(0.03,0.03))+
  scale_y_continuous(limits=c(0,101), expand=c(0.01,0.01))+
  ggtitle("Oxford Nanopore")
pC

outfile=paste(plotdir,"Figure_for_presentation_nonB_vs_PacBio.png", sep="")
ggsave(outfile,plot = pB,scale = 1,dpi = 300,limitsize = TRUE,width=4,height=4)
outfile=paste(plotdir,"Figure_for_presentation_nonB_vs_ONT.png", sep="")
ggsave(outfile,plot = pC,scale = 1,dpi = 300,limitsize = TRUE,width=4,height=4)
