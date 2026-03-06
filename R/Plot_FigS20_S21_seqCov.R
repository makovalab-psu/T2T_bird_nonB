################################################################################
### R code for plotting non-B motif content versus sequencing depth 
### written by Linnéa Smeds Sept 11, 2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(viridis)
library(ggh4x)
library(patchwork)
library(purrr)
library(ggpubr)
plotdir="plots/"


# Input files 
hififile="coverage/bTaeGut7v0.4_MT_rDNA.hifi.nonB_and_seqCov.bed"
ontfile="coverage/bTaeGut7v0.4_MT_rDNA.ont.nonB_and_seqCov.bed"
groupfile="helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt"

# Reading in the data 
grouptib<-groupfile %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Chr=V1, Length=V2, Group=V3)

DATAH<-hififile %>% read.table(header=TRUE) %>% as_tibble() %>% inner_join(grouptib) %>% 
  filter(Group=="macro" | Group=="micro") %>%  
  mutate(Group=case_when(Group=="macro" ~ "Macro",
                         Group=="micro" ~ "Micro",
                         TRUE ~ Group)) 


DATAO<-ontfile %>% read.table(header=TRUE) %>% as_tibble() %>% inner_join(grouptib) %>% 
  filter(Group=="macro" | Group=="micro") %>% 
  mutate(Group=case_when(Group=="macro" ~ "Macro",
                         Group=="micro" ~ "Micro",
                         TRUE ~ Group)) 



DATAH$Group <- factor(DATAH$Group, levels=c("Macro","Micro","Dot"))
DATAH$NonB <- factor(DATAH$NonB, levels=c("Any","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))
DATAO$Group <- factor(DATAO$Group, levels=c("Macro","Micro","Dot"))
DATAO$NonB <- factor(DATAO$NonB, levels=c("Any","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))


# Set colors
vcolors=viridis(7)
type_col=c("Macro"="#440154","Micro"="#9460a1", "Dot"="#e3cfe8")

################################################################################
# Plot 
# Make function to use for each seq type separately 

make_plot <- function(data, seqtype = NULL) {
  p <- ggplot(data, aes(x = SeqCov, y = NonBCov*100)) +
    geom_point(color="lightgray", alpha = 0.1, size = 1, show.legend = FALSE) +
    stat_smooth(method = "lm", se = TRUE, level=0.95, linetype = "dashed", color="darkred", fill="red", show.legend=FALSE) + # Customize method, standard error, color, and linetype
    facet_grid(rows=vars(Group), cols=vars(NonB)) +
    scale_fill_manual(values = type_col) +
    scale_color_manual(values = type_col) +
    stat_cor(aes(label=after_stat(rr.label)), color="red3", # Display only the R-squared label
             label.x = 95, label.y = 97, hjust=1, show.legend=FALSE)+ # Adjust position as needed
    stat_cor(method="spearman", aes(label = after_stat(sub("R", "ρ", r.label))), color="blue", # Display only the rho label
             label.x = 95, label.y = 87, hjust=1, show.legend=FALSE)+ # Adjust position as needed
    
    labs(x=paste(seqtype, "coverage (X)", sep=" "), y="Non-B motif coverage (%)")+
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      strip.text.x = element_text(size = 12),
      strip.background = element_blank(),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 12),
      axis.ticks.length = unit(4, "pt"),
      panel.spacing = unit(0.4, "lines")
    )+
    scale_x_continuous(limits=c(0,100), expand=c(00.1,0.01))+
    scale_y_continuous(limits=c(0,101), expand=c(0.01,0.01))
  p
  return(p)
}

p1<-make_plot(DATAH, "Hifi")
p1
outfile=paste(plotdir,"FigS20_nonB_vs_coverage.hifi.macromicro.png", sep="")
ggsave(outfile,plot = p1,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=6)


p2<-make_plot(DATAO, "ONT")
p2
outfile=paste(plotdir,"FigS21_nonB_vs_coverage.ONT.macromicro.png", sep="")
ggsave(outfile,plot = p2,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=6)



