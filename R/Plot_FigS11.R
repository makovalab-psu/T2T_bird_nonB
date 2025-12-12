################################################################################
### R code for plotting tandem repeat unit length histogram in zebra finch. 
# Zebra finch T2T genome. 
### written by Linnéa Smeds August 27, 2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(viridis)
library(ggh4x)
library(patchwork)
library(purrr)
plotdir="plots/"

# Input files 
file="repeats/introns_TRF.lengths.txt"
groupfile="groups.txt"

# Reading in the data 
tib<-file %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Chr=V1, RepLength=V2)
grouptib<-groupfile %>% read.table(header=TRUE) %>% as_tibble() 

# Combine 
DATA <-tib %>% full_join(grouptib) %>% 
  mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                  Group=="micro" ~ "Micro",
                                  Group=="microdot" ~ "Dot",
                                  TRUE ~ Group)) 
DATA$Classification <- factor(DATA$Classification, levels=c("Macro","Micro","Dot"))

type_col=c("#440154","#9460a1", "#e3cfe8")
################################################################################
# PLOT

p <- ggplot(DATA%>%filter(RepLength<100), aes(x = RepLength, fill = Classification, color = Classification)) + 
  geom_histogram(binwidth=1, alpha = 0.6, show.legend = FALSE) +
  facet_wrap(Classification ~., nrow=3) +  # single-column facet
  scale_fill_manual(values = type_col) +
  scale_color_manual(values = type_col) +
  labs(x="Tandem repeat unit length (bp)", y="Counts")+
  theme(
    panel.grid = element_line(color="gray95"),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    strip.text.x = element_text(size = 14),
    strip.background = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.ticks.length = unit(4, "pt"),
    #    panel.spacing.y = unit(1, "lines")
  )
p

outfile=paste(plotdir,"FigS9_TRunitLength_Introns.png", sep="")
ggsave(outfile,plot = p,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=10)

