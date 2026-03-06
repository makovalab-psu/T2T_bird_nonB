################################################################################
### R code for plotting non-B motif content vs GC content
### written by Linnéa Smeds 30-June-2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(patchwork) 
library(ggpmisc)
plotdir="plots/"


# Input files 
densfile="coverage/bTaeGut7v0.4_MT_rDNA.per_chrom.tsv"
gcfile="stats/bTaeGut7v0.4_MT_rDNA.GC_per_chrom.txt"
groupfile="helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt"

# Reading in the data 
dens_tib<-densfile %>% read.table(header=TRUE) %>% as_tibble()
gc_tib<- gcfile %>% read.table(header=TRUE) %>% as_tibble() 
grouptib <- groupfile %>% read.table(header=FALSE) %>% 
  as_tibble() %>% rename(Chr=V1, Length=V2, Group=V3) %>% 
    mutate(Group=case_when(Group=="macro" ~ "Macro",
                           Group=="micro" ~ "Micro",
                           Group=="dot" ~ "Dot",
                           TRUE ~ Group))

DATA <- dens_tib %>% inner_join(gc_tib) %>% inner_join(grouptib)
DATA$NonB <- factor(DATA$NonB, levels=c("Any","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))
DATA$Group <- factor(DATA$Group, levels=c("Macro", "Micro", "Dot"))

type_col=c("#440154","#9460a1", "#e3cfe8")
################################################################################
# PLOT Non-B CCONTENT VS GC CONTENT  

p<-ggplot(DATA, aes(x=GC_cont*100, y=Coverage*100, fill=Group, color=Group))+
  geom_point(aes(shape=Group), show.legend=TRUE, alpha=0.6, size=3)+
  labs(x="GC content (%)", y="Non-B motif content (%)", scale="free_y") +
  scale_fill_manual(values=type_col)+
  scale_color_manual(values=type_col)+
  scale_shape_manual(values=c(21,24,22))+
  facet_wrap(NonB~., ncol=3, scales = "free_y")+
  theme(panel.background = element_rect(fill = 'white', colour="black"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        panel.grid = element_blank(),
        panel.spacing.y=unit(1, "lines"),
        legend.background = element_rect(color = NA),
        legend.position="bottom",
        #legend.justification = "right",
        legend.title = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(size=18),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))+
  guides(fill = guide_legend(nrow = 1))
p

outfile=paste(plotdir,"Fig_S3_NonB_GC_ZF.png", sep="")
ggsave(outfile,plot = p,scale = 1,dpi = 300,limitsize = TRUE,width=8,height=9)


