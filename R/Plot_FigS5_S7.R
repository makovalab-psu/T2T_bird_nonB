################################################################################
### R code for plotting non-B motif overlaps in Chicken
### written by Linnéa Smeds 11-July-2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(ggupset)
require(ggtext)
require(patchwork)
require(viridis)
library(dplyr)
library(purrr)
options(scipen = 999, decimals=2)
plotdir="plots/"

# Reading in the data, create tibbles 
fileMulti="overlap/chicken.v23.merged.summary.txt"
filePair="overlap/chicken.v23.merged.pairwise.txt"

# Colorscheme
viridis_colors <- viridis(5)

tibM<-fileMulti %>% read.table(header=TRUE) %>% as_tibble() %>% 
  rename(Comb=NonB, Bp=Overlap) %>% mutate(Types=strsplit(Comb, "-")) %>% 
  mutate(PrintName=case_when(Region=="autosomes" ~ "Autosomes",
                             Region=="chrW" ~ "Chromosome W",
                             Region=="chrZ" ~ "Chromosome Z",
                             Region=="macro" ~ "Macro",
                             Region=="micro" ~ "Micro",
                             Region=="microdot" ~ "Dot")) 
tibM$PrintName=factor(tibM$PrintName, levels=c("Macro", "Micro", "Dot", "Autosomes", "Chromosome Z", "Chromosome W"))

tibP<-filePair %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(textcol=if_else(Frac>0.5,"white","black")) %>% 
  mutate(PrintName=case_when(Region=="autosomes" ~ "Autosomes",
                       Region=="chrW" ~ "Chromosome W",
                       Region=="chrZ" ~ "Chromosome Z",
                       Region=="macro" ~ "Macro",
                       Region=="micro" ~ "Micro",
                       Region=="microdot" ~ "Dot"))
tibP$NB1=factor(tibP$NB1, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))
tibP$NB2=factor(tibP$NB2, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))
tibP$PrintName=factor(tibP$PrintName, levels=c("Macro", "Micro", "Dot", "Autosomes", "Chromosome Z", "Chromosome W"))


################################################################################
# PART1 upset plot
# Only plot groups with more than 10kb intersect 

SUBM <- tibM %>% filter(Bp>10000) %>% 
  filter(str_detect(Region, "cro"))%>% arrange(desc(Bp)) %>% 
#  filter(!str_detect(Region, "cro"))%>% arrange(desc(Bp)) %>% 
  mutate(vcolor=viridis_colors[lengths(Types)])

SUBM <- SUBM %>%
  mutate(n_types = map_int(Types, length))


overlap_levels <- sort(unique(SUBM$n_types))
viridis_colors <- viridis(length(overlap_levels))

# Get the number of dots for each color 
counts<-SUBM %>% select(Comb, n_types) %>% group_by(Comb, n_types) %>%
  unique() %>% ungroup() %>% group_by(n_types) %>% select(n_types) %>%
  summarize(n=n()) %>% select(n) %>%deframe()

repeats <- counts * seq_along(counts)
color_vector <- rep(viridis_colors, times = repeats)

 
pup<-ggplot(SUBM, aes(x=Types)) +
  geom_col(aes(y=Bp/1000000, fill=vcolor), alpha=0.8, show.legend = FALSE) +
  facet_wrap(~PrintName, nrow=6, scales="free_y", strip.position="right")+
  scale_fill_identity() + #manual(values=viridis_colors)+
#  scale_color_identity() + #manual(values=viridis_colors)+
  scale_x_upset(order_by="degree")+
  theme_combmatrix(
   combmatrix.panel.point.color.fill = color_vector,
   combmatrix.panel.line.color = color_vector,
   combmatrix.panel.point.color.empty = "gray90",
   combmatrix.panel.striped_background.color.two = "grey95",
   combmatrix.panel.point.size = 2,
  ) +
  labs(x="") +
  theme(panel.background = element_rect(fill = 'white', colour="gray"),
   #     panel.grid.major = element_line(colour = 'gray95'),
        strip.background = element_rect(fill = 'white', colour="white"),
        plot.title = element_text(size=18),
        strip.text = element_markdown(size=14, face="bold"),
      #  panel.spacing = unit(1, "lines"),
        axis.text.y=element_text(size=12),
        axis.title.y = element_text(size=14,vjust = 0.5), # Adjust 'r' to a smaller value
        axis.ticks.x=element_blank())+
  scale_y_continuous(name='Mb')+
  labs(title = bold('B') ~ '')
pup

################################################################################
# PART2 heatmap with pairwise overlap

SUBP<-tibP %>% filter(str_detect(Region, "cro"))
#SUBP<-tibP %>% filter(!str_detect(Region, "cro"))

ppw<-ggplot(SUBP, aes(x=NB2, y=forcats::fct_rev(NB1), fill=Frac)) +
  geom_tile(color="white", show.legend = FALSE) +
  geom_text(aes(label=round(Frac*100, digits=1), color=textcol), size=4, show.legend=FALSE) +
  facet_wrap(~PrintName, nrow=6, scales='free', strip.position="right") +
  labs(x="", y="") +
  scale_fill_gradient(low = "#f2f6ff", high = viridis_colors[2], name="Overlap", na.value = 'white') +
  scale_color_manual(values=c("black","white")) +
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.grid.minor = element_line(colour = 'white'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.spacing = unit(0, "lines"),
        legend.position="bottom",
        plot.title = element_text(size=18),
        legend.title=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.text=element_text(size=12),
        strip.background = element_rect(fill='white'),
        strip.text = element_text(size=14, face="bold"))+
  labs(title = bold('A') ~ '')
ppw


################################################################################
# COMBINE Upsetplot with Heatmap

pcomb <- ppw + pup + plot_layout(widths = c(1, 2.5)) 
pcomb

outfile=paste(plotdir,"Chicken_Overlap_MacroMicroDot_viridis.png", sep="")
#outfile=paste(plotdir,"Chicken_Overlap_AutoZW_viridis.png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 600,limitsize = TRUE,width=12,height=9)

# Make separate panels to paste in illustrator
#outfile=paste(plotdir,"Overlap_MacroMicroDot_A.png", sep="")
#outfile=paste(plotdir,"Overlap_AutoZW_A.png", sep="")
#ggsave(outfile,plot = ppw,scale = 1,dpi = 600,limitsize = TRUE,width=4,height=9)

outfile=paste(plotdir,"Overlap_MacroMicroDot_B.png", sep="")
#outfile=paste(plotdir,"Overlap_AutoZW_B.png", sep="")
ggsave(outfile,plot = pup,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=9)



