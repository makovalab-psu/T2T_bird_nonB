################################################################################
### R code for plotting non-B motif density in Zebra finch
### written by Linnéa Smeds 30-June-2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(patchwork) 
library(ggpmisc)
library(ggpubr)
library(viridis)
library(ggh4x)
plotdir="plots/"
setwd("/Users/lbs5874/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Documents/Projects/ZebraFinch")


# INPUT FILES
densfile1="coverage/bTaeGut7v0.4_MT_rDNA.per_chrom.tsv"
gwdensfile1="coverage/bTaeGut7v0.4_MT_rDNA.per_genome.tsv"
densfile2="coverage/bTaeGut7v0.4_MT_rDNA.different.per_chrom.tsv"
gwdensfile2="coverage/bTaeGut7v0.4_MT_rDNA.different.per_genome.tsv"
windfile="coverage/bTaeGut7v0.4_MT_rDNA.different.merged.100kb.txt"
cenfile="ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff"
ABfile="ref/bTaeGut7v0.4_MT_rDNA.Cooltools.E1.200kbp.flipped.dip.collated.v0.1.bed"

# DATA TIBBLES
centib<-cenfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  select(V1,V4,V5) %>% rename(Chr=V1, Start=V4, Stop=V5) %>%
  mutate(Compartment="CEN")

DATAC<-windfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(midpoint=Start+(Stop-Start)/2) %>%
  group_by(NonB) %>%
  mutate(label=case_when(NonB=="G4" ~ "g4discovery",
                        NonB=="G4quadron" ~ "Quadron",
                        NonB=="Z" ~ "Z-DNA Hunter m2",
                        NonB=="ZDNAm1" ~ "Z-DNA Hunter m1",
                        NonB=="Zgfa" ~ "gfa",
                        NonB=="Zseeker" ~ "ZSeeker")) %>%
  mutate(Type=if_else(NonB=="G4" | NonB=="G4quadron", "G4", "Z" ))


DATAC$label <- factor(DATAC$label, levels=c("g4discovery", "Quadron", "gfa", "Z-DNA Hunter m1", "Z-DNA Hunter m2","ZSeeker"))

MERGED<- ABfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  select(V1,V2,V3,V4) %>% rename(Chr=V1, Start=V2, Stop=V3,Compartment=V4) %>%
  bind_rows(centib)


# DEFINE COLORS 
vcol <- viridis(7)
viridis_colors=c(vcol[6],vcol[7])
################################################################################
################################################################################ 
# PLOT COVERAGE ALONG EXAMPLE CHROMOSOMES 

# Choose 1 macro, 1 micro, 1 dot from the data 
SUB <- DATAC %>% filter(Chr=="chr8_mat" | Chr=="chr36_mat" | Chr=="chr18_mat" ) %>% 
  mutate(ChrName=case_when(Chr=="chr8_mat" ~ "Chr8",
                           Chr=="chr18_mat" ~ "Chr18",
                           Chr=="chr36_mat" ~ "Chr36",
                           TRUE ~ Chr)) %>% 
  mutate(max=max(Dens)/1000)

# Make a tibble with only Chr Start and stop
CHR<-SUB %>% group_by(ChrName) %>% mutate(Start=0, Stop=max(Stop)) %>%
  select(ChrName, Start, Stop) %>% unique() %>% ungroup()

# Add AB regions 
ABC<-MERGED %>% filter(Chr=="chr8_mat"| Chr=="chr36_mat" | Chr=="chr18_mat" ) %>% 
  mutate(ChrName=case_when(Chr=="chr8_mat" ~ "Chr8",
                           Chr=="chr18_mat" ~ "Chr18",
                           Chr=="chr36_mat" ~ "Chr36",
                           TRUE ~ Chr)) %>% 
  select(Chr,Start,Stop,ChrName,Compartment)

# Apply to all data frames
chr_order <- c("Chr8", "Chr18", "Chr36")
SUB$ChrName <- factor(SUB$ChrName, levels = chr_order)
CHR$ChrName <- factor(CHR$ChrName, levels = chr_order)
ABC$ChrName  <- factor(ABC$ChrName,  levels = chr_order)

p<-ggplot(SUB, aes(x=midpoint/1000000, y=Dens/1000, color=Type, fill=Type)) +
  geom_line(show.legend = FALSE) +
  geom_area(alpha = 0.7, show.legend = TRUE) +
  facet_grid(label~ChrName, scales="free", space="free_x")+
  scale_color_manual(values=viridis_colors)+
  scale_fill_manual(values=viridis_colors)+
  theme(panel.background = element_rect(fill = 'white', colour="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=18, face="bold"),
        strip.text.y = element_text(size=12, hjust=0.5),
        strip.text.x =element_text(size=12, hjust=0),
        panel.spacing.y = unit(0.2, "lines"),
        panel.spacing.x = unit(0.5, "lines"),
        plot.margin  = margin(0,0,0,0),
        axis.text.y=element_text(size=8),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title = element_text(size=12, vjust = 0.5), # Adjust 'r' to a smaller value
  )+
  scale_y_continuous(
    name = 'Coverage (%)',
    expand = c(0, 0),
    breaks = function(x) pretty(x, n = 2)  # Only ~2 ticks per panel
  )+
  scale_x_continuous(name='', expand = c(0, 0))
p

# Add a dummy tract to plot the AB compartments and centromeres 
p_dummy <- ggplot() +
  geom_rect(data=CHR,
            aes(xmin=Start/1e6, xmax=Stop/1e6, ymin=0.2, ymax=1),
            fill="white", color="lightgray") +
  geom_rect(data=ABC, aes(xmin=Start/1e6, xmax=Stop/1e6, ymin=0.2, ymax=1, fill=Compartment), color=NA,show.legend = TRUE) +
  facet_grid(.~ChrName, scales="free_x", space="free_x") +
  scale_fill_manual(values=c("#E38AAA7F","#6F9DD07F", "red"))+
  geom_point(
    data = ABC%>%filter(Compartment=="CEN"),
    aes(x = (Start + Stop)/2 / 1e6, y = 1.2),
    shape = 25,  # filled triangle pointing down
    size = 1,
    fill = "red",
    color = "red",
    inherit.aes = FALSE)+
  #theme_void() +
  theme(
    legend.position = "none",
    #legend.title = element_blank(),
    panel.background = element_rect(fill="white", colour="white"),
    panel.spacing.x = unit(0.5, "lines"),
    strip.text.x =element_blank(),
    axis.title.x = element_text(size=14, vjust = 0.5), 
    axis.text.x  = element_text(size=10),
    axis.ticks.x = element_line(color="black", linewidth=0.3),
    axis.ticks.length.x = unit(4, "pt"),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    plot.margin  = margin(-5,20,10,30)
  )+
  scale_x_continuous(name='Position (Mb)', expand = c(0, 0))+
  coord_cartesian(ylim = c(0.2, 1.1), clip = "off")
p_dummy

pcomb <- p / plot_spacer() / p_dummy +
  plot_layout(heights = c(5, -0.1, 0.15), guides = "collect") &
  theme(legend.position = "bottom")
pcomb


################################################################################
# MERGE PLOT 
outfile=paste(plotdir,"FigS2_different_annotation.png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 300, bg="white", limitsize = TRUE,width=12,height=10)



