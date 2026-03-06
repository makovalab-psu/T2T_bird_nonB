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
ABfile="ref/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.bed"

# DATA TIBBLES
dens_tib1<-densfile1 %>% read.table(header=TRUE) %>% as_tibble()
dens_tib2<-densfile2 %>% read.table(header=TRUE) %>% as_tibble()
gwdens_tib1<-gwdensfile1 %>% read.table(header=FALSE) %>% as_tibble() %>% rename(NonB=V1, Bp=V2, Density=V3)
gwdens_tib2<-gwdensfile2 %>% read.table(header=FALSE) %>% as_tibble() %>% rename(NonB=V1, Bp=V2, Density=V3)



gwdens_tib$NonB <- factor(gwdens_tib$NonB, levels=c("ALL","APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))

DATAB <- dens_tib %>%
    inner_join(typetib)
DATAB$NonB <- factor(DATAB$NonB, levels=c("ALL","APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))
DATAB$Group <- factor(DATAB$Group, levels=c("Macro", "Micro", "Dot"))

# DATA TIBBLES FOR C
centib<-cenfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  select(V1,V4,V5) %>% rename(Chr=V1, Cen_start=V4, Cen_stop=V5)

DATAC<-windfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(midpoint=Start+(Stop-Start)/2) %>% inner_join(centib) %>% 
  group_by(NonB) %>%
  mutate(label=case_when(NonB=="G4" ~ "g4discovery",
                        NonB=="G4quadron" ~ "Quadron",
                        NonB=="Z" ~ "Z-DNA Hunter m2",
                        NonB=="ZDNAm1" ~ "Z-DNA Hunter m1",
                        NonB=="Zgfa" ~ "gfa",
                        NonB=="Zseeker" ~ "ZSeeker")) %>%
  mutate(Type=if_else(NonB=="G4" | NonB=="G4quadron", "G4", "Z" ))




DATAC$label <- factor(DATAC$label, levels=c("g4discovery", "Quadron", "gfa", "Z-DNA Hunter m1", "Z-DNA Hunter m2","ZSeeker"))
ABDOT<- ABfile %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Chr=V1, Start=V2, Stop=V3,Comp=V4)


# DEFINE COLORS 
vcol <- viridis(7)
viridis_colors=c(vcol[6],vcol[7])
################################################################################
################################################################################ 
# PLOT FIGURE 1C, COVERAGE ALONG EXAMPLE CHROMOSOMES 

# Choose 1 macro, 1 micro, 1 dot from the data 
SUB <- DATAC %>% filter(Chr=="chr8_mat" | Chr=="chr36_mat" | Chr=="chr18_mat" ) %>% 
  mutate(ChrName=case_when(Chr=="chr8_mat" ~ "Chr8",
                           Chr=="chr18_mat" ~ "Chr18",
                           Chr=="chr36_mat" ~ "Chr36",
                           TRUE ~ Chr)) %>% 
  mutate(max=max(Dens)/1000)
# Get the centromeres 
CEN<-SUB %>% select(ChrName, Cen_start, Cen_stop, max, NonB) %>% distinct() %>%
  mutate(Cen_mid=(Cen_start+Cen_stop)/2)

# Make a tibble with only Chr Start and stop
CHR<-SUB %>% group_by(ChrName) %>% mutate(Start=0, Stop=max(Stop)) %>%
  select(ChrName, Start, Stop) %>% unique() %>% ungroup()

# Add AB regions 
AB<-ABDOT %>% filter(Chr=="chr8_mat"| Chr=="chr36_mat" | Chr=="chr18_mat" ) %>% 
  mutate(ChrName=case_when(Chr=="chr8_mat" ~ "Chr8",
                           Chr=="chr18_mat" ~ "Chr18",
                           Chr=="chr36_mat" ~ "Chr36",
                           TRUE ~ Chr)) %>% 
  select(Chr,Start,Stop,ChrName,Comp)

# Apply to all data frames
chr_order <- c("Chr8", "Chr18", "Chr36")
SUB$ChrName <- factor(SUB$ChrName, levels = chr_order)
CEN$ChrName <- factor(CEN$ChrName, levels = chr_order)
CHR$ChrName <- factor(CHR$ChrName, levels = chr_order)
AB$ChrName  <- factor(AB$ChrName,  levels = chr_order)

pc<-ggplot(SUB, aes(x=midpoint/1000000, y=Dens/1000, color=Type, fill=Type)) +
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
pc

# Add a dummy tract to plot the AB compartments and centromeres 
p_dummy <- ggplot() +
  geom_rect(data=CHR,
            aes(xmin=Start/1e6, xmax=Stop/1e6, ymin=0.2, ymax=1),
            fill="white", color="lightgray") +
  geom_rect(data=AB, aes(xmin=Start/1e6, xmax=Stop/1e6, ymin=0.2, ymax=1, fill=Comp), color=NA, alpha=0.6,show.legend = FALSE) +
  geom_rect(data=CEN, aes(xmin=Cen_start/1000000, xmax=Cen_stop/1e6, ymin=0.2, ymax=1), fill='red', color=NA,show.legend = FALSE) +
  facet_grid(.~ChrName, scales="free_x", space="free_x") +
  scale_fill_manual(values=c("#E38AAA","#6F9DD0"))+
  theme_void() +
  theme(
    panel.background = element_rect(fill="white", colour="white"),
    panel.spacing.x = unit(0.5, "lines"),
    strip.text.x =element_blank(),
    axis.title.x = element_text(size=14, vjust = 0.5), 
    axis.text.x  = element_text(size=10),
    axis.ticks.x = element_line(color="black", linewidth=0.3),
    axis.ticks.length.x = unit(4, "pt"),
    plot.margin  = margin(-5,20,10,30)
  )+
  scale_x_continuous(name='Position (Mb)', expand = c(0, 0))
p_dummy


leg <- ggpubr::get_legend(pc, position="top")

pcomb<-cowplot::plot_grid(pc, p_dummy, leg, ncol=1,
                          rel_heights=c(1, 0.07, 0.09))
pcomb

################################################################################
# MERGE PLOT 
outfile=paste(plotdir,"FigS23_different_annotation.png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 300, bg="white", limitsize = TRUE,width=12,height=10)



