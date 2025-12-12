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
plotdir="plots/"
setwd("/Users/lbs5874/Documents/Projects/ZebraFinch/")

# FIGURE 1A ARE NON-B SYMBOLS, ADDED IN ADOBE ILLUSTRATOR 

# INPUT FILES FOR B
densfile="densities/bTaeGut7v0.4_MT_rDNA.nonB_per_chrom.txt"
gwdensfile="densities/bTaeGut7v0.4_MT_rDNA.nonB_genome_wide.txt"
macrfile="macro.txt"
micrfile="micro.txt"
dotfile="microdot.txt"

# INPUT FILES FOR C
windfile="densities/bTaeGut7v0.4_MT_rDNA.merged.100kb.txt"
cenfile="ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff"
ABfile="ref/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.bed"

# DATA TIBBLES FOR B 
dens_tib<-densfile %>% read.table(header=TRUE) %>% as_tibble()
gwdens_tib<-gwdensfile %>% read.table(header=FALSE) %>% as_tibble() %>% rename(NonB=V1, Bp=V2, Density=V3)
gwdens_tib$NonB <- factor(gwdens_tib$NonB, levels=c("ALL","APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))
macrtib<-macrfile %>% read.table(header=FALSE) %>% 
  as_tibble() %>% rename(Chr=V1) %>% mutate(Group="Macro")
micrtib<-micrfile %>% read.table(header=FALSE) %>% 
  as_tibble() %>% rename(Chr=V1) %>% mutate(Group="Micro")
dottib<- dotfile %>% read.table(header=FALSE) %>% 
  as_tibble() %>% rename(Chr=V1) %>% mutate(Group="Dot")
typetib<- macrtib %>% bind_rows(micrtib) %>% bind_rows(dottib)

DATAB <- dens_tib %>%
    inner_join(typetib)
DATAB$NonB <- factor(DATAB$NonB, levels=c("ALL","APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))
DATAB$Group <- factor(DATAB$Group, levels=c("Macro", "Micro", "Dot"))

# DATA TIBBLES FOR C
centib<-cenfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  select(V1,V4,V5) %>% rename(Chr=V1, Cen_start=V4, Cen_stop=V5)

DATAC<-windfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(midpoint=Start+(Stop-Start)/2) %>% inner_join(centib) %>% 
  group_by(NonB)
DATAC$NonB <- factor(DATAC$NonB, levels=c("APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))
ABDOT<- ABfile %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Chr=V1, Start=V2, Stop=V3,Comp=V4)


# DEFINE COLORS 
type_col=c("#440154","#9460a1", "#e3cfe8")
viridis_colors <- viridis(8)
################################################################################
# PLOT FIGURE 1B, COVERAGE IN DIFFERENT GROUPS 

# Define pariwise comparisons for stats test 
comparisons <- list(
  c("Macro", "Micro"),
  c("Micro", "Dot"),
  c("Macro", "Dot")
)

pb<-ggplot(DATAB, aes(x=Group, y=Density*100, fill=Group, color=Group))+
  geom_hline(data = gwdens_tib, aes(yintercept = Density*100), linetype="dashed", color = "red3")+
  geom_boxplot(show.legend=TRUE, alpha=0.9)+
  facet_wrap(vars(NonB), ncol=9, scales="free_y")+
  scale_fill_manual(values=type_col)+
  scale_color_manual(values=type_col)+
  stat_compare_means(
    comparisons = comparisons, label = "p.signif", p.adjust.method = "BH", # FDR correction
    label.y.npc = "top", tip.length = 0.005,
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    )
  ) +
  labs(x="", y="Coverage (%)", scale="free_y") +
  theme(panel.background = element_rect(fill = 'white', colour="black"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        panel.grid = element_blank(),
        legend.margin = margin(0,0,0,0),
        plot.margin = margin(1,1,0,1),
        legend.background = element_blank(),
        legend.text=element_text(size=14),
        legend.position="bottom",
        axis.ticks.x=element_blank(),
        legend.justification = "center",
        legend.key=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size=14),
        plot.title=element_text(size=18, face="bold"),
  )+
  guides(fill = guide_legend(nrow = 1))+
  ggtitle("B")
pb


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

pc<-ggplot(SUB, aes(x=midpoint/1000000, y=Dens/1000, color=NonB, fill=NonB)) +
  #  geom_rect(data=CEN, aes(xmin=Cen_start/1000000, ymin=0.0, ymax=max+0.01, xmax=Cen_stop/1000000), inherit.aes = FALSE, fill="gray", color="gray", linewidth=1, show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  geom_area(alpha = 0.7, show.legend = FALSE) +
  facet_grid(NonB~ChrName, scales="free", space="free_x")+
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
        axis.title = element_text(size=12, vjust = 0.5), # Adjust 'r' to a smaller value
  )+
  scale_y_continuous(
    name = 'Coverage (%)',
    expand = c(0, 0),
    breaks = function(x) pretty(x, n = 2)  # Only ~2 ticks per panel
  )+
  scale_x_continuous(name='', expand = c(0, 0))+
  ggtitle("C")
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
    plot.margin  = margin(5,0,-10,0)
  )+
  scale_x_continuous(name='Position (Mb)', expand = c(0, 0))
p_dummy

p_comb<- pc/ plot_spacer() / p_dummy  + plot_layout(heights = c(5,-0.6, 0.2), guides = "collect") 
p_comb

################################################################################
# MERGE PLOT B AND C 

p_final<- (pb / p_comb)+plot_layout(heights=c(1,1))
p_final

outfile=paste(plotdir,"Fig1BC_Density.png", sep="")
outfile=paste(plotdir,"Fig1BC_Density.svg", sep="")
ggsave(outfile,plot = p_final,scale = 1,dpi = 600,limitsize = TRUE,width=12,height=9)



