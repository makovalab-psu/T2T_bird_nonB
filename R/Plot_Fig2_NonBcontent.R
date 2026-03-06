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
setwd("/Users/lbs5874/OneDrive - The Pennsylvania State University/Documents/Projects/ZebraFinch/")

# READ IN SPECIES LIST
spfile="helpfiles/species_list.txt"
sptib<-spfile %>% read.table(header=FALSE) %>% as_tibble()

# DEFINE COLORS 
type_col=c("#440154","#9460a1", "#e3cfe8")
viridis_colors <- viridis(7)
################################################################################

# DEFINE COLORS 
type_col=c("#440154","#9460a1", "#e3cfe8", "lightgray")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function for plotting 

# ONE FACET FOR EACH NONB TYPE 
make_boxplot <- function(data1, data2) {
  
  # Define pairwise comparisons for stats test
  groups <- as.character(unique(droplevels(data1$Group)))
  comparisons <- combn(groups, 2, simplify = FALSE)
  
  p<-ggplot(data1, aes(x=Group, y=Coverage*100, fill=Group, color=Group))+
    geom_hline(data = data2, aes(yintercept = Coverage*100), linetype="dashed", color = "red3")+
    geom_boxplot(show.legend=TRUE, alpha=0.9)+
    facet_wrap(vars(NonB), ncol=9, scales="free_y")+
    scale_fill_manual(values=type_col)+
    scale_color_manual(values=type_col)+
    stat_compare_means(
      comparisons = comparisons, label = "p.signif", p.adjust.method = "BH", # FDR correction
      label.y.npc = "top", tip.length = 0.005,
      bracket.nudge.y = -1,
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
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
    guides(fill = guide_legend(nrow = 1))+
    ggtitle("")
  return(p)
}

# ONE FACET FOR EACH SPECIES 
make_spboxplot <- function(data1, data2) {
  
  # Define pairwise comparisons for stats test
  groups <- unique(as.character(data1$Group))
  comparisons <- combn(groups, 2, simplify = FALSE)
  
  p<-ggplot(data1, aes(x=Group, y=Coverage*100, fill=Group, color=Group))+
    geom_hline(data = data2, aes(yintercept = Coverage*100), linetype="dashed", color = "red3")+
    geom_boxplot(show.legend=TRUE, alpha=0.9)+
    facet_wrap(~Species, ncol=8)+
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
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.17)))+
    guides(fill = guide_legend(nrow = 1))+
    ggtitle("")
  return(p)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PREP BOXPLOTS 

# Make empty list 
plots <- vector("list", length(sptib$V3))
names(plots) <- sptib$V3
COMB <-tibble(Species=character(0), Group=character(0), NonB=character(0), Coverage=numeric(0))
GWCOMB <-tibble(Species=character(0), NonB=character(0), Coverage=numeric(0))

# Loop over files in species 
for (i in 1:length(sptib$V3)){
  show(paste("Running pipeline for:",sptib$V3[i]))
  covfile=paste("coverage/",sptib$V3[i],".per_chrom.tsv", sep="")
  gwcovfile=paste("coverage/",sptib$V3[i],".per_genome.tsv", sep="")
  grfile=paste("helpfiles/",sptib$V3[i],".groups.txt", sep="")
  
  # Make tibbles 
  covtib<-covfile %>% read.table(header=TRUE) %>% as_tibble()
  gwcovtib<-gwcovfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
    rename(NonB=V1, Bp=V2, Coverage=V3) 
  gwcovtib$NonB <- factor(gwcovtib$NonB, levels=c("Any","APR", "DR", "STR", "IR","TRI", "G4", "Z"))
  gwcovtib <- gwcovtib %>% drop_na()
  typetib <- grfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
    rename(Chr=V1, Length=V2, Group=V3) %>%
    mutate(Group=case_when(Group=="macro" ~ "Macro",
                           Group=="micro" ~ "Micro",
                           Group=="dot" ~ "Dot",
                           Group=="unplaced" ~ "Unplaced",
                           TRUE ~ Group))
  
  # Joined tibble
  DATA <- covtib %>%
    inner_join(typetib)
  DATA$NonB <- factor(DATA$NonB, levels=c("Any","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))
  DATA$Group <- factor(DATA$Group, levels=unique(DATA$Group))
  DATA<-DATA %>% drop_na() %>% mutate(Species=sptib$V1[i])
  COMB<-bind_rows(COMB,DATA%>%select(-Bp,-Length))
  GWCOMB<-bind_rows(GWCOMB,gwcovtib%>%select(-Bp)%>%mutate(Species=sptib$V1[i]))
  
  p<-make_boxplot(DATA, gwcovtib)
  p
  plots[[i]] <- p
}


for (nm in names(plots)) {
  ggsave(
    filename = paste0("plots/", nm, "_NonB_boxplot.png"),
    plot = plots[[nm]],
    scale = 1,dpi = 300,limitsize = TRUE,width=12,height=4
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIGURE 2A, TAKE FROM ABOVE 

pa<-plots[[1]] +ggtitle("A")
pa


#( ALL SPECIES IN ONE, ONLY USING TOTAL Non-B CONTENT )
ALL <- COMB%>%filter(NonB=="Any")%>%
  mutate(Species=factor(Species, levels=c("zebra_finch",
                                          "ural_owl",
                                          "bandtailed_pigeon",
                                          "annas_hummingbird",
                                          "great_bustard",
                                          "chicken",
                                          "peking_duck",
                                          "emu")))
ALL$Group <- factor(ALL$Group, levels=unique(ALL$Group))

ALLGW <- GWCOMB%>%filter(NonB=="Any") %>% 
  mutate(Species=factor(Species, levels=c("zebra_finch",
                                          "ural_owl",
                                          "bandtailed_pigeon",
                                          "annas_hummingbird",
                                          "great_bustard",
                                          "chicken",
                                          "peking_duck",
                                          "emu")))


pall<-make_spboxplot(ALL, ALLGW)
ggsave(filename = "plots/8sp_NonB_boxplot.png",
       plot = pall,
       scale = 1,dpi = 300,limitsize = TRUE,width=12,height=5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT FIGURE 2C (remove ZF from above)

OTHER<-ALL %>% filter(Species!="zebra_finch")
OTHERGW<-ALLGW %>% filter(Species!="zebra_finch")
pc<-make_spboxplot(OTHER, OTHERGW)
pc=pc+ggtitle("C")+theme(strip.text=element_blank())
pc             



################################################################################ 
# PLOT FIGURE 2B, COVERAGE ALONG EXAMPLE CHROMOSOMES 

# INPUT FILES FOR B
windfile="coverage/bTaeGut7v0.4_MT_rDNA.merged.100kb.txt"
cenfile="ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff"
ABfile="ref/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.bed"

# DATA TIBBLES FOR B
centib<-cenfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  select(V1,V4,V5) %>% rename(Chr=V1, Cen_start=V4, Cen_stop=V5)

DATA<-windfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(midpoint=Start+(Stop-Start)/2) %>% inner_join(centib) %>% 
  group_by(NonB)
DATA$NonB <- factor(DATA$NonB, levels=c("APR", "DR", "STR", "IR", "TRI", "G4", "Z"))
ABDOT<- ABfile %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Chr=V1, Start=V2, Stop=V3,Comp=V4)



# Choose 1 macro, 1 micro, 1 dot from the data 
SUB <- DATA %>% filter(Chr=="chr8_mat" | Chr=="chr36_mat" | Chr=="chr18_mat" ) %>% 
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

p1<-ggplot(SUB, aes(x=midpoint/1000000, y=Dens/1000, color=NonB, fill=NonB)) +
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
  ggtitle("B")
p1

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

pb<- p1/ plot_spacer() / p_dummy  + plot_layout(heights = c(5,-0.5, 0.2), guides = "collect") 
pb

################################################################################
# MERGE PLOT B AND C 

p_final<- (pa / pb / plot_spacer() / pc)+plot_layout(heights=c(1, 1.5, 0.2, 0.8))
p_final

outfile=paste(plotdir,"Fig2_nonBcontent.png", sep="")
ggsave(outfile,plot = p_final,scale = 1,dpi = 300,limitsize = TRUE,width=9,height=12)
outfile=paste(plotdir,"Fig2_nonBcontent.svg", sep="")
ggsave(outfile,plot = p_final,scale = 1,dpi = 600,limitsize = TRUE,width=9,height=12)



