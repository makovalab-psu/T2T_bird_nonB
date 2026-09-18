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

# READ IN SPECIES LIST
spfile="helpfiles/species_list.txt"
sptib<-spfile %>% read.table(header=FALSE) %>% as_tibble()

# DEFINE COLORS 
type_col=c("#440154","#9460a1", "#e3cfe8")
viridis_colors <- viridis(7)
################################################################################

# DEFINE COLORS 
#type_col=c("#440154","#9460a1", "#e3cfe8", "lightgray")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function for plotting 

# ONE FACET FOR EACH NONB TYPE 
make_boxplot <- function(data1, data2, data3) {
  
  # Define pairwise comparisons for stats test
  groups <- as.character(unique(droplevels(data1$Group)))
  comparisons <- combn(groups, 2, simplify = FALSE)
  
  p<-ggplot(data1, aes(x=Group, y=Coverage*100, fill=Group))+
    geom_boxplot(show.legend=TRUE, alpha=0.9, color="black")+
    geom_boxplot(data = data3, show.legend=FALSE, alpha=0.9, color="grey", fill="lightgrey")+
    geom_hline(data = data2, aes(yintercept = Coverage*100), linetype="dashed", color = "red3")+
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
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)), limits = c(0, NA))+
    guides(fill = guide_legend(nrow = 1))+
    ggtitle("")
  return(p)
}

# ONE FACET FOR EACH SPECIES 
make_spboxplot <- function(data1, data2) {
  
  # Define pairwise comparisons for stats test
  groups <- unique(as.character(data1$Group))
  comparisons <- combn(groups, 2, simplify = FALSE)
  
  p<-ggplot(data1, aes(x=Group, y=Coverage*100, fill=Group))+
    geom_hline(data = data2, aes(yintercept = Coverage*100), linetype="dashed", color = "red3")+
    geom_boxplot(show.legend=TRUE, alpha=0.9, color="black")+
    facet_wrap(~Names, ncol=9)+
    scale_fill_manual(values=type_col)+
    scale_color_manual(values=type_col)+
    stat_compare_means(
      comparisons = comparisons, label = "p.signif", p.adjust.method = "BH", # FDR correction
      tip.length = 0.005, label.y.npc = "top", 
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")
      )
    ) +
    labs(x="", y="Coverage (%)", scale="free_y") +
    theme(panel.background = element_rect(fill = 'white', colour="black"),
          strip.background = element_rect(fill = 'white', colour="white"),
          strip.text = element_text(size=10),
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
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)))+
    coord_cartesian(ylim = c(0,85))+
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
    rename(NonB=V1, Bp=V2, Coverage=V3) %>% mutate(NonB=if_else(NonB=="Any", "All", NonB))
  gwcovtib$NonB <- factor(gwcovtib$NonB, levels=c("All","APR", "DR", "STR", "IR","TRI", "G4", "Z"))
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
    inner_join(typetib) %>% mutate(NonB=if_else(NonB=="Any", "All", NonB))
  DATA$NonB <- factor(DATA$NonB, levels=c("All","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))
  DATA$Group <- factor(DATA$Group, levels=unique(DATA$Group))
  DATA<-DATA %>% drop_na() %>% mutate(Species=sptib$V1[i])
  COMB<-bind_rows(COMB,DATA%>%select(-Bp,-Length))
  GWCOMB<-bind_rows(GWCOMB,gwcovtib%>%select(-Bp)%>%mutate(Species=sptib$V1[i]))
  
  if(sptib$V1[i]=="zebra_finch"){
    shuffile=paste("coverage/shuffledZFMerged.per_chrom.tsv", sep="")
    SHUF<-shuffile %>% read.table(header=TRUE) %>% as_tibble() %>% inner_join(typetib)
    SHUF$Group <- factor(SHUF$Group, levels=unique(SHUF$Group))
    SHUF$NonB <- factor(SHUF$NonB, levels=c("All","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))
  }
  else{
    SHUF<-tibble(Chr=character(), NonB=factor(), Bp=numeric(), Coverage=numeric(), Length=numeric(), Group=factor())
  }
  
  p<-make_boxplot(DATA, gwcovtib, SHUF)
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

ggsave(filename = paste0("plots/ZF_NonB_boxplot_with_shuffleShadows.png"),
  plot = pa, scale = 1,dpi = 300,limitsize = TRUE,width=12,height=4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#( ALL SPECIES IN ONE, ONLY USING TOTAL Non-B CONTENT )
ALL <- COMB %>% filter(NonB=="All") %>%
  mutate(Names=case_when(Species=="zebra_finch" ~ "Zebra finch",
                          Species=="ural_owl" ~ "Ural owl",
                          Species=="bandtailed_pigeon" ~ "Bandtailed\npigeon",
                          Species=="annas_hummingbird" ~ "Anna's\nhummingbird",
                          Species=="great_bustard" ~ "Great bustard",
                          Species=="chicken" ~ "Chicken",
                          Species=="golden_pheasant" ~ "Golden\npheasant",
                          Species=="peking_duck" ~ "Peking duck",
                          Species=="emu" ~ "Emu",
                          TRUE ~ Species)) %>%
  mutate(Names=factor(Names, levels=c("Zebra finch",
                                          "Ural owl",
                                          "Bandtailed\npigeon",
                                          "Anna's\nhummingbird",
                                          "Great bustard",
                                          "Chicken",
                                          "Golden\npheasant",
                                          "Peking duck",
                                          "Emu")))
ALL$Group <- factor(ALL$Group, levels=unique(ALL$Group))

ALLGW <- GWCOMB %>% filter(NonB=="All") %>%
  mutate(Names=case_when(Species=="zebra_finch" ~ "Zebra finch",
                         Species=="ural_owl" ~ "Ural owl",
                         Species=="bandtailed_pigeon" ~ "Bandtailed\npigeon",
                         Species=="annas_hummingbird" ~ "Anna's\nhummingbird",
                         Species=="great_bustard" ~ "Great bustard",
                         Species=="chicken" ~ "Chicken",
                         Species=="golden_pheasant" ~ "Golden\npheasant",
                         Species=="peking_duck" ~ "Peking duck",
                         Species=="emu" ~ "Emu",
                         TRUE ~ Species)) %>%
  mutate(Names=factor(Names, levels=c("Zebra finch",
                                      "Ural owl",
                                      "Bandtailed\npigeon",
                                      "Anna's\nhummingbird",
                                      "Great bustard",
                                      "Chicken",
                                      "Golden\npheasant",
                                      "Peking duck",
                                      "Emu")))

#PLOT ALL FOR PRESENTATION
pall<-make_spboxplot(ALL, ALLGW)
ggsave(filename = "plots/9sp_NonB_boxplot.png",
       plot = pall,
       scale = 1,dpi = 300,limitsize = TRUE,width=12,height=4)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT FIGURE 2C (remove ZF from above)

OTHER<-ALL %>% filter(Species!="zebra_finch")
OTHERGW<-ALLGW %>% filter(Species!="zebra_finch")
pc<-make_spboxplot(OTHER, OTHERGW)
pc=pc+ggtitle("C")
pc             


#ggsave(filename = paste0("plots/Fig2C_nonB_in_other_sp.png"),
#       plot = pc, scale = 1,dpi = 300,limitsize = TRUE,width=12,height=4)

################################################################################ 
# PLOT FIGURE 2B, COVERAGE ALONG EXAMPLE CHROMOSOMES 

# INPUT FILES FOR B
windfile="coverage/bTaeGut7v0.4_MT_rDNA.merged.100kb.txt"
repfile="coverage/bTaeGut7v0.4_MT_rDNA.repeats.100kb.bed"
genefile="coverage/bTaeGut7v0.4_MT_rDNA.coding_full_genes.100kb.bed"
gcfile="coverage/bTaeGut7v0.4_MT_rDNA.GC.100kb.bed"
cenfile="ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff"
ABfile="ref/bTaeGut7v0.4_MT_rDNA.Cooltools.E1.200kbp.flipped.dip.collated.v0.1.bed"

# DATA TIBBLES FOR B
centib<-cenfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  select(V1,V4,V5) %>% rename(Chr=V1, Start=V4, Stop=V5) %>% mutate(Compartment="CEN")
DATAB1<-windfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(midpoint=Start+(Stop-Start)/2) %>% 
  group_by(NonB)
DATAB1$NonB <- factor(DATAB1$NonB, levels=c("APR", "DR", "STR", "IR", "TRI", "G4", "Z"))

ABDOT<- ABfile %>% read.table(header=FALSE) %>% as_tibble() %>%
  select(V1,V2,V3,V4) %>% rename(Chr=V1, Start=V2, Stop=V3,Compartment=V4)
# Merge AB and Centromere regions 
MERGED<-ABDOT %>% bind_rows(centib)

# Read in and merge other coverage files
genetib<-genefile %>% read.table(header=FALSE) %>% as_tibble() %>% mutate(Type="Genes", Cov=V4/1000)
gctib<-gcfile %>% read.table(header=FALSE) %>% as_tibble() %>% mutate(Type="GC", Cov=V4*100)
reptib<-repfile %>% read.table(header=FALSE) %>% as_tibble() %>% mutate(Type="Repeats", Cov=V4/1000)
DATAB2<-genetib %>% bind_rows(gctib) %>% bind_rows(reptib) %>% 
    rename(Chr=V1, Start=V2, Stop=V3,Orig=V4) %>% mutate(midpoint=Start+(Stop-Start)/2) 
DATAB2$Type <- factor(DATAB2$Type, levels=c("GC", "Genes", "Repeats"))


# Choose 1 macro, 1 micro, 1 dot from the data 
SUBB1 <- DATAB1 %>% filter(Chr=="chr8_mat" | Chr=="chr36_mat" | Chr=="chr18_mat" ) %>% 
  mutate(ChrName=case_when(Chr=="chr8_mat" ~ "Chr8",
                           Chr=="chr18_mat" ~ "Chr18",
                           Chr=="chr36_mat" ~ "Chr36",
                           TRUE ~ Chr)) %>% 
  mutate(max=max(Dens)/1000)

SUBB2 <- DATAB2 %>% filter(Chr=="chr8_mat" | Chr=="chr36_mat" | Chr=="chr18_mat" ) %>% 
  mutate(ChrName=case_when(Chr=="chr8_mat" ~ "Chr8",
                           Chr=="chr18_mat" ~ "Chr18",
                           Chr=="chr36_mat" ~ "Chr36",
                           TRUE ~ Chr)) %>% 
  mutate(max=max(Cov)/1000)

# Make a tibble with only Chr Start and stop
CHR<-SUBB1 %>% group_by(ChrName) %>% mutate(Start=0, Stop=max(Stop)) %>%
  select(ChrName, Start, Stop) %>% unique() %>% ungroup()

# Add compartment regions 
ABC<-MERGED %>% filter(Chr=="chr8_mat"| Chr=="chr36_mat" | Chr=="chr18_mat" ) %>% 
  mutate(ChrName=case_when(Chr=="chr8_mat" ~ "Chr8",
                           Chr=="chr18_mat" ~ "Chr18",
                           Chr=="chr36_mat" ~ "Chr36",
                           TRUE ~ Chr)) %>% 
  select(Chr,Start,Stop,ChrName,Compartment)

# Apply to all data frames
chr_order <- c("Chr8", "Chr18", "Chr36")
SUBB1$ChrName <- factor(SUBB1$ChrName, levels = chr_order)
SUBB2$ChrName <- factor(SUBB2$ChrName, levels = chr_order)
CHR$ChrName <- factor(CHR$ChrName, levels = chr_order)
ABC$ChrName  <- factor(ABC$ChrName,  levels = chr_order)

p1<-ggplot(SUBB1, aes(x=midpoint/1000000, y=Dens/1000, color=NonB, fill=NonB)) +
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

# Other coverage features coverage 
p2<-ggplot(SUBB2, aes(x=midpoint/1000000, y=Cov)) +
  geom_line(show.legend = FALSE, color="grey") +
  geom_area(alpha = 0.7, show.legend = FALSE, fill="grey") +
  facet_grid(Type~ChrName, scales="free", space="free_x")+
  theme(panel.background = element_rect(fill = 'white', colour="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.title = element_blank(),
        strip.text.y = element_text(size=8, hjust=0.5, angle=-45),
        strip.text.x = element_blank(),
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
  scale_x_continuous(name='', expand = c(0, 0))
p2


# Add a dummy tract to plot the AB compartments and centromeres 
p_dummy <- ggplot() +
  geom_rect(data=CHR,
            aes(xmin=Start/1e6, xmax=Stop/1e6, ymin=0.2, ymax=1),
            fill="white", color="lightgray") +
  geom_rect(data=ABC, aes(xmin=Start/1e6, xmax=Stop/1e6, ymin=0.2, ymax=1, fill=Compartment), color=NA,show.legend = TRUE) +
  geom_point(
    data = ABC%>%filter(Compartment=="CEN"),
    aes(x = (Start + Stop)/2 / 1e6, y = 1.5),
    shape = 25,  # filled triangle pointing down
    size = 1,
    fill = "red",
    color = "red",
    inherit.aes = FALSE)+
  facet_grid(.~ChrName, scales="free_x", space="free_x") +
  scale_fill_manual(values=c("#E38AAA7f","#6F9DD07f", "red"))+
  theme_void() +
  theme(
    panel.background = element_rect(fill="white", colour="white"),
    panel.spacing.x = unit(0.5, "lines"),
    legend.position = "bottom",
    legend.justification = "right",
    strip.text.x =element_blank(),
    axis.title.x = element_text(size=14, vjust = 0.5), 
    axis.text.x  = element_text(size=10),
    axis.ticks.x = element_line(color="black", linewidth=0.3),
    axis.ticks.length.x = unit(4, "pt"),
    plot.margin  = margin(5,0,-10,0)
  )+
  scale_x_continuous(name='Position (Mb)', expand = c(0, 0))+
  coord_cartesian(ylim = c(0.2, 1.1), clip = "off")
p_dummy

pb<- p1/ p2 / plot_spacer() / p_dummy  + 
  plot_layout(heights = c(5, 2,-0.5, 0.2))
pb

################################################################################
# MERGE PLOT A, B AND C 

p_final<- (pa / pb / plot_spacer() / pc)+plot_layout(heights=c(1, 1.8, 0.15, 0.8))
p_final

outfile=paste(plotdir,"Fig2_nonBcontent.png", sep="")
ggsave(outfile,plot = p_final,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=13)
outfile=paste(plotdir,"Fig2_nonBcontent.svg", sep="")
ggsave(outfile,plot = p_final,scale = 1,dpi = 600,limitsize = TRUE,width=10,height=13)


#outfile=paste(plotdir,"Fig2B_nonBcontent.png", sep="")
#ggsave(outfile,plot = pb,scale = 1,dpi = 300,limitsize = TRUE,width=9,height=4)
#outfile=paste(plotdir,"Fig2C_nonBcontent.png", sep="")
#ggsave(outfile,plot = pc,scale = 1,dpi = 300,limitsize = TRUE,width=9,height=4)

