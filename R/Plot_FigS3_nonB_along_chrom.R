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


# DEFINE COLORS 
viridis_colors <- viridis(7)
################################################################################
################################################################################ 
# PLOT FIGURE S3, COVERAGE ALONG CHROMOSOMES 

# INPUT FILES FOR B
windfile="coverage/bTaeGut7v0.4_MT_rDNA.merged.100kb.txt"
cenfile="ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff"
ABfile="ref/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.bed"

# DATA TIBBLES FOR B
centib<-cenfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  select(V1,V4,V5) %>% rename(Chr=V1, Cen_start=V4, Cen_stop=V5)

patterns<-c("_mat", "_pat")
replacements<-c("", "")

DATA<-windfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(midpoint=Start+(Stop-Start)/2) %>% inner_join(centib) %>% 
  mutate(PrintName=str_replace_all(Chr, "_mat", ""), 
         PrintName=str_replace_all(PrintName, "_pat", ""), 
         PrintName=str_replace_all(PrintName, "chr", "Chr")) %>%
  group_by(NonB) 
DATA$NonB <- factor(DATA$NonB, levels=c("APR", "DR", "STR", "IR", "TRI", "G4", "Z"))
AB<- ABfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Chr=V1, Start=V2, Stop=V3,Comp=V4) %>% 
  mutate(PrintName=str_replace_all(Chr, "_mat", ""), 
         PrintName=str_replace_all(PrintName, "_pat", ""), 
         PrintName=str_replace_all(PrintName, "chr", "Chr")) 


################################################################################
make_landscape_plot<-function(data1, data2, list) {
  # Get the centromeres 
  CEN<-data1 %>% select(PrintName, Cen_start, Cen_stop, max, NonB) %>% distinct() %>%
    mutate(Cen_mid=(Cen_start+Cen_stop)/2)
  
  # Make a tibble with only Chr Start and stop
  CHR<-data1 %>% group_by(PrintName) %>% mutate(Start=0, Stop=max(Stop)) %>%
    select(PrintName, Start, Stop) %>% unique() %>% ungroup()

  # Apply factorization to all data frames
  CEN$PrintName <- factor(CEN$PrintName, levels = list)
  CHR$PrintName <- factor(CHR$PrintName, levels = list)
  
  p1<-ggplot(data1, aes(x=midpoint/1000000, y=Dens/1000, color=NonB, fill=NonB)) +
    #  geom_rect(data=CEN, aes(xmin=Cen_start/1000000, ymin=0.0, ymax=max+0.01, xmax=Cen_stop/1000000), inherit.aes = FALSE, fill="gray", color="gray", linewidth=1, show.legend = FALSE) +
    geom_line(show.legend = FALSE) +
    geom_area(alpha = 0.7, show.legend = FALSE) +
    facet_grid(NonB~PrintName, scales="free", space="free_x")+
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
    scale_x_continuous(name='', expand = c(0, 0))
  p1
  
  # Add a dummy tract to plot the AB compartments and centromeres 
  p_dummy <- ggplot() +
    geom_rect(data=CHR,
              aes(xmin=Start/1e6, xmax=Stop/1e6, ymin=0.2, ymax=1),
              fill="white", color="lightgray") +
    geom_rect(data=data2, aes(xmin=Start/1e6, xmax=Stop/1e6, ymin=0.2, ymax=1, fill=Comp), color=NA, alpha=0.6,show.legend = FALSE) +
    geom_rect(data=CEN, aes(xmin=Cen_start/1000000, xmax=Cen_stop/1e6, ymin=0.2, ymax=1), fill='red', color=NA,show.legend = FALSE) +
    facet_grid(.~PrintName, scales="free_x", space="free_x") +
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
  
  p<- p1/ plot_spacer() / p_dummy  + plot_layout(heights = c(5,-0.5, 0.2), guides = "collect") 
  return(p)

} 

################################################################################
# make lists 
l1<-c("chr1_mat", "chr6_mat")
l2<-c("chr1A_mat", "chr4_mat")
l3<-c("chr2_mat")
l4<-c("chr3_mat", "chr7_mat")
l5<-c("chr4_mat", "chr5_mat")
l6<-c("chr4A_mat", "chr8_mat", "chr9_mat", "chr10_mat", "chr11_mat", "chr12_mat")
l7<-c("chr13_mat", "chr14_mat", "chr15_mat", "chr16_mat", "chr17_mat", "chr18_mat", "chr19_mat")
l8<-c("chr20_mat", "chr21_mat", "chr22_mat", "chr23_mat", "chr24_mat", "chr25_mat", "chr26_mat", "chr27_mat")
l9<-c("chr28_mat", "chr29_mat","chr30_mat","chr31_mat", "chr32_mat", "chr33_mat", "chr34_mat", "chr35_mat", "chr36_mat", "chr37_mat")
l10<-c("chrW_mat", "chrZ_pat")
 
ALL<-list(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)
plots <- vector("list", length(ALL))

for (i in 1:length(ALL)) {
  l=ALL[[i]]
    
  namelist<-gsub("_mat", "",l)
  namelist<-gsub("_pat", "",namelist)
  namelist<-gsub("chr", "Chr",namelist)
  
  # Subsample  
  SUB <- DATA %>% filter(Chr %in% l) %>% 
    mutate(max=max(Dens)/1000)
  
  SUB$PrintName<-factor(SUB$PrintName, levels = namelist)
  
  # Add AB regions 
  ABSUB<- AB %>% filter(Chr %in% l) %>% 
    select(Chr,Start,Stop,PrintName,Comp)
  ABSUB$PrintName<-factor(ABSUB$PrintName, levels = namelist)
  
  plots[[i]] <- make_landscape_plot(SUB, ABSUB, namelist)
  
}

################################################################################
# GROUP FIGURES TOGETHER AND PLOT 

page1 <- (wrap_elements(plots[[1]]) / wrap_elements(plots[[2]]) / wrap_elements(plots[[3]])) + plot_layout(heights=c(1,1,1))
page2 <- (wrap_elements(plots[[4]]) / wrap_elements(plots[[5]]) / wrap_elements(plots[[6]])) + plot_layout(heights=c(1,1,1))
page3 <- (wrap_elements(plots[[7]]) / wrap_elements(plots[[8]]) / wrap_elements(plots[[9]])) + plot_layout(heights=c(1,1,1))
page4 <- wrap_elements(plots[[10]])


outfile=paste(plotdir,"FigS3_part1_nonBcontent.png", sep="")
ggsave(outfile,plot = page1,scale = 1,dpi = 300,limitsize = TRUE,width=9,height=12)
outfile=paste(plotdir,"FigS3_part2_nonBcontent.png", sep="")
ggsave(outfile,plot = page2,scale = 1,dpi = 300,limitsize = TRUE,width=9,height=12)
outfile=paste(plotdir,"FigS3_part3_nonBcontent.png", sep="")
ggsave(outfile,plot = page3,scale = 1,dpi = 300,limitsize = TRUE,width=9,height=12)
outfile=paste(plotdir,"FigS3_part4_nonBcontent.png", sep="")
ggsave(outfile,plot = page4,scale = 1,dpi = 300,limitsize = TRUE,width=9,height=4)


