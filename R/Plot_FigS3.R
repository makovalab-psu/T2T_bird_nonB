################################################################################
### R code for plotting a comparison between zebra finch and chicken 
### written by Linnéa Smeds 
### September 4, 2025
################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(patchwork) 
plotdir="plots/"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INPUT FILES 
zdensfile="densities/bTaeGut7v0.4_MT_rDNA.nonB_per_chrom.txt"
zscaffile="ref/bTaeGut7v0.4_MT_rDNA.fa.fai"
zgcfile="stats/bTaeGut7v0.4_MT_rDNA.GC_per_chrom.txt"
zgroupfile="groups.txt"
cdensfile="densities/chicken.v23.nonB_per_chrom.txt"
cscaffile="ref/chicken.v23.fa.fai"
cgcfile="stats/chicken.v23.GC_per_chrom.txt"
cgroupfile="chicken.v23.groups.txt"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN DATA 
zdens_tib<-zdensfile %>% read.table(header=TRUE) %>% as_tibble()
zscaf_tib<- zscaffile %>% read.table(header=FALSE) %>% as_tibble() %>% 
            select(V1,V2) %>% rename(Chr=V1, Length=V2)
zgc_tib<- zgcfile %>% read.table(header=TRUE) %>% as_tibble() 
zgrouptib<-zgroupfile %>% read.table(header=TRUE) %>%
  mutate(Group=case_when(Group=="macro" ~ "Macro",
                         Group=="micro" ~ "Micro",
                         Group=="microdot" ~ "Dot",
                         TRUE ~ Group))
cdens_tib<-cdensfile %>% read.table(header=TRUE) %>% as_tibble()
cscaf_tib<- cscaffile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  select(V1,V2) %>% rename(Chr=V1, Length=V2)
cgc_tib<- cgcfile %>% read.table(header=TRUE) %>% as_tibble() 
cgrouptib<-cgroupfile %>% read.table(header=TRUE) %>%
  mutate(Group=case_when(Group=="macro" ~ "Macro",
                         Group=="micro" ~ "Micro",
                         Group=="microdot" ~ "Dot",
                         TRUE ~ Group))

tmpz<-zdens_tib %>% inner_join(zscaf_tib) %>% inner_join(zgc_tib) %>% 
  inner_join(zgrouptib) %>% mutate(Species="Zebra finch")
tmpc<-cdens_tib %>% inner_join(cscaf_tib) %>% inner_join(cgc_tib) %>% 
  inner_join(cgrouptib) %>% mutate(Species="Chicken")

DATA <- tmpz %>% bind_rows(tmpc)
DATA$NonB <- factor(DATA$NonB, levels=c("ALL","APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))
DATA$Group <- factor(DATA$Group, levels=c("Macro", "Micro", "Dot"))
DATA$Species <- factor(DATA$Species, levels=c("Zebra finch", "Chicken"))

# DEFINE COLORS 
type_col=c("#440154","#9460a1", "#e3cfe8")
################################################################################
# PLOT PANEL A, NON-B VS CHROMOSOME LENGTH  

# Just use "ALL"
SUBSET<-DATA %>% filter(NonB=="ALL")

pa<-ggplot(SUBSET, aes(x=Length/1000000, y=Density*100, fill=Group, color=Group))+
  geom_point(aes(shape=Group), show.legend=FALSE, alpha=0.6, size=3)+
  geom_segment(data=SUBSET%>%filter(Chr=="chrW"), aes(x = Length/1000000+15, y = Density*100+4, xend = Length/1000000+3, yend = Density*100+1), # Start and end coordinates of the arrow
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), # Adjust arrowhead length
               color = "red", # Set arrow color
               linewidth = 0.2) +
  labs(x="Chromosome Size (Mb)", y="Non-B motif content (%)") +
  scale_fill_manual(values=type_col)+
  scale_color_manual(values=type_col)+
  scale_shape_manual(values=c(21,24,22))+
  scale_x_log10()+
  facet_wrap(Species~.)+
  theme(panel.background = element_rect(fill = 'white', colour="black"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=14),
        panel.grid = element_blank(),
        panel.spacing=unit(1, "lines"),
        legend.background = element_rect(color = NA),
        legend.key = element_rect(color = NA),
        legend.position="right",
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        plot.title = element_text(size=18, face="bold"),
   )+
  ggtitle("A")+
  guides(fill = guide_legend(nrow = 3))
pa

#outfile=paste(plotdir,"FigS3A_NonB_density_vs_ChrLength.png", sep="")
#ggsave(outfile,plot = pc,scale = 1,dpi = 600,limitsize = TRUE,width=9,height=8)

################################################################################
# PLOT PANEL B, PLOT NON-B VS GC CONTENT  

pb<-ggplot(SUBSET, aes(x=GC_cont*100, y=Density*100, fill=Group, color=Group))+
  geom_point(aes(shape=Group), show.legend=TRUE, alpha=0.6, size=3)+
  geom_segment(data=SUBSET%>%filter(Chr=="chrW"), aes(x = GC_cont*100-4, y = Density*100+5, xend = GC_cont*100-1, yend = Density*100+1), # Start and end coordinates of the arrow
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), # Adjust arrowhead length
               color = "red", # Set arrow color
               linewidth = 0.2) +
  labs(x="GC content (%)", y="Non-B motif content (%)") +
  scale_fill_manual(values=type_col)+
  scale_color_manual(values=type_col)+ 
  scale_shape_manual(values=c(21,24,22))+
  facet_wrap(Species~.)+
  theme(panel.background = element_rect(fill = 'white', colour="black"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=14),
        panel.grid = element_blank(),
        panel.spacing=unit(1, "lines"),
        legend.background = element_rect(color = NA),
        legend.key = element_rect(color = NA),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.margin = margin(t=0,l=0,b=0,r=0),
        plot.title = element_text(size=18, face="bold"),
  )+
  guides(fill = guide_legend(nrow = 1))+
  ggtitle("B")
pb

combined_plot <- pa / pb
combined_plot
outfile=paste(plotdir,"FigS3_CompareChZF.png", sep="")
ggsave(outfile,plot = combined_plot,scale = 1,dpi = 600,limitsize = TRUE,width=9,height=10)
