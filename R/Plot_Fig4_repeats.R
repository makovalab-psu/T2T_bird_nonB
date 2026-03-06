################################################################################
### R code for plotting non-B motif enrichment in repeats in Zebra finch.
### written by Linnéa Smeds, 2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(viridis)
require(patchwork)
require(ggtext)
library(dplyr)
library(ggplot2)
library(cowplot)

plotdir="plots/"
setwd("/Users/lbs5874/OneDrive - The Pennsylvania State University/Documents/Projects/ZebraFinch/")

################################################################################
################################################################################
# PART A, General repeat enrichment 

# Input files 
repfile="repeats/bTaeGut7v0.4_MT_rDNA.enrichment.tsv"
lenfile="repeats/TE_TRF_SAT_lengths.ordered.txt"
repgroupfile="repeats/bTaeGut7v0.4_MT_rDNA.enrichment.group.tsv"
lenfile1="repeats/macro.TE_TRF_SAT_lengths.txt"
lenfile2="repeats/micro.TE_TRF_SAT_lengths.txt"
lenfile3="repeats/dot.TE_TRF_SAT_lengths.txt"

# Reading in the data 
rep_tib1<-repfile %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(Group="GW")
rep_tib2<-repgroupfile %>% read.table(header=TRUE) %>% as_tibble()
len_tib<-lenfile %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Order=V1, Repeat=V2, GW=V3)
len_tib1<-lenfile1 %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Repeat=V1, macro=V2)
len_tib2<-lenfile2 %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Repeat=V1, micro=V2)
len_tib3<-lenfile3 %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Repeat=V1, dot=V2)

LEN<-len_tib %>% full_join(len_tib1) %>% full_join(len_tib2) %>% 
  full_join(len_tib3) %>% mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>%
  mutate(Lenfilter=if_else(GW>=10000, "Yes", "No")) %>%
  pivot_longer(cols=c(GW,macro, micro, dot), names_to="Group", values_to="Length")

DATA<- rep_tib1 %>% full_join(rep_tib2) %>% inner_join(LEN, by=c("Repeat","Group")) %>%
  mutate(Class=case_when(str_detect(Repeat, "DNA") ~ "TEs",
                         str_detect(Repeat, "LINE") ~ "TEs",
                         str_detect(Repeat, "SINE") ~ "TEs",
                         str_detect(Repeat, "MITE") ~ "TEs",
                         str_detect(Repeat, "LTR") ~ "TEs",
                         str_detect(Repeat, "Helitron") ~ "TEs",
                         str_detect(Repeat, "Retroposon") ~ "TEs",
                         str_detect(Repeat, "mer") ~ "TRs",
                         str_detect(Repeat, "Tgut") ~ "Satellites",
                         str_detect(Repeat, "TEL") ~ "Satellites",
                         TRUE ~ "Other"))  %>%
  mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                  Group=="micro" ~ "Micro",
                                  Group=="dot" ~ "Dot",
                                  TRUE ~ Group)) %>%
  mutate(Repeat=case_when(Repeat=="TEL" ~ "<b>Telomeric</b>",
                          Repeat=="Tgut191A" ~ "<b>Tgut191A</b>",
                          Repeat=="Tgut717A" ~ "<b>Tgut716A</b>",
                          Repeat=="LTR/Ngaro" ~ "<b>LTR/Ngaro</b>",
                         TRUE ~ Repeat)) %>%
  mutate(PrintName=paste(Repeat," (",sprintf("%.2f",Length/1000000),"Mb)", sep="")) %>%
  mutate(Enrichment_gw=if_else(Enrichment_gw==0, 0.0001, Enrichment_gw)) %>%
  mutate(logval=log2(Enrichment_gw)) %>% 
  mutate(nonB=if_else(nonB=="ALL", "Any", nonB)) %>% 
  mutate(PrintName=reorder(PrintName,Order))
  
  
DATA$nonB <- factor(DATA$nonB, levels=c("APR", "DR", "STR", "IR","TRI", "G4", "Z","Any"))
#DATA$nonB <- factor(DATA$nonB, levels=c("APR", "DR", "DR50-5", "G4", "g4discovery", "IR", "IR30-10", "IR_SM15", "MR", "TRI", "MR0-10", "MR0-15", "MR01-8", "STR", "Z", "Zseeker", "ZDNAHunter", "Any"))
DATA$Class <- factor(DATA$Class, levels=c("TEs", "TRs", "Satellites", "RNA", "Other"))
DATA$Classification<-factor(DATA$Classification, levels=c("Macro", "Micro", "Dot"))
DATA<-DATA %>% filter(!is.na(nonB))

################################################################################
# PLOT REPEAT ENRICHMENT; ONLY Repeat types with more than 10kb in total, and skip 
# "Satellite" as a joint class 
# Also tried plotting without logscale, but that doesn't work at all, the 
# highest enrichment value is way too high. 

SUB1<-DATA %>% filter(Repeat!="Satellite/Satellite" & Group!="GW" & Lenfilter=="Yes")
#outfile=paste(plotdir,"Fig3_NonB_enrichment_repeats_log2_5kb_",group,".png", sep="")
outfile=paste(plotdir,"Fig4_NonB_enrichment_repeats_log2_facet.png", sep="")

scalestart=min(SUB1 %>% filter(Enrichment_gw>0) %>%select(logval))
scaleend=max(SUB1 %>% filter(Enrichment_gw>0) %>%select(logval))
scalebreak=-1*scalestart/(scaleend-scalestart)

pA<-ggplot(SUB1, aes(reorder(Repeat, Order), forcats::fct_rev(nonB))) +
  geom_tile(aes(fill=logval), color = "white",lwd = 0.3,linetype = 1) +
  facet_grid(Classification~Class, space="free", scales="free_x") +
  scale_fill_gradientn(name="log2(Enrichment)",
                       colors=c("darkblue", "white","red"),
                      # colors=c("#075AFF", "white","#FF0000"),
                       values=c(0, scalebreak, 1), na.value="grey")+
  labs(x="", y="") +
  theme(axis.text.x = element_markdown(angle = 55, vjust = 1, hjust=1, size=8),
        legend.position="right",
        legend.text=element_text(size=8),
        legend.title=element_text(size=10, hjust=0.5, angle=90),
        legend.title.position="left",
        plot.title=element_text(size=18, face="bold"),
        legend.key.width  = unit(0.5, "lines"),
        legend.key.height = unit(3, "lines"),
        legend.margin = margin(t=0,r=0,b=0,l=0, unit="pt"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text.x = element_text(colour = 'black', size=11, face="bold"),
        strip.text.y = element_text(colour = 'black', size=11),
        panel.background = element_rect(fill = 'white', colour="white"),
        panel.spacing.y=unit(2, "lines"),
     )+
  ggtitle("A")+
  scale_y_discrete(position="left")
pA

################################################################################
################################################################################
################################################################################
# PART B,C, Intron TR enrichment 


# Input files 
enfile="repeats/bTaeGut7v0.4_MT_rDNA.intron_enrichment.group.tsv"
lenfile="repeats/bTaeGut7v0.4_MT_rDNA.introns_TR.lengthsummary.tsv"
repfile="repeats/bTaeGut7v0.4_MT_rDNA.introns_TR.lengths.tsv"

# Reading in the data 
en_tib<-enfile %>% read.table(header=TRUE) %>% as_tibble()
lentib<-lenfile %>% read.table(header=TRUE) %>% as_tibble() 
reptib<-repfile %>% read.table(header=TRUE) %>% as_tibble() 

# Category and Compartment data together 
DATA <- en_tib  %>% 
  mutate(PrintName=case_when(Subset=="TRF" ~ "tandem repeat",  
                             Subset=="noTRF" ~ "outside tandem repeat",
                             TRUE ~ Subset)) %>%
  mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                  Group=="micro" ~ "Micro",
                                  Group=="dot" ~ "Dot",
                                  TRUE ~ Group)) 
DATA$NonB <- factor(DATA$NonB, levels=c("ALL","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))
DATA$PrintName <- factor(DATA$PrintName, levels=c("tandem repeat","outside tandem repeat"))

# Length with fractions of Bp
LEN <- lentib %>% group_by(Group) %>% 
  mutate(totbp=sum(Length), frac = Length/totbp) %>%
  ungroup() %>% mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                                Group=="micro" ~ "Micro",
                                                Group=="dot" ~ "Dot",
                                                TRUE ~ Group)) %>%
  mutate(ClassCol = ifelse(Subset == "noTRF", "No", as.character(Classification)))
LEN$Classification<-factor(LEN$Classification, levels=c("Macro", "Micro", "Dot"))
LEN$ClassCol<-factor(LEN$ClassCol, levels=c("No", "Macro", "Micro", "Dot"))


# Repeats 
REP <-reptib %>% mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                                 Group=="micro" ~ "Micro",
                                                 Group=="dot" ~ "Dot",
                                                 TRUE ~ Group)) 
REP$Classification<-factor(REP$Classification, levels=c("Macro", "Micro", "Dot"))


# Set colors 
vcolors=viridis(7)
classification_fill <- c("Macro" = "#896192","Micro" = "#B99AC1", "Dot" = "#E8DCEB", "No" = "white")
classification_colors <- c("Macro" = "#896192","Micro" = "#B99AC1", "Dot" = "#E8DCEB","No" = "lightgray")
################################################################################
# PLOT FUNCTIONS 
################################################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO MAKE ENRICHMENT BAR PLOTS 
make_plot <- function(data, y_max = NULL, classif = NULL) {
  p<-ggplot(data, aes(x=PrintName, y=Enrichment_gw, fill=NonB, color=NonB)) + 
    geom_bar(stat="identity", position=position_dodge(), alpha=0.6) +
    scale_fill_manual(values=vcolors) +
    scale_color_manual(values=vcolors) +
    geom_hline(yintercept=1, linetype="dashed", color = "red") +
    scale_x_discrete(expand = c(0.08, 0.08)) +
    labs(x = NULL, y = "Fold enrichment") +  
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      legend.position = "none", # We'll keep legend only on one
      axis.text.x = element_text(size=12),
      panel.grid = element_blank(),
      axis.ticks.y=element_line(),
      legend.title=element_blank(),
      plot.title=element_text(size=18, face="bold"),
    )+
    ggtitle("B")
  if (!is.null(y_max)) {
    p <- p + 
      coord_cartesian(ylim = c(0, y_max)) +
      scale_y_continuous(expand = c(0.0, 0.0))
    
  }
  if (!is.null(classif)) {
    p <- p +
      annotate("text",
               x = Inf,
               #x = median(as.numeric(factor(data$PrintName))),  # midpoint of x-axis
               y = y_max,  
               label = classif,
               hjust = 1.2, vjust = 1.6,
               size = 5 , alpha = 0.65)  # subtle label
  }
  return(p)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO MAKE PIE CHARTS 
make_pie <- function(data) {
  p<-ggplot(data, aes(x = 1, y = frac, fill = ClassCol, color = ClassCol)) +
    geom_bar(stat = "identity", width = 1, alpha=0.6, show.legend = FALSE) +  # single bar per pie
    coord_polar(theta = "y", clip = "off") +  # pie chart
    geom_text(data=data%>%filter(Subset=="TRF"), aes(x = 1, y = frac/2, label = paste(sprintf("%.1f",frac*100),"%", sep="")),
              color = "black", size = 4, fontface = "italic", vjust=-1)+
    facet_wrap(Classification~., nrow=3)+
    theme_minimal() +
    labs(x = NULL, y = NULL) +  
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      strip.text=element_blank(),
      axis.title=element_blank(),
      plot.title=element_text(size=18, face="bold"),
    ) +
    scale_fill_manual(values = classification_fill) +
    scale_color_manual(values = classification_colors)
    #ggtitle("C")
  return(p)
} 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO PLOT REPEAT LENGTH DISTRIBUTION

make_histo <- function(data) {
  
  p <- ggplot(data%>%filter(RepUnitLen<100), aes(x = RepUnitLen, fill = Classification, color = Classification)) + 
    geom_histogram(binwidth=1, alpha = 0.6, show.legend = FALSE) +
    facet_wrap(Classification ~., nrow=3,  strip.position = "right") +  # single-column facet
    scale_fill_manual(values = classification_fill) +
    scale_color_manual(values = classification_colors) +
    labs(x="Tandem repeat unit length (bp)", y="Counts")+
    theme(
      panel.grid = element_line(color="gray95"),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      strip.text.y = element_text(size = 14),
      strip.background = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.ticks.length = unit(4, "pt"),
      plot.title=element_text(size=18, face="bold"),
      #    panel.spacing.y = unit(1, "lines")
    )+
    ggtitle("C")
  return(p)
}

################################################################################
# EXTRACT SUBSETS OF DATA, CHROMOSOME CATEGORIES 
################################################################################

# EXTRACT CHROMOSOME CATEGORIES 
# Enrichment, first panels 
data1_macro <- DATA %>% filter(Classification == "Macro") %>%
  filter(NonB!="ALL")
data1_micro <- DATA %>% filter(Classification == "Micro") %>%
  filter(NonB!="ALL")
data1_dot <- DATA %>% filter(Classification == "Dot") %>% 
  filter(NonB!="ALL")

# Length, for pie
data2_macro <- LEN %>% filter(Classification == "Macro")
data2_micro <- LEN %>% filter(Classification == "Micro")
data2_dot <- LEN %>% filter(Classification == "Dot")

################################################################################
#  MAKE 7 SEPARATE PLOTS (THREE FOR EACH B, C)
################################################################################
# Create the individual plots
p1 <- make_plot(data1_macro,25) + theme(axis.text.x = element_blank()) +
  theme(legend.position = "top",  legend.justification = "right", axis.title.y= element_blank()) + 
  guides(fill = guide_legend(nrow = 1))
p2 <- make_plot(data1_micro,25) + theme(axis.text.x = element_blank(), plot.title=element_blank())
p3 <- make_plot(data1_dot,25) + theme(axis.title.y= element_blank(), plot.title=element_blank())
p4<-make_pie(data2_macro)
p5<-make_pie(data2_micro)
p6<-make_pie(data2_dot)
p7<-make_histo(REP)

# COMBINE 
pB1<- p1 + inset_element(p4, left = 0.7, bottom = 0.2, right = 0.99, top = 0.99)
pB2<- p2 + inset_element(p5, left = 0.7, bottom = 0.2, right = 0.99, top = 0.99)
pB3<- p3 + inset_element(p6, left = 0.7, bottom = 0.2, right = 0.99, top = 0.99)
pB <- (pB1 / pB2 / pB3)

# Combine plots with patchwork with different widths
pBC <- (pB | p7) + plot_layout(widths = c(2.5, 1))
pBC


################################################################################
################################################################################
################################################################################
# PART D, CD SPECTRA, Written by Jacob Sieg


####Define factor levels and palettes####
Buffer_levels = c("100 mM KCl 140 mM LiCl 20 mM MOPS pH 7.2", "140 mM LiCl 20 mM MOPS pH 7.2") #Intentionally hardcoded
Buffer_labels = c("100 mM KCl 140 mM LiCl", "140 mM LiCl")
Buffer_palette =c("#DB5829", "#F4A637")
Sample_levels = c("Blank", "Tgut368A", "CR1_LINE", "CR1_A", "CR1_B") #Intentionally hardcoded
Sample_labels = c("Blank", "Tgut368A", "CR1 LINE", "CR1 A", "CR1 B")
Wavelength_limits = c(200, 320)
####Read in data####
df_CD = read.csv("js3025_CD_of_Linneas_bird_satellite_pG4s/js3025_calculations_spreadsheet - Preprocessed_CD_Data.csv")
### Set factor levels ###
df_CD$Buffer = factor(df_CD$Buffer, levels = Buffer_levels, labels = Buffer_labels)

####Function that plots data for each Sample####
plot_data = function(x = "CR1_A", panelname){
  df_CD_x = df_CD %>% filter(Sample == x)
  
  P_CD = ggplot(df_CD_x, aes(x = Wavelength_nm, y = Molar_ellipticity, color = Buffer, linetype = Buffer)) +
    geom_line() +
    theme_classic() +
    scale_color_manual(values = Buffer_palette) +
    scale_x_continuous(limits = Wavelength_limits, breaks = seq(Wavelength_limits[1], Wavelength_limits[2], by = 20)) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.text = element_text(size =10),
          legend.title = element_text(size = 10), 
          margins=margin(t=0,b=0,l=0,r=0),
          plot.title = element_text(size = 18, face="bold"),
          legend.position = "bottom") +
    xlab("") +
    ylab("")+
    ggtitle(panelname)
  return(P_CD)
}

####Plot and save data####
#"Tgut368A" "CR1 LINE", "CR1 A", "CR1 B"
pD1 = plot_data(x = "Tgut368A", "")+
  theme(legend.position = "none")
pD2 = plot_data(x = "CR1_LINE", "")+
  theme(legend.position = "none")
pD3 = plot_data(x = "CR1_A", "")+
  theme(legend.position = "none")
pD4 = plot_data(x = "CR1_B", "")+
  theme(legend.position = "none")

# Get Legend
legend_plot <- get_legend(
  plot_data(x = "CR1_A", "")+
    theme(legend.position = "bottom"))
legend_plot <- wrap_elements(legend_plot)

# Merge
pD=((pD1+pD2+pD3+pD4) + plot_layout(widths=c(1,1,1,1)))


# Make Y axis 
y_axis <- ggplot() +
  labs(y = "Molar ellipticity") +
  theme_void() +
  theme(
    axis.title.y = element_text(angle = 90, size = 14, vjust=0),
    plot.margin = margin(r = 0)
  )
# Make X axis 
x_axis <- ggplot() +
  labs(x = "Wavelength (nm)") +
  theme_void() +
  theme(
    axis.title.x = element_text(size = 14, vjust=0),
    plot.margin = margin(t = 0)
  )
# Add y first
pD_with_y <- (y_axis | pD) +
  plot_layout(widths = c(0.01, 1))
pD_with_x <- pD_with_y /
  x_axis +
  plot_layout(heights = c(1, 0.01))
pD_with_x 
# Add legend 
pD <- (wrap_elements(pD_with_x) / legend_plot) + 
  plot_layout(heights=c(1,0.05))

pD<-pD + plot_annotation(title = "D", theme=theme(plot.title = element_text(size = 18, face="bold")))
pD

################################################################################
# COMBINE THE THREE PARTS INTO ONE 
final=(pA/pBC/pD)+plot_layout(heights=c(3,3.5,3))
final
outfile=paste(plotdir,"Fig4_repeats.svg", sep="")
ggsave(outfile,plot = final,scale = 1,dpi = 600,limitsize = TRUE,width=11,height=15)
outfile=paste(plotdir,"Fig4_repeats.png", sep="")
ggsave(outfile,plot = final,scale = 1,dpi = 300,limitsize = TRUE,width=11,height=15)
# Add sequences in adobe











################################################################################
#OLD CODE, PLOT LENGTH DISTRIBUTION PER CLASS INSTEAD OF EXACT LENGTH
make_bar <- function(data, y_max = NULL) {
  p<-ggplot(data, aes(x = Repeat, y=Number, fill = Classification, color = Classification)) + 
    geom_bar(stat="identity", position=position_dodge(), alpha=0.6, show.legend = FALSE) +
    scale_fill_manual(values=classification_colors) +
    scale_color_manual(values=classification_colors) +
    scale_x_discrete(expand = c(0.08, 0.08)) +
    facet_wrap(Classification~., nrow=3, strip.position = "right")+
    labs(x ="", y="Annotated TR in introns") +  
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      strip.background = element_blank(),
      axis.text.x = element_text(size=10, angle=45, vjust=1, hjust=1),
      axis.title.y=element_text(size=12),
      panel.grid = element_blank(),
      strip.text = element_text(size=14),
      axis.ticks.x=element_blank(),
      plot.title=element_text(size=18, face="bold"),
    )+
    ggtitle("D")
  if (!is.null(y_max)) {
    p <- p + 
      coord_cartesian(ylim = c(0, y_max)) +
      scale_y_continuous(expand = c(0.0, 0.0))
    
  }
  return(p)
}
p7<-make_bar(data3_macro,160000)+ theme(axis.text.x = element_blank(), axis.title.y= element_blank()) 
p8<-make_bar(data3_micro,160000)+ theme(axis.text.x = element_blank(), plot.title=element_blank()) 
p9<-make_bar(data3_dot,160000)+ theme(axis.title.y = element_blank(), plot.title=element_blank()) 

