################################################################################
### R code for plotting non-B motif enrichment and methylation distributions in 
# functional regions in the Zebra finch T2T genome.
### written by Linnéa Smeds Sept 15, 2025.

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(viridis)
require(patchwork)
library(cowplot)  # for ggdraw()
library(grid)     # for textGrob()
require(ggh4x)
plotdir="plots/"
setwd("/Users/lbs5874/Documents/Projects/ZebraFinch/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input files 
# PANEL A
enfile="functional/enrichment_groups.v4.tsv"
varfile<-"functional/minmax.50perc.100rep.96CI.txt"
# PANEL B
g4file="functional/G4_strand_density_perGene.v5.tsv"
# PANEL C
methfile="methylation/functional/Group.allCpG.merged.txt"
intergenfile="methylation/functional/Group.allCpG.intergenic.txt"
# PANEL D
cenfile="functional/chicken.v23.enrichment_groups.tsv"
csigfile="functional/minmax.chicken.v23.50perc.100rep.96CI.txt"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reading in the data 
# PANEL A
vartib<-varfile %>% read.table(header=TRUE) %>% as_tibble()
DATAA <- enfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  inner_join(vartib) %>% 
  mutate(PrintName=case_when(Class=="nonprotcoding" ~ "lncRNA",  
                            Class=="promoter" ~ "1 kb upstream\nfrom TSS\n(putative promoters)",
                            Class=="all_repeats" ~ "Repeats",
                            Class=="intronic" ~ "Intronic",
                            Class=="intergenic" ~ "Intergenic",
                            Class=="UTR5" ~ "5'UTRs",
                            Class=="UTR3" ~ "3'UTRs",
                            TRUE ~ Class)) %>%
  mutate(Group=case_when(Group=="macro" ~ "Macro",
                                  Group=="micro" ~ "Micro",
                                  Group=="microdot" ~ "Dot",
                                  TRUE ~ Group)) 
DATAA$nonB <- factor(DATAA$nonB, levels=c("ALL","APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))
DATAA$PrintName <- factor(DATAA$PrintName, levels=c("1 kb upstream\nfrom TSS\n(putative promoters)","5'UTRs","CDS", "3'UTRs", "Intronic", "lncRNA","Intergenic"))

# PANEL B
DATAB <- g4file %>% read.table(header=TRUE) %>% as_tibble() %>%
  add_row(Class="intergenic", Group="Macro", Strand="Coding", Gene="dummy", GenLen=0, G4Len=0, Density=0) %>% 
  add_row(Class="intergenic", Group="Micro", Strand="Coding", Gene="dummy", GenLen=0, G4Len=0, Density=0) %>% 
  add_row(Class="intergenic", Group="Dot", Strand="Coding", Gene="dummy", GenLen=0, G4Len=0, Density=0) %>% 
  add_row(Class="intergenic", Group="Macro", Strand="Template", Gene="dummy", GenLen=0, G4Len=0, Density=0) %>% 
  add_row(Class="intergenic", Group="Micro", Strand="Template", Gene="dummy", GenLen=0, G4Len=0, Density=0) %>% 
  add_row(Class="intergenic", Group="Dot", Strand="Template", Gene="dummy", GenLen=0, G4Len=0, Density=0) %>% 
  mutate(PrintName=case_when(Class=="nonprotcoding" ~ "lncRNA",  
                             Class=="promoter" ~ "1 kb upstream\nfrom TSS",
                             Class=="all_repeats" ~ "Repeats",
                             Class=="intronic" ~ "Intronic",
                             Class=="UTR5" ~ "5'UTRs",
                             Class=="UTR3" ~ "3'UTRs",
                             Class=="intergenic" ~ "Intergenic",
                             TRUE ~ Class)) %>%
  mutate(Group=case_when(Group=="macro" ~ "Macro",
                                  Group=="micro" ~ "Micro",
                                  Group=="microdot" ~ "Dot",
                                  TRUE ~ Group)) 
DATAB$PrintName <- factor(DATAB$PrintName, levels=c("1 kb upstream\nfrom TSS","5'UTRs","CDS", "3'UTRs", "Intronic", "lncRNA","Intergenic"))

# PANEL C
methtib<-methfile %>% read.table(header=TRUE) %>% as_tibble() 
intergentib<-intergenfile %>% read.table(header=TRUE) %>% as_tibble() 
#sumtib<-sumfile %>% read.table(header=TRUE) %>% as_tibble() 


# Summarize gene regions with median and mean per transcript 
DATAC1 <-methfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  group_by(Group,Class,Type,Trx) %>% 
  summarize(median=median(Score), mean=mean(Score)) %>% 
  mutate(PrintName=case_when(Class=="nonprotcoding" ~ "lncRNA",  
                             Class=="promoter" ~ "1 kb upstream\nfrom TSS",
                             Class=="all_repeats" ~ "Repeats",
                             Class=="intronic" ~ "Intronic",
                             Class=="intergenic" ~ "Intergenic",
                             Class=="UTR5" ~ "5'UTRs",
                             Class=="UTR3" ~ "3'UTRs",
                             TRUE ~ Class)) %>%
  mutate(Group=case_when(Group=="macro" ~ "Macro",
                          Group=="micro" ~ "Micro",
                          Group=="microdot" ~ "Dot",
                            TRUE ~ Group)) %>% 
  mutate(PrintType=case_when(Type=="coding" ~ "Coding",
                             Type=="template" ~ "Template",
                             Type=="background" ~ "Background",
                             TRUE ~ Group)) %>% 
  mutate(Region=case_when(Type=="coding" ~ "inside G4",
                          Type=="template" ~ "inside G4",
                          Type=="background" ~ "outside G4",
                          TRUE ~ NA))
DATAC1$Group <- factor(DATAC1$Group, levels=c("Macro","Micro","Dot"))
DATAC1$PrintType <- factor(DATAC1$PrintType, levels=c("Background","Coding","Template"))
DATAC1$PrintName <- factor(DATAC1$PrintName, levels=c("1 kb upstream\nfrom TSS","5'UTRs","CDS", "3'UTRs", "Intronic", "lncRNA", "Intergenic"))
DATAC1$Region <- factor(DATAC1$Region, levels=c("inside G4","outside G4"))

# For intergenic we use all sites 
DATAC2<-intergentib %>%
  mutate(Group=case_when(Group=="macro" ~ "Macro",
                          Group=="micro" ~ "Micro",
                          Group=="microdot" ~ "Dot",
                          TRUE ~ Group)) %>% 
  mutate(Region=case_when(Type=="ignorant" ~ "inside G4",
                          Type=="background" ~ "outside G4",
                          TRUE ~ NA)) %>% mutate(PrintName="Intergenic", median=Score)
DATAC2$Group <- factor(DATAC2$Group, levels=c("Macro","Micro","Dot"))

# PANEL D
csigtib<-csigfile %>% read.table(header=TRUE) %>% as_tibble()
DATAD <- cenfile %>% read.table(header=TRUE) %>% as_tibble() %>% inner_join(csigtib) %>%
  mutate(PrintName=case_when(Class=="promoter" ~ "1 kb upstream\nfrom TSS\n(putative promoters)",
                             Class=="all_repeats" ~ "Repeats",
                             Class=="introns" ~ "Intronic",
                             Class=="intergenic" ~ "Intergenic",
                             Class=="UTR5" ~ "5'UTRs",
                             Class=="UTR3" ~ "3'UTRs",
                             TRUE ~ Class)) %>%
  mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                  Group=="micro" ~ "Micro",
                                  Group=="microdot" ~ "Dot",
                                  Group=="Fullgenome" ~ "Genome-wide",
                                  TRUE ~ Group)) 
DATAD$nonB <- factor(DATAD$nonB, levels=c("ALL","APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))
DATAD$PrintName <- factor(DATAD$PrintName, levels=c("1 kb upstream\nfrom TSS\n(putative promoters)","5'UTRs","CDS", "3'UTRs", "Intronic","Intergenic"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set colors 
vcolors=viridis(8)
type_col=c("#440154","#9460a1", "#e3cfe8")
################################################################################
# CODE FOR FIGURE 2A, ENRICHMENT IN FUNCTIONAL REGIONS 
################################################################################
# (Make three separate plots and combine with patchwork)

data_macro <- DATAA %>% filter(Group == "Macro") %>%
  filter(nonB!="ALL") %>% filter(PrintName!="Repeats")
data_micro <- DATAA %>% filter(Group == "Micro") %>%
  filter(nonB!="ALL") %>% filter(PrintName!="Repeats")
data_dot <- DATAA %>% filter(Group == "Dot") %>% 
  filter(nonB!="ALL") %>% filter(PrintName!="Repeats")

# Function to make consistent bar plots
make_plot <- function(data, y_max = NULL, classif = NULL) {
  p<-ggplot(data, aes(x=PrintName, y=Enrichment_gw, fill=nonB, color=nonB)) + 
    geom_bar(stat="identity", position=position_dodge(), alpha=0.6) +
    geom_errorbar(stat="identity", position=position_dodge(), aes(x=PrintName, ymin=Min, ymax=Max), colour="grey40")+
    scale_fill_manual(values=vcolors) +
    scale_color_manual(values=vcolors) +
    geom_hline(yintercept=1, linetype="dashed", color = "red") +
    scale_x_discrete(expand = c(0.08, 0.08)) +
    labs(x = NULL, y = "Fold enrichment") +  
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      legend.position = "none", # We'll keep legend only on one
      axis.text.x = element_text(size=10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      axis.title.y=element_text(size=12),
      axis.ticks.y=element_line()
    )
  if (!is.null(y_max)) {
    p <- p + 
      coord_cartesian(ylim = c(0, y_max)) +
      scale_y_continuous(expand = c(0.0, 0.0), breaks = seq(0, y_max, by = 3))
    
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

# Create the individual plots
pa1 <- make_plot(data_macro, 6, "Macro")  + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top",  
        legend.justification = "right",
        plot.title=element_text(size=18, face="bold", hjust=0)) + 
  ggtitle("A")+ guides(fill = guide_legend(nrow = 1))
pa2 <- make_plot(data_micro, 8, "Micro") + theme(axis.text.x = element_blank())
pa3 <- make_plot(data_dot, 16.1, "Dot") + theme(axis.title.y = element_blank())

# Combine plots with patchwork with different heights
pa<-(pa1 / pa2 / pa3) + plot_layout(heights = c(6, 8, 16.1))
pa

################################################################################
# CODE FOR FIGURE 2B, G4s ON CODING AND TEMPLATE STRANDS
################################################################################
# PLOT G4 STRAND COVERAGE
# (Make three separate plots and combine with patchwork)

data_macro2 <- DATAB %>% filter(Group == "Macro") 
data_micro2 <- DATAB %>% filter(Group == "Micro") 
data_dot2 <- DATAB %>% filter(Group == "Dot") 

make_bplot_with_stats <- function(data, y_max, classif = NULL) {
  # Summary of the data
  data_summary <- data %>%
    group_by(PrintName, Strand) %>%
    summarise(mean_density = mean(Density, na.rm = TRUE), .groups = "drop")
  # Paired wilcoxons test
  wilcox_results <- data %>%
    # first collapse duplicates per gene/strand/PrintName
    group_by(PrintName, Gene, Strand) %>%
    summarise(Density = mean(Density, na.rm = TRUE), .groups = "drop") %>%
    # now pivot to wide so each gene has Coding & Template in columns
    pivot_wider(names_from = Strand, values_from = Density) %>%
    group_by(PrintName) %>%
    summarise(
      p_value = {
        tmp <- cur_data() %>% filter(!is.na(Coding), !is.na(Template))
        if (nrow(tmp) > 0) {
          wilcox.test(tmp$Coding, tmp$Template, paired = TRUE, exact = FALSE)$p.value
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr"),
           p_signif = case_when(
                      is.na(p_adj)      ~ "",
                      p_adj <= 0.001    ~ "***",
                      p_adj <= 0.01     ~ "**",
                      p_adj <= 0.05     ~ "*",
                      TRUE                ~ "ns"
  ))
  #And plot 
  p<-ggplot(data_summary, aes(x = PrintName, y = mean_density*100,
                           fill = Strand, color = Strand)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.6) +
    scale_fill_manual(values = c(vcolors[7],"skyblue")) +
    scale_color_manual(values = c(vcolors[7],"blue")) +
    # Add significance text
    geom_text(
      data = wilcox_results,
      aes(x = PrintName, y = y_max - 1.5, 
          label = p_signif),
      inherit.aes = FALSE,
      vjust = 0
    )+
    scale_x_discrete(expand = c(0.08, 0.08)) +
    labs(x = NULL, y = "G4 Coverage (%)") +  
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      legend.position = "none", # We'll keep legend only on one
      axis.text = element_text(size=10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      title = element_blank(),
      axis.ticks.y=element_line(),
      axis.title.y=element_text(size=12),
    )
  if (!is.null(y_max)) {
    p <- p + 
      coord_cartesian(ylim = c(0, y_max)) +
      scale_y_continuous(expand = c(0.0, 0.0), breaks = seq(0, y_max, by = 2))
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

# Get max values using a function 
get_max <- function(data) {
  (data %>% group_by(PrintName, Strand, Group) %>%
    summarise(mean_density = mean(Density, na.rm = TRUE), .groups="drop") %>% 
    group_by(Group) %>% summarize(max=max(mean_density)))$max+0.02
}
  
max1<-get_max(data_macro2)*100
max2<-get_max(data_micro2)*100
max3<-get_max(data_dot2)*100

# Create the individual plots
pb1 <- make_bplot_with_stats(data_macro2, max1, "Macro")  + 
  theme(axis.text.x = element_blank(),
        legend.position = "top",  
        legend.justification = "right",
        axis.title.y = element_blank(),
        plot.title=element_text(size=18, face="bold", hjust=0)) + 
  ggtitle("B")+
  guides(fill = guide_legend(nrow = 1))
pb2 <- make_bplot_with_stats(data_micro2, max2, "Micro") + theme(axis.text.x = element_blank())
pb3 <- make_bplot_with_stats(data_dot2, max3, "Dot") + theme(axis.title.y = element_blank())

# Combine plots with patchwork with different heights
pb<-(pb1 / pb2 / pb3) + 
  plot_layout(heights = c(max1, max2, max3))
pb

################################################################################
# CODE FOR FIGURE 2C, METHYLATION IN AND OUTSIDE G4s 
################################################################################

# Calculate pvals in a function 
get_pval <- function(data) {
  pvals <- data %>%
    group_by(Group, PrintName, Class) %>%
    summarize(
      p = tryCatch(
        wilcox.test(median ~ Region)$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p, method = "fdr"),
      signif = case_when(
        is.na(p_adj)        ~ "n.s.",
        p_adj <= 0.001      ~ "***",
        p_adj <= 0.01       ~ "**",
        p_adj <= 0.05       ~ "*",
        TRUE            ~ "n.s."
      ),
      x = 0.5,          # Fixed x position for label
      y = Inf           # Top of the plot
    )
  return(pvals)
}

# Make a function that plots methylation inside vs outside of G4s 
make_densityplot <- function(data, sigdata, showleg = TRUE, showx = TRUE, showy=TRUE, ytext=FALSE) {
  p <- ggplot(data, aes(x = median/100, fill = Region, color = Region)) + 
    geom_density(alpha = 0.6, show.legend = showleg) +
    facet_nested(Group ~ PrintName, scales = "free_y") +  # single-column facet
    scale_fill_manual(values = c(vcolors[7],"gray90")) +
    scale_color_manual(values = c(vcolors[7],"gray30")) +
    geom_text(
      data = sigdata,
      aes(x = x, y = y, label = signif),
      inherit.aes = FALSE,
      vjust = 1.5,       # Adjust vertical position from the top edge
      size = 4
    )+
    labs(x="Methylation level", y="Density")+
    scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      plot.title=element_text(size=18, face="bold"), 
      strip.text.x = element_text(size = 10),
      strip.text.y = if(ytext==TRUE) element_text(size = 12) else element_blank(),
      strip.background = element_blank(),
      legend.key = element_rect(fill = NA, color = NA),
      axis.title.x = if(showx==TRUE) element_text(size = 12) else element_blank(),
      axis.title.y = if(showy==TRUE) element_text(size = 12) else element_blank(),
      axis.text.y = element_text(size = 10),
      axis.ticks.length = unit(4, "pt"),
      #    panel.spacing.y = unit(1, "lines")
    ) 
  return(p)
}

# Calculate pvalues 
pvals1<-get_pval(DATAC1)
pvals2<-get_pval(DATAC2)

# Plot each class panel and combine 
pc1<-make_densityplot(DATAC1%>%filter(Class=="promoter"), pvals1%>%filter(Class=="promoter"),FALSE, FALSE, TRUE)+ 
  theme(plot.title=element_text(size=18, face="bold")) + 
  ggtitle("C")
pc2<-make_densityplot(DATAC1%>%filter(Class=="UTR5"), pvals1%>%filter(Class=="UTR5"), FALSE, FALSE, FALSE)
pc3<-make_densityplot(DATAC1%>%filter(Class=="CDS"), pvals1%>%filter(Class=="CDS"), FALSE, FALSE, FALSE)
pc4<-make_densityplot(DATAC1%>%filter(Class=="UTR3"), pvals1%>%filter(Class=="UTR3"), TRUE, TRUE, FALSE)
pc5<-make_densityplot(DATAC1%>%filter(Class=="intronic"), pvals1%>%filter(Class=="intronic"), FALSE, FALSE, FALSE)
pc6<-make_densityplot(DATAC1%>%filter(Class=="nonprotcoding"), pvals1%>%filter(Class=="nonprotcoding"), FALSE, FALSE, FALSE)
pc7<-make_densityplot(DATAC2, pvals2, FALSE, FALSE, FALSE, TRUE)

pc<-pc1+pc2+pc3+pc4+pc5+pc6+pc7 + plot_layout(widths=c(1,1,1,1,1,1,1))
pc

################################################################################
# CODE FOR PANEL D, FUNCTIONAL ENRICHMENT IN CHICKEN
################################################################################
# (Make three separate plots and combine with patchwork)

data_macro4 <- DATAD %>% filter(Classification == "Macro") %>%
  filter(nonB!="ALL") %>% filter(PrintName!="Repeats")
data_micro4 <- DATAD %>% filter(Classification == "Micro") %>%
  filter(nonB!="ALL") %>% filter(PrintName!="Repeats")
data_dot4 <- DATAD %>% filter(Classification == "Dot") %>% 
  filter(nonB!="ALL") %>% filter(PrintName!="Repeats")

# Function to make consistent bar plots
make_plot <- function(data, y_max = NULL, classif = NULL) {
  p<-ggplot(data, aes(x=PrintName, y=Enrichment_gw, fill=nonB, color=nonB)) + 
    geom_bar(stat="identity", position=position_dodge(), alpha=0.6) +
    geom_errorbar(stat="identity", position=position_dodge(), aes(x=PrintName, ymin=Min, ymax=Max), colour="grey40")+
    scale_fill_manual(values=vcolors) +
    scale_color_manual(values=vcolors) +
    geom_hline(yintercept=1, linetype="dashed", color = "red") +
    scale_x_discrete(expand = c(0.08, 0.08)) +
    labs(x = NULL, y = "Fold Enrichment") +  
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      legend.position = "none", # We'll keep legend only on one
      axis.text.x = element_text(size=10),
      axis.title.y =element_text(size=12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      title = element_blank(),
      axis.ticks.y=element_line()
    )
  if (!is.null(y_max)) {
    p <- p + 
      coord_cartesian(ylim = c(0, y_max)) +
      scale_y_continuous(expand = c(0.0, 0.0), breaks = seq(0, y_max, by = 3))
    
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

pd1 <- make_plot(data_macro4, 6, "Macro")  + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top",  
        legend.justification = "right",
        plot.title=element_text(size=18, face="bold", hjust=0)) + 
  ggtitle("D") + guides(fill = guide_legend(nrow = 1))
pd2 <- make_plot(data_micro4, 8, "Micro") + theme(axis.text.x = element_blank())
pd3 <- make_plot(data_dot4, 19.5, "Dot") + theme(axis.title.y = element_blank())

# Combine plots with patchwork with different heights
pd<-(pd1 / pd2 / pd3) + plot_layout(heights = c(6, 8, 19.5))
pd

################################################################################
# ADD FOUR PANELS TOGETHER 

#
p_final<- (wrap_elements(pa)/ wrap_elements(pb) / wrap_elements(pc) /wrap_elements(pd))
p_final
outfile=paste(plotdir,"Fig2.png", sep="")
ggsave(outfile,plot = p_final,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=20)
outfile=paste(plotdir,"Fig2.svg", sep="")
ggsave(outfile,plot = p_final,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=20)

#p_final<- ((pa+pb) + plot_layout(width=c(3,3))) / pc
#outfile=paste(plotdir,"Fig2ABC_square_biggerA.png", sep="")
#ggsave(outfile,plot = p_final,scale = 1,dpi = 600,limitsize = TRUE,width=12,height=12)

