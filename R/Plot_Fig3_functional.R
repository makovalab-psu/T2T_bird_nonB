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
setwd("/Users/lbs5874/OneDrive - The Pennsylvania State University/Documents/Projects/ZebraFinch/")

# SPECIES LIST 
spfile="helpfiles/species_list.txt"
sptib<-spfile %>% read.table(header=FALSE) %>% as_tibble()

# Set colors 
vcolors=viridis(7)
type_col=c("#440154","#9460a1", "#e3cfe8")

################################################################################
# FUNCTION TO MAKE ENRICHMENT PLOT FOR PANEL 2A AND SUPPLEMENTARY FIGURES
################################################################################

# Function to make consistent bar plots
make_enrich_plot <- function(data, y_max = NULL, classif = NULL) {
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAKE THE ENRICHMENT PLOT FOR EACH SPECIES 
# Prepare list for saving plots 
plots <- vector("list", length(sptib$V3))
names(plots) <- sptib$V1
# INFILE FOR ALL SPECIES 
enfile="functional/8sp.enrichment.groups.tsv"
## PLOT TITLES FOR COMBINED ENRICHMENT PLOTS
header_a=c("A", "A", "C", "E", "G", "I", "K", "M")

# Loop over all species 
for (i in 1:length(sptib$V1)){
  show(paste("Running pipeline for:",sptib$V1[i]))
  
  # READ IN DATA 
  minmaxfile=paste("functional/",sptib$V3[i],".RegionSampling.minmax.50perc.100rep.96CI.txt", sep="")
  minmaxtib<-minmaxfile %>% read.table(header=TRUE) %>% as_tibble()
  # Read in the enrichment data
  ENRICH <- enfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
    filter(Species==sptib$V1[i]) %>% inner_join(minmaxtib) %>% 
    mutate(PrintName=case_when(Class=="lncrna" ~ "lncRNA",  
                               Class=="promoter" ~ "1 kb upstream\nfrom TSS\n(putative promoters)",
                               Class=="introns" ~ "Intronic",
                               Class=="intergenic" ~ "Intergenic",
                               Class=="UTR5" ~ "5'UTRs",
                               Class=="UTR3" ~ "3'UTRs",
                               TRUE ~ Class)) %>%
    mutate(Group=case_when(Group=="macro" ~ "Macro",
                           Group=="micro" ~ "Micro",
                           Group=="dot" ~ "Dot",
                           TRUE ~ Group)) 
  ENRICH$nonB <- factor(ENRICH$nonB, levels=c("ALL","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))
  ENRICH$PrintName <- factor(ENRICH$PrintName, levels=c("1 kb upstream\nfrom TSS\n(putative promoters)","5'UTRs","CDS", "3'UTRs", "Intronic", "lncRNA", "Intergenic"))
  
  max_macro<-ceiling(max(ENRICH$Max[ENRICH$Group=="Macro"]))+1
  max_micro<-ceiling(max(ENRICH$Max[ENRICH$Group=="Micro"]))+1
  max_dot<-ceiling(max(ENRICH$Max[ENRICH$Group=="Dot"]))+1
  
  
  p1 <- make_enrich_plot(ENRICH%>%filter(Group=="Macro"), max_macro, "Macro")  + 
    theme(axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "top",  
          legend.justification = "right",
          plot.title=element_text(size=18, face="bold", hjust=0)) + 
    guides(fill = guide_legend(nrow = 1)) #+ggtitle("A")
  p2 <- make_enrich_plot(ENRICH%>%filter(Group=="Micro"), max_micro, "Micro") + theme(axis.text.x = element_blank())
  p3 <- make_enrich_plot(ENRICH%>%filter(Group=="Dot"), max_dot, "Dot") + theme(axis.title.y = element_blank())
  
  # Combine plots with patchwork with different heights
  pa<-(p1 / p2 / p3) + plot_layout(heights = c(max_macro, max_micro, max_dot))+
    plot_annotation(
      title = header_a[i],
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0)
      )
    )
  pa
  
  plots[[i]] <- pa
}

# PLOT EACH SPECIES 
for (nm in names(plots)) {
  ggsave(
    filename = paste0("plots/", nm, "_NonB_FunctionalEnrichment.png"),
    plot = plots[[nm]],
    scale = 1,dpi = 300,limitsize = TRUE,width=9,height=9
  )
}


################################################################################
# FIGURE 2B: FUNCTION FOR PLOTTING G4 ON TEMPLATE VS CODING STRAND 
################################################################################

make_bplot_with_stats <- function(data, data2, y_max, classif = NULL) {
  # Summary of the data
  data_summary <- data %>%
    group_by(PrintName, Strand) %>%
    summarise(mean_cov = mean(Coverage, na.rm = TRUE), .groups = "drop")
  # Paired wilcoxons test
  wilcox_results <- data %>%
    # first collapse duplicates per class/strand/PrintName
    group_by(PrintName, Gene, Strand) %>%
    summarise(mean_cov = mean(Coverage, na.rm = TRUE), .groups = "drop") %>%
    # now pivot to wide so each gene has Coding & Template in columns
    pivot_wider(names_from = Strand, values_from = mean_cov) %>%
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
  p<-ggplot(data_summary, aes(x = PrintName, y = Coverage*100,
                              fill = Strand, color = Strand)) +
    geom_bar(data=data2, stat = "identity", position = position_dodge(width = 0.9), alpha = 0.6) +
    scale_fill_manual(values = c(vcolors[6],"skyblue")) +
    scale_color_manual(values = c(vcolors[6],"blue")) +
    # Add significance text
    geom_text(
      data = wilcox_results,
      aes(x = PrintName, y = y_max -1.2, 
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT FOR ALL SPECIES 

g4_plots <- vector("list", length(sptib$V1))
names(g4_plots) <- sptib$V1
# INFILEs FOR ALL SPECIES 
strfile="functional/8sp.G4_strand.perGene.tsv"
sumfile="functional/8sp.G4_strand.groups.tsv"
# HEADER FOR COMBINED PLOTS
header_b=c("B", "B", "D", "F", "H", "J", "L", "N")



# Get max values using a function 
get_max <- function(data, group) {
  ceiling(100*(data %>% filter(Group==group) %>% group_by(Class, Strand) %>%
      summarize(max=max(Coverage),na.rm = TRUE, .groups="drop") %>%
        summarize(max=max(max)))$max)+1.5
}

# Loop over all species 
for (i in 1:length(sptib$V1)){
  show(paste("Running G4 Coding vs Template pipeline for:",sptib$V1[i]))
  
  # READ IN DATA 
  G4SUM <- sumfile %>% read.table(header=TRUE) %>% as_tibble() %>%
    filter(Species==sptib$V1[i]) %>%
    mutate(PrintName=case_when(Class=="lncrna" ~ "lncRNA",  
                               Class=="promoter" ~ "1 kb upstream\nfrom TSS",
                               Class=="introns" ~ "Intronic",
                               Class=="UTR5" ~ "5'UTRs",
                               Class=="UTR3" ~ "3'UTRs",
                               Class=="intergenic" ~ "",
                               TRUE ~ Class)) %>%
    mutate(Group=case_when(Group=="macro" ~ "Macro",
                           Group=="micro" ~ "Micro",
                           Group=="dot" ~ "Dot",
                           TRUE ~ Group)) 
  G4SUM$PrintName <- factor(G4SUM$PrintName, levels=c("1 kb upstream\nfrom TSS","5'UTRs","CDS", "3'UTRs", "Intronic", "lncRNA",""))
  
  G4STRAND <- strfile %>% read.table(header=TRUE) %>% as_tibble() %>%
    filter(Species==sptib$V1[i]) %>%
    add_row(Class="intergenic", Group="Macro", Strand="Coding", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
    add_row(Class="intergenic", Group="Micro", Strand="Coding", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
    add_row(Class="intergenic", Group="Dot", Strand="Coding", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
    add_row(Class="intergenic", Group="Macro", Strand="Template", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
    add_row(Class="intergenic", Group="Micro", Strand="Template", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
    add_row(Class="intergenic", Group="Dot", Strand="Template", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
    mutate(PrintName=case_when(Class=="lncrna" ~ "lncRNA",  
                               Class=="promoter" ~ "1 kb upstream\nfrom TSS",
                               Class=="introns" ~ "Intronic",
                               Class=="UTR5" ~ "5'UTRs",
                               Class=="UTR3" ~ "3'UTRs",
                               Class=="intergenic" ~ "",
                               TRUE ~ Class)) %>%
    mutate(Group=case_when(Group=="macro" ~ "Macro",
                           Group=="micro" ~ "Micro",
                           Group=="dot" ~ "Dot",
                           TRUE ~ Group)) 
  G4STRAND$PrintName <- factor(G4STRAND$PrintName, levels=c("1 kb upstream\nfrom TSS","5'UTRs","CDS", "3'UTRs", "Intronic", "lncRNA",""))
  
  g4max_macro<-get_max(G4SUM, "Macro")
  g4max_micro<-get_max(G4SUM, "Micro")
  g4max_dot<-get_max(G4SUM, "Dot")
  
  
  # Create the individual plots
  pb1 <- make_bplot_with_stats(G4STRAND%>%filter(Group=="Macro"), G4SUM%>%filter(Group=="Macro"), g4max_macro, "Macro")  + 
    theme(axis.text.x = element_blank(),
          legend.position = "top",  
          legend.justification = "right",
          axis.title.y = element_blank(),
          plot.title=element_text(size=18, face="bold", hjust=0)) + 
    guides(fill = guide_legend(nrow = 1)) #+ggtitle("B")
  pb2 <- make_bplot_with_stats(G4STRAND%>%filter(Group=="Micro"), G4SUM%>%filter(Group=="Micro"), g4max_micro, "Micro") + theme(axis.text.x = element_blank())
  pb3 <- make_bplot_with_stats(G4STRAND%>%filter(Group=="Dot"), G4SUM%>%filter(Group=="Dot"), g4max_dot, "Dot") + theme(axis.title.y = element_blank())
  
  # Combine plots with patchwork with different heights
  pb<-(pb1 / pb2 / pb3) + 
    plot_layout(heights = c(g4max_macro, g4max_micro, g4max_dot))+
    plot_annotation(
      title = header_b[i],
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0)
      )
    )
  pb
  
  g4_plots[[i]] <- pb
}


# PLOT EACH SPECIES 
for (nm in names(g4_plots)) {
  ggsave(
    filename = paste0("plots/", nm, "_NonB_G4_strand_coverage.png"),
    plot = g4_plots[[nm]],
    scale = 1,dpi = 300,limitsize = TRUE,width=9,height=9
  )
}

################################################################################
# CODE FOR FIGURE 2C, METHYLATION IN AND OUTSIDE G4s 
################################################################################

# Make function to calculate p-values 
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
    scale_fill_manual(values = c(vcolors[6],"gray90")) +
    scale_color_manual(values = c(vcolors[6],"gray30")) +
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the data 

# Infiles 
methfile="methylation/bTaeGut7v0.4_MT_rDNA.allCpG.group.txt"
intergenfile="methylation/bTaeGut7v0.4_MT_rDNA.allCpG.intergenic.group.txt"

# Make Tibbles
methtib<-methfile %>% read.table(header=TRUE) %>% as_tibble() 
intergentib<-intergenfile %>% read.table(header=TRUE) %>% as_tibble() 
#sumtib<-sumfile %>% read.table(header=TRUE) %>% as_tibble() 

# Summarize gene regions with median and mean per transcript 
DATAC1 <-methfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  group_by(Group,Class,Type,Trx) %>% 
  summarize(median=median(Score), mean=mean(Score)) %>% 
  mutate(PrintName=case_when(Class=="lncrna" ~ "lncRNA",  
                             Class=="promoter" ~ "1 kb upstream\nfrom TSS",
                             Class=="introns" ~ "Intronic",
                             Class=="intergenic" ~ "Intergenic",
                             Class=="UTR5" ~ "5'UTRs",
                             Class=="UTR3" ~ "3'UTRs",
                             TRUE ~ Class)) %>%
  mutate(Group=case_when(Group=="macro" ~ "Macro",
                         Group=="micro" ~ "Micro",
                         Group=="dot" ~ "Dot",
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
                         Group=="dot" ~ "Dot",
                         TRUE ~ Group)) %>% 
  mutate(Region=case_when(Type=="ignorant" ~ "inside G4",
                          Type=="background" ~ "outside G4",
                          TRUE ~ NA)) %>% mutate(PrintName="Intergenic", median=Score)
DATAC2$Group <- factor(DATAC2$Group, levels=c("Macro","Micro","Dot"))



# Calculate pvalues 
pvals1<-get_pval(DATAC1)
pvals2<-get_pval(DATAC2)

# Plot each class panel and combine 
pc1<-make_densityplot(DATAC1%>%filter(Class=="promoter"), pvals1%>%filter(Class=="promoter"),FALSE, FALSE, TRUE)+ 
  theme(plot.title=element_text(size=18, face="bold")) #+ggtitle("C")

pc2<-make_densityplot(DATAC1%>%filter(Class=="UTR5"), pvals1%>%filter(Class=="UTR5"), FALSE, FALSE, FALSE)
pc3<-make_densityplot(DATAC1%>%filter(Class=="CDS"), pvals1%>%filter(Class=="CDS"), FALSE, FALSE, FALSE)
pc4<-make_densityplot(DATAC1%>%filter(Class=="UTR3"), pvals1%>%filter(Class=="UTR3"), TRUE, TRUE, FALSE)
pc5<-make_densityplot(DATAC1%>%filter(Class=="introns"), pvals1%>%filter(Class=="introns"), FALSE, FALSE, FALSE)
pc6<-make_densityplot(DATAC1%>%filter(Class=="lncrna"), pvals1%>%filter(Class=="lncrna"), FALSE, FALSE, FALSE)
pc7<-make_densityplot(DATAC2, pvals2, FALSE, FALSE, FALSE, TRUE)

pc<-pc1+pc2+pc3+pc4+pc5+pc6+pc7 + plot_layout(widths=c(1,1,1,1,1,1,1)) +
  plot_annotation(
    title = "C",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0)
    )
  )
pc



################################################################################
# COMBINE THE PANELS INTO FINAL FIGURES
################################################################################

# FIGURE 2
pa<-plots[["zebra_finch"]] 
pb<-g4_plots[["zebra_finch"]] 

p_final<- wrap_elements(pa)/ wrap_elements(pb) / wrap_elements(pc)
p_final
outfile=paste(plotdir,"Fig3_functional.png", sep="")
ggsave(outfile,plot = p_final,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=15)
outfile=paste(plotdir,"Fig3_functional.svg", sep="")
ggsave(outfile,plot = p_final,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=15)

# SUPPLEMENTARY, ONE SIDE PER SPECIES 
for (nm in names(g4_plots)) {
  if(nm!="zebra_finch") {
    pa<-plots[[nm]] 
    pb<-g4_plots[[nm]] 
    p2_final<- wrap_elements(pa)/ wrap_elements(pb)
    ggsave(
      filename = paste0("plots/", nm, "_Functional_combined.png"),
      plot = p2_final,
      scale = 1,dpi = 300,limitsize = TRUE,width=9,height=13
    )
  }
}


