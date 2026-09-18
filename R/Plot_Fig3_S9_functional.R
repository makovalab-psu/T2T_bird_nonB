################################################################################
### R code for plotting non-B motif enrichment in functional regions in the 
### Zebra finch T2T genome.
# written by Linnéa Smeds Sept 15, 2025.

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

# SPECIES LIST 
spfile="helpfiles/species_list.txt"
sptib<-spfile %>% read.table(header=FALSE) %>% as_tibble()

# Set colors 
vcolors=viridis(7)
type_col=c("#440154","#9460a1", "#e3cfe8")

################################################################################
# FUNCTION TO MAKE ENRICHMENT PLOT FOR PANEL 3A AND SUPPLEMENTARY FIGURES
################################################################################

# Help function, make triangle to mark chromosome category average
make_left_triangle <- function(data, x_col = "base_x", y_col = "Enrichment",
                               size_x = 0.1, size_y = NULL, id_col = NULL) {
  # default vertical size scaled relative to the y-range actually present
  if (is.null(size_y)) {
    size_y <- diff(range(data[[y_col]], na.rm = TRUE)) * 0.03
    if (size_y == 0 || is.na(size_y)) size_y <- 0.3
  }
  if (is.null(id_col)) {
    data$.tri_id <- seq_len(nrow(data))
  } else {
    data$.tri_id <- data[[id_col]]
  }
  
  # For each point, create 3 vertices: tip pointing LEFT, base pointing RIGHT
  purrr::pmap_dfr(data, function(...) {
    row <- list(...)
    cx <- row[[x_col]]
    cy <- row[[y_col]]
    tibble::tibble(
      x = c(cx - size_x, cx + size_x, cx + size_x),  # tip, top-right, bottom-right
      y = c(cy,          cy + size_y, cy - size_y),
      .tri_id = row[[".tri_id"]],
      NonB = row[["NonB"]]
    )
  })
}


# Helper: for each bar (PrintName x nonB), flag whether its [Min, Max] error-bar
# range fails to overlap the matching baseline value for that nonB type (from
# data2). Returns `data` with an added `.sig` logical column and a `.dodge_x`
# numeric x-position matching where position_dodge() would have placed that bar,
# so the asterisk lands directly above the correct bar even with multiple nonB
# categories dodged side-by-side within each PrintName group.
flag_baseline_sig <- function(data, data2, x_col = "PrintName",
                              bar_key = "nonB", base_key = "NonB",
                              base_val = "Enrichment", dodge_width = 0.9) {
  base_lookup <- data2 %>%
    dplyr::select(dplyr::all_of(c(base_key, base_val))) %>%
    dplyr::rename(.base_val = dplyr::all_of(base_val)) %>%
    dplyr::rename(!!bar_key := dplyr::all_of(base_key))
  
  # Determine dodge x-offsets: same ordering ggplot's position_dodge() uses.
  # IMPORTANT: position_dodge() dodges based on the groups actually PRESENT
  # in this panel's data, not the full set of factor levels declared globally
  # (e.g. nonB has 8 levels overall, but a given Group=="Macro"/"Micro"/"Dot"
  # subset may only contain a handful of them). Using the global level count
  # here would divide the dodge width by the wrong denominator and misalign
  # offsets whenever a panel is missing some categories -- so we drop unused
  # levels first to match what ggplot itself actually dodges.
  used_levels <- levels(droplevels(data[[bar_key]]))
  n_grp <- length(used_levels)
  grp_idx <- as.numeric(factor(data[[bar_key]], levels = used_levels))
  offset <- (grp_idx - (n_grp + 1) / 2) * (dodge_width / n_grp)
  
  data %>%
    dplyr::left_join(base_lookup, by = bar_key) %>%
    dplyr::mutate(
      .sig = .data$Min > .data$.base_val,
      .dodge_x = as.numeric(.data[[x_col]]) + offset
    )
}



# Function to make consistent bar plots
make_enrich_plot <- function(data, data2, y_max = NULL, classif = NULL) {

  n_levels <- length(levels(data$PrintName))
  base_x <- n_levels + 0.7
  
  set.seed(1)
  data2 <- data2 %>%
    dplyr::mutate(base_x = base_x + runif(dplyr::n(), -0.08, 0.08))
  tri_data <- make_left_triangle(data2, x_col = "base_x", y_col = "Enrichment")
  
  # Flag bars whose error range does NOT overlap their matching baseline value
  sig_data <- flag_baseline_sig(data, data2) %>%
    dplyr::filter(.sig)
  
  p<-ggplot(data, aes(x=PrintName, y=Enrichment_gw, fill=nonB, color=nonB)) + 
    geom_bar(stat="identity", position=position_dodge(), alpha=0.6) +
    geom_errorbar(stat="identity", position=position_dodge(), aes(x=PrintName, ymin=Min, ymax=Max), colour="grey40")+
    scale_fill_manual(values=vcolors) +
    scale_color_manual(values=vcolors) +
 #   geom_hline(yintercept=1, linetype="dashed", color = "red") +
 #   geom_point(data=data2, aes(x=base_x, y=Enrichment, color=NonB, fill=NonB), alpha=0.8, shape = "\u25C0", inherit.aes = FALSE,
#               size=4, position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    #scale_x_discrete(expand = c(0.08, 0.08)) +
    geom_polygon(data = tri_data, aes(x = x, y = y, group = .tri_id, fill = NonB, color=NonB),
                 inherit.aes = FALSE, alpha=0.8, linewidth = 0.3) +
    
    geom_text(data = sig_data, aes(x = .dodge_x, y = Max, label = "*"),
              inherit.aes = FALSE, size = 6, vjust = 0.3, nudge_y = {
                if (!is.null(y_max)) y_max * 0.015 else diff(range(data$Max, na.rm = TRUE)) * 0.03
              }) +
    
    scale_x_discrete(expand = expansion(add = c(0.08, 0.8))) +
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
      axis.ticks.y=element_line(),
     # plot.margin = margin(t = 5.5, r = 12, b = 5.5, l = 5.5) # extra right margin for points
    ) +
    # Allow the points to render even if they sit right at/near the panel edge
    coord_cartesian(clip = "off", ylim = if (!is.null(y_max)) c(0, y_max) else NULL)
  if (!is.null(y_max)) {
    p <- p + 
      scale_y_continuous(expand = c(0.0, 0.0), breaks = seq(0, y_max, by = 3))
  }
  if (!is.null(classif)) {
    p <- p +
      annotate("text", x = Inf, y = y_max, label = classif,
               hjust = 1.2, vjust = 1.6, size = 4 , alpha = 0.65)  # subtle label
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
basefile="coverage/9sp.group_enrichment.tsv"

## PLOT TITLES FOR COMBINED ENRICHMENT PLOTS
header_a=c("A", "", "E", "G", "A", "C", "K", "I")

# Loop over all species except golden phesant
for (i in 1:8){
  show(paste("Running pipeline for:",sptib$V1[i]))
  
  # READ IN DATA 
  minmaxfile=paste("functional/",sptib$V3[i],".RegionSampling.minmax.50perc.100rep.96CI.txt", sep="")
  minmaxtib<-minmaxfile %>% read.table(header=TRUE) %>% as_tibble()
  # Read in the enrichment data
  ENRICH <- enfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
    filter(Species==sptib$V1[i]) %>% inner_join(minmaxtib) %>% 
    filter(Class!="lncrna") %>% 
    mutate(PrintName=case_when(#Class=="lncrna" ~ "lncRNA",  
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
  ENRICH$PrintName <- factor(ENRICH$PrintName, levels=c("1 kb upstream\nfrom TSS\n(putative promoters)","5'UTRs","CDS", "3'UTRs", "Intronic", "Intergenic")) #"lncRNA", 
  
  BASE <- basefile %>% read.table(header=TRUE) %>% as_tibble() %>% 
    filter(Species==sptib$V1[i]) %>% filter(NonB!="All")
  max_macro<-ceiling(max(ENRICH$Max[ENRICH$Group=="Macro"]))+1
  max_micro<-ceiling(max(ENRICH$Max[ENRICH$Group=="Micro"]))+1
  max_dot<-ceiling(max(ENRICH$Max[ENRICH$Group=="Dot"]))+1
  
  if(sptib$V1[i]=="chicken") {
    p1 <- make_enrich_plot(ENRICH%>%filter(Group=="Macro"), BASE%>%filter(Group=="macro"), max_macro, "Macro")  + 
      theme(axis.text.x = element_blank(),
            axis.title.y = element_blank())
  }
  else {
    p1 <- make_enrich_plot(ENRICH%>%filter(Group=="Macro"), BASE%>%filter(Group=="macro"), max_macro, "Macro")  + 
      theme(axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "top",  
            legend.justification = "right") + 
      guides(fill = guide_legend(nrow = 1)) #+ggtitle("A")
  }
  p2 <- make_enrich_plot(ENRICH%>%filter(Group=="Micro"), BASE%>%filter(Group=="micro"), max_micro, "Micro") + theme(axis.text.x = element_blank())
  p3 <- make_enrich_plot(ENRICH%>%filter(Group=="Dot"), BASE%>%filter(Group=="dot"), max_dot, "Dot") + theme(axis.title.y = element_blank())
  
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
    filename = paste0("plots/", nm, "_NonB_FunctionalEnrichment_nolncRNA.png"),
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
  
  bar_heights <- data2 %>%
    group_by(PrintName) %>%
    summarise(bar_top = max(Coverage, na.rm = TRUE) * 100, .groups = "drop")
  wilcox_results <- wilcox_results %>%
    left_join(bar_heights, by = "PrintName") %>%
    mutate(text_y = bar_top + (y_max * 0.03))  # small offset above the taller bar
  
  #And plot 
  p<-ggplot(data_summary, aes(x = PrintName, y = Coverage*100,
                              fill = Strand, color = Strand)) +
    geom_bar(data=data2, stat = "identity", position = position_dodge(width = 0.9), alpha = 0.6) +
    scale_fill_manual(values = c(vcolors[6],"skyblue")) +
    scale_color_manual(values = c(vcolors[6],"blue")) +
    # Add significance text
    geom_text(
      data = wilcox_results,
      aes(x = PrintName, y = text_y, 
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
      axis.text.y = element_text(size=10),
      axis.text.x = element_text(size=8, angle=45, hjust=1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      title = element_blank(),
      axis.ticks.y=element_line(),
      axis.title.y=element_text(size=12),
   #   plot.margin=margin(t=0,r=0,b=0,l=0),
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
               size = 4 , alpha = 0.65)  # subtle label
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
header_b=c("B", "", "F", "H", "B", "D", "L", "J")



# Get max values using a function 
get_max <- function(data, group) {
  ceiling(100*(data %>% filter(Group==group) %>% group_by(Class, Strand) %>%
      summarize(max=max(Coverage),na.rm = TRUE, .groups="drop") %>%
        summarize(max=max(max)))$max)+1.5
}

# Loop over all species 
for (i in 1:8){
  show(paste("Running G4 Coding vs Template pipeline for:",sptib$V1[i]))
  
  # READ IN DATA 
  G4SUM <- sumfile %>% read.table(header=TRUE) %>% as_tibble() %>%
    filter(Species==sptib$V1[i]) %>%
    filter(Class!="lncrna" & Class!="intergenic") %>%
    mutate(PrintName=case_when( #Class=="lncrna" ~ "lncRNA",  
                               Class=="promoter" ~ "Put. prom.",
                               Class=="introns" ~ "Intronic",
                               Class=="UTR5" ~ "5'UTRs",
                               Class=="UTR3" ~ "3'UTRs",
                               TRUE ~ Class)) %>%
    mutate(Group=case_when(Group=="macro" ~ "Macro",
                           Group=="micro" ~ "Micro",
                           Group=="dot" ~ "Dot",
                           TRUE ~ Group)) 
  G4SUM$PrintName <- factor(G4SUM$PrintName, levels=c("Put. prom.","5'UTRs","CDS", "3'UTRs", "Intronic")) #"lncRNA",
  
  G4STRAND <- strfile %>% read.table(header=TRUE) %>% as_tibble() %>%
    filter(Species==sptib$V1[i]) %>%
    filter(Class!="lncrna") %>% 
 #   add_row(Class="intergenic", Group="Macro", Strand="Coding", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
#    add_row(Class="intergenic", Group="Micro", Strand="Coding", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
#    add_row(Class="intergenic", Group="Dot", Strand="Coding", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
#    add_row(Class="intergenic", Group="Macro", Strand="Template", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
#    add_row(Class="intergenic", Group="Micro", Strand="Template", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
#    add_row(Class="intergenic", Group="Dot", Strand="Template", Coverage=0, Gene="dummy", GenLen=0, G4Len=0) %>% 
    mutate(PrintName=case_when(#Class=="lncrna" ~ "lncRNA",  
                               Class=="promoter" ~ "Put. prom.",
                               Class=="introns" ~ "Intronic",
                               Class=="UTR5" ~ "5'UTRs",
                               Class=="UTR3" ~ "3'UTRs",
                               TRUE ~ Class)) %>%
    mutate(Group=case_when(Group=="macro" ~ "Macro",
                           Group=="micro" ~ "Micro",
                           Group=="dot" ~ "Dot",
                           TRUE ~ Group)) 
  G4STRAND$PrintName <- factor(G4STRAND$PrintName, levels=c("Put. prom.","5'UTRs","CDS", "3'UTRs", "Intronic")) #"lncRNA",
  
  g4max_macro<-get_max(G4SUM, "Macro")
  g4max_micro<-get_max(G4SUM, "Micro")
  g4max_dot<-get_max(G4SUM, "Dot")
  
  
  # Create the individual plots
  if(sptib$V1[i]=="chicken") {
    pb1 <- make_bplot_with_stats(G4STRAND%>%filter(Group=="Macro"), G4SUM%>%filter(Group=="Macro"), g4max_macro, "Macro")  + 
      theme(axis.text.x = element_blank(),
            axis.title.y = element_blank())
  }
  else{
    pb1 <- make_bplot_with_stats(G4STRAND%>%filter(Group=="Macro"), G4SUM%>%filter(Group=="Macro"), g4max_macro, "Macro")  + 
      theme(axis.text.x = element_blank(),
            legend.position = "top",  
            legend.justification = "right",
            axis.title.y = element_blank(),
            plot.title=element_text(size=18, face="bold", hjust=0)) + 
      guides(fill = guide_legend(nrow = 1)) #+ggtitle("B")
  }
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
    filename = paste0("plots/", nm, "_NonB_G4_strand_coverage_nolncrna.png"),
    plot = g4_plots[[nm]],
    scale = 1,dpi = 300,limitsize = TRUE,width=9,height=9
  )
}

################################################################################
################################################################################
# COMBINE THE PANELS INTO FINAL FIGURES
################################################################################

# FIGURE 3
pa1<-plots[["zebra_finch"]] 
pa2<-plots[["chicken"]] 
pb1<-g4_plots[["zebra_finch"]] 
pb2<-g4_plots[["chicken"]] 


pleft<- wrap_elements(pa1) / wrap_elements(pa2)
pright<- wrap_elements(pb1) / wrap_elements(pb2) 
p_final<- (wrap_elements(pleft) + wrap_elements(pright)) + plot_layout(ncol=2, widths=c(2,1))
p_final

outfile=paste(plotdir,"Fig3_functional_withchicken.png", sep="")
ggsave(outfile,plot = p_final,scale = 1,dpi = 300,limitsize = TRUE,width=12,height=12)
outfile=paste(plotdir,"Fig3_functional_withchicken.svg", sep="")
ggsave(outfile,plot = p_final,scale = 1,dpi = 300,limitsize = TRUE,width=12,height=12)


# ------------------------------------------------------------------------
# SUPPLEMENTARY FIG9, ONE SIDE PER SPECIES 
for (nm in names(g4_plots)) {
  if(nm!="zebra_finch" & nm!="chicken") {
    pa<-plots[[nm]] 
    pb<-g4_plots[[nm]] 
    p2_final<- wrap_elements(pa) + wrap_elements(pb) +  plot_layout(ncol=2, widths=c(2,1))
    ggsave(
      filename = paste0("plots/FigS9_", nm, "_Functional_combined_nolncrna.png"),
      plot = p2_final,
      scale = 1,dpi = 300,limitsize = TRUE,width=9,height=7
    )
  }
}

################################################################################
# FOR POSTER
outfile=paste(plotdir,"Fig3A_for_poster.png", sep="")
ggsave(outfile,plot = pa,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=6)

# Plot combined C (not divided into chromosome groups) for poster:

# Make function to calculate p-values 
get_pval <- function(data) {
  pvals <- data %>%
    group_by(PrintName,Class) %>%
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
    facet_nested(~PrintName, scales = "free_y") +  # single-column facet
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
      legend.justification = "right",
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


# Summarize gene regions with median and mean per transcript 
DATAC1 <-methfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  group_by(Group,Class,Type,Trx) %>% 
  summarize(median=median(Score), mean=mean(Score)) %>% 
  filter(Class!="lncrna") %>% 
  mutate(PrintName=case_when(#Class=="lncrna" ~ "lncRNA",  
    Class=="promoter" ~ "1 kb upstream\nfrom TSS",
    Class=="introns" ~ "Intronic",
    Class=="intergenic" ~ "Intergenic",
    Class=="UTR5" ~ "5'UTRs",
    Class=="UTR3" ~ "3'UTRs",
    TRUE ~ Class)) %>%
  mutate(PrintType=case_when(Type=="coding" ~ "Coding",
                             Type=="template" ~ "Template",
                             Type=="background" ~ "Background",
                             TRUE ~ Group)) %>% 
  mutate(Region=case_when(Type=="coding" ~ "inside G4",
                          Type=="template" ~ "inside G4",
                          Type=="background" ~ "outside G4",
                          TRUE ~ NA))
DATAC1$PrintType <- factor(DATAC1$PrintType, levels=c("Background","Coding","Template"))
DATAC1$PrintName <- factor(DATAC1$PrintName, levels=c("1 kb upstream\nfrom TSS","5'UTRs","CDS", "3'UTRs", "Intronic", "Intergenic")) #"lncRNA",
DATAC1$Region <- factor(DATAC1$Region, levels=c("inside G4","outside G4"))

# For intergenic we use all sites 
DATAC2<-intergentib %>%
  mutate(Region=case_when(Type=="ignorant" ~ "inside G4",
                          Type=="background" ~ "outside G4",
                          TRUE ~ NA)) %>% mutate(PrintName="Intergenic", median=Score)


# Calculate pvalues 
pvals1<-get_pval(DATAC1)
pvals2<-get_pval(DATAC2)

# Plot each class panel and combine 
pc1<-make_densityplot(DATAC1%>%filter(Class=="promoter"), pvals1%>%filter(Class=="promoter") ,FALSE, FALSE, TRUE)+ 
  theme(plot.title=element_blank()) 
pc2<-make_densityplot(DATAC1%>%filter(Class=="UTR5"), pvals1%>%filter(Class=="UTR5"), FALSE, FALSE, FALSE)
pc3<-make_densityplot(DATAC1%>%filter(Class=="CDS"), pvals1%>%filter(Class=="CDS"), FALSE, FALSE, FALSE)
pc4<-make_densityplot(DATAC1%>%filter(Class=="UTR3"), pvals1%>%filter(Class=="UTR3"), TRUE, TRUE, FALSE)
pc5<-make_densityplot(DATAC1%>%filter(Class=="introns"), pvals1%>%filter(Class=="introns"), FALSE, FALSE, FALSE)
#pc6<-make_densityplot(DATAC1%>%filter(Class=="lncrna"), pvals1%>%filter(Class=="lncrna"), FALSE, FALSE, FALSE)
pc7<-make_densityplot(DATAC2, pvals2, FALSE, FALSE, FALSE, TRUE)

pc<-pc1+pc2+pc3+pc4+pc5+pc7 + plot_layout(widths=c(1,1,1,1,1,1))
pc




ggsave(outfile,plot = pc,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=3.5)
outfile=paste(plotdir,"Fig3C_for_poster.png", sep="")





