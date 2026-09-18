################################################################################
### R code for plotting methylation, CD and G4-seq
# written by Linnéa Smeds 2026, code streamlined with Claude Sonnet 5 (Anthropic).

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(viridis)
require(patchwork)
require(ggtext)
library(ggplot2)
library(data.table)
library(jsonlite)

plotdir="plots/"

# Set colors 
vcolors=viridis(7)
type_col=c("#440154","#9460a1", "#e3cfe8")

################################################################################
################################################################################
################################################################################
# CODE FOR FIGURE 5A, METHYLATION IN AND OUTSIDE G4s 

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

# Single-call density plot function for the FULL combined facet grid
# (Group x PrintName, all 18 panels at once -> uniform sizing, one set
# of margins/axes/legend instead of 6 stitched-together sub-plots).
#
# KEY FIX: facet_grid2(scales="free_y", independent="y") gives EVERY
# individual panel its own y-axis range (unlike facet_grid/facet_nested,
# where "free_y" can only vary *by row*, not per-panel). This is what
# lets us keep readable y-scales for all 18 panels without having to
# build them as 6 separate stitched-together plots.
make_densityplot_full <- function(data, sigdata) {
  p <- ggplot(data, aes(x = median/100, fill = Region, color = Region)) + 
    geom_density(alpha = 0.6) +
    facet_grid2(Group ~ PrintName, scales = "free_y", independent = "y",
                axes = "x", remove_labels = "x") +
    scale_fill_manual(values = c(vcolors[6],"gray90")) +
    scale_color_manual(values = c(vcolors[6],"gray30")) +
    geom_text(
      data = sigdata,
      aes(x = x, y = y, label = signif),
      inherit.aes = FALSE,
      vjust = 1.5,
      size = 4
    ) +
    labs(x="Methylation level", y="Density") +
    scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      panel.spacing = unit(3, "pt"),           # tight spacing between facets
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      plot.title = element_text(size=18, face="bold"),
      strip.text.x = element_text(size = 10),
      strip.text.y = element_text(size = 12),
      strip.background = element_blank(),
      legend.key = element_rect(fill = NA, color = NA),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(size = 8),     # per-panel y text now visible everywhere -> keep small
      axis.ticks.length = unit(3, "pt"),
      plot.margin = margin(2, 2, 2, 2)          # trim default 5.5pt margins
    )
  return(p)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the data 

methfile="methylation/bTaeGut7v0.4_MT_rDNA.allCpG.group.txt"
intergenfile="methylation/bTaeGut7v0.4_MT_rDNA.allCpG.intergenic.group.txt"

methtib<-methfile %>% read.table(header=TRUE) %>% as_tibble() 
intergentib<-intergenfile %>% read.table(header=TRUE) %>% as_tibble() 

# Summarize gene regions with median and mean per transcript 
DATAA1 <-methfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  group_by(Group,Class,Type,Trx) %>% 
  summarize(median=median(Score), mean=mean(Score)) %>% 
  filter(Class!="lncrna") %>% 
  mutate(PrintName=case_when(
    Class=="promoter" ~ "Put. prom.",
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
DATAA1$Group <- factor(DATAA1$Group, levels=c("Macro","Micro","Dot"))
DATAA1$PrintType <- factor(DATAA1$PrintType, levels=c("Background","Coding","Template"))
DATAA1$Region <- factor(DATAA1$Region, levels=c("inside G4","outside G4"))

# For intergenic we use all sites 
DATAA2<-intergentib %>%
  mutate(Group=case_when(Group=="macro" ~ "Macro",
                         Group=="micro" ~ "Micro",
                         Group=="dot" ~ "Dot",
                         TRUE ~ Group)) %>% 
  mutate(Region=case_when(Type=="ignorant" ~ "inside G4",
                          Type=="background" ~ "outside G4",
                          TRUE ~ NA)) %>% mutate(PrintName="Intergenic", median=Score)
DATAA2$Group <- factor(DATAA2$Group, levels=c("Macro","Micro","Dot"))

# Calculate pvalues 
pvals1<-get_pval(DATAA1)
pvals2<-get_pval(DATAA2)

# Combine everything into a tidy dataset
levels_order <- c("Put. prom.","5'UTRs","CDS","3'UTRs","Intronic","Intergenic")

DATAA_all <- bind_rows(
  DATAA1 %>% filter(Class %in% c("promoter","UTR5","CDS","UTR3","introns")) %>%
    select(Group, Class, PrintName, median, Region),
  DATAA2 %>% select(Group, Class, PrintName, median, Region)
) %>%
  mutate(PrintName = factor(PrintName, levels = levels_order))

pvals_all <- bind_rows(
  pvals1 %>% filter(Class %in% c("promoter","UTR5","CDS","UTR3","introns")),
  pvals2
) %>%
  mutate(PrintName = factor(PrintName, levels = levels_order))

pA <- make_densityplot_full(DATAA_all, pvals_all) +
  plot_annotation(
    title = "A",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0))
  )
pA

################################################################################
################################################################################
################################################################################
# CODE FOR FIG 5B, CD SPECTRA, Written by Jacob Sieg


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
    labs(x = "Wavelength (nm)", y = "Molar ellipticity") +   # real labels -> lets patchwork collect them
    ggtitle(panelname)
  return(P_CD)
}
####Plot and save data####
#"Tgut368A" "CR1 LINE", "CR1 A", "CR1 B"
pB1 = plot_data(x = "Tgut368A", "") + theme(legend.position = "none")
pB2 = plot_data(x = "CR1_LINE", "") + theme(legend.position = "none")
pB3 = plot_data(x = "CR1_A", "") + theme(legend.position = "none")
pB4 = plot_data(x = "CR1_B", "") + theme(legend.position = "none")

# One row, shared axis titles + collected axis text 
pB_row <- (pB1 + pB2 + pB3 + pB4) +
  plot_layout(nrow = 1, axis_titles = "collect", axes = "collect")

# Get legend
legend_plot <- get_legend(
  plot_data(x = "CR1_A", "") + theme(legend.position = "bottom")
)

# wrap plots and legend
pB <- (wrap_elements(full = pB_row) / wrap_elements(full = legend_plot)) +
  plot_layout(heights = c(1, 0.1))

pB <- pB + plot_annotation(
  title = "B",
  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0,
                                          margin = margin(b = 2)),
                plot.margin = margin(2,2,2,2))
)
pB

#### OLD 
pB1 = plot_data(x = "Tgut368A", "")+
  theme(legend.position = "none")
pB2 = plot_data(x = "CR1_LINE", "")+
  theme(legend.position = "none")
pB3 = plot_data(x = "CR1_A", "")+
  theme(legend.position = "none")
pB4 = plot_data(x = "CR1_B", "")+
  theme(legend.position = "none")

# Get Legend
legend_plot <- get_legend(
  plot_data(x = "CR1_A", "")+
    theme(legend.position = "bottom"))
legend_plot <- wrap_elements(legend_plot)

# Merge
pB=((pB1+pB2+pB3+pB4) + plot_layout(widths=c(1,1,1,1)))


# Make Y axis 
y_axis <- ggplot() +
  labs(y = "Molar ellipticity") +
  theme_void() +
  theme(
    axis.title.y = element_text(angle = 90, size = 10, vjust=0),
    plot.margin = margin(r = 0)
  )
# Make X axis 
x_axis <- ggplot() +
  labs(x = "Wavelength (nm)") +
  theme_void() +
  theme(
    axis.title.x = element_text(size = 10, vjust=0),
    plot.margin = margin(t = 0)
  )
# Add y first
pB_with_y <- (y_axis | pB) +
  plot_layout(widths = c(0.01, 1))
pB_with_x <- pB_with_y /
  x_axis +
  plot_layout(heights = c(1, 0.01))
pB_with_x 
# Add legend 

pB <- (wrap_elements(full = pB_with_x) / legend_plot) + 
  plot_layout(heights=c(1,0.05))

pB <- pB + plot_annotation(title = "B", theme=theme(plot.title = element_text(size = 18, face="bold")))
pB


################################################################################
################################################################################
################################################################################
# CODE FOR FIGURE 5C, CHIP SEQ DATA IN CHICKEN


# Input files
mat_file <- "experimental/matrix/SRR9603969_vs_SRR9603970.TSS.matrix.gz"
groups_file <- "helpfiles/chicken.v23.groups.txt"

# Read the group data 
grouptib <- groups_file %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(chr=V1, length=V2, group=V3) %>% mutate(group = str_to_title(str_to_lower(group)))

## Read metadata header (first line, starts with '@')
header_line <- read_lines(mat_file, n_max = 1)
meta <- fromJSON(sub("^@", "", header_line))

upstream   <- meta$upstream[[1]]
downstream <- meta$downstream[[1]]
binsize    <- meta$`bin size`[[1]]
sample_lab <- meta$sample_labels[[1]]

n_bins <- (upstream + downstream) / binsize

## Read the matrix body as a tibble
# 6 metadata columns (chr, start, end, name, score, strand) + N signal bins
meta_names   <- c("chr", "start", "end", "name", "score", "strand")
signal_names <- paste0("bin_", seq_len(n_bins))

mat <- read_tsv(
  mat_file,
  skip       = 1,
  col_names  = c(meta_names, signal_names),
  na         = c("nan", "NA", "NaN"),
  col_types  = cols(
    chr    = col_character(),
    start  = col_double(),
    end    = col_double(),
    name   = col_character(),
    score  = col_character(),
    strand = col_character(),
    .default = col_double()
  )
) %>%
  mutate(gene_id = row_number())

stopifnot(ncol(mat) - length(meta_names) - 1 == n_bins)  # -1 for gene_id


## Bin coordinates relative to TSS 
bin_starts <- seq(-upstream, downstream - binsize, by = binsize)
bin_mid    <- bin_starts + binsize / 2
bin_lookup <- tibble(bin = signal_names, position = bin_mid)

## Join chromosome groups to genes
group_levels <- c("Macro", "Micro", "Dot")
mat <- mat %>%
  left_join(grouptib %>% select(chr, group), by = "chr") %>%
  mutate(group = factor(group, levels = group_levels))

missing_chr <- mat %>% filter(is.na(group)) %>% distinct(chr) %>% pull(chr)
if (length(missing_chr) > 0) {
  message("Note: no group assigned for: ", paste(missing_chr, collapse = ", "),
          " (excluded from group-based plots)")
}

## Long format: one row per gene x bin 
mat_long <- mat %>%
  filter(!is.na(group)) %>%
  select(gene_id, chr, group, all_of(signal_names)) %>%
  pivot_longer(all_of(signal_names), names_to = "bin", values_to = "value") %>%
  left_join(bin_lookup, by = "bin")

## Per-group mean profile (for combined profile plot) --------
profile_df <- mat_long %>%
  group_by(group, position) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    se       = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups  = "drop"
  ) %>%
  left_join(
    mat_long %>% distinct(gene_id, group) %>% count(group, name = "n_genes"),
    by = "group"
  ) %>%
  mutate(group_lab = fct_reorder(paste0(group, " (n=", n_genes, ")"), as.integer(group)))

## Combined profile plot, without legend
pC_base <- ggplot(profile_df, aes(x = position, y = mean_val, color = group_lab, fill = group_lab)) +
  geom_ribbon(aes(ymin = mean_val - se, ymax = mean_val + se), alpha = 0.3, color = NA) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values=type_col)+
  scale_color_manual(values=type_col)+
  labs(
    x = "Distance from TSS (bp)",
    y = "Mean ChIP / Control ratio (RPKM)",
    title = "C",
  ) +
  theme_minimal(base_size = 12)+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    plot.title=element_text(size=18, face="bold"), 
    strip.background = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.ticks.length = unit(4, "pt"),
    plot.margin = margin(2, 2, 2, 2)
  )

pC_noleg <- pC_base + theme(legend.position = "none")

pC_legend <- get_legend(
  pC_base +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.35, "cm"),
      legend.spacing.y = unit(1, "pt")
    ) +
    guides(fill = guide_legend(nrow = 3), color = guide_legend(nrow = 1))
)
pC_legend <- wrap_elements(full = pC_legend)

# Combine
pC <- (pC_noleg / pC_legend) + plot_layout(heights = c(1, 0.13))
pC



################################################################################
################################################################################
################################################################################
# COMBINE THE THREE PARTS INTO ONE
top_panels <- wrap_elements(pA) + wrap_elements(pC) + plot_layout(widths=c(3.5,1))
top_panels


final=wrap_elements(top_panels) / plot_spacer() /  wrap_elements(pB)  + plot_layout(heights=c(1.5,-0.3,1))
final
outfile=paste(plotdir,"Fig5_G4.svg", sep="")
ggsave(outfile,plot = final,scale = 1,dpi = 600,limitsize = TRUE,width=10,height=6)
outfile=paste(plotdir,"Fig5_G4.png", sep="")
ggsave(outfile,plot = final,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=6)

# The combined plot has a lot ow white space, try two separate:
outfile=paste(plotdir,"Fig5_TOP.svg", sep="")
ggsave(outfile,plot = top_panels,scale = 1,dpi = 600,limitsize = TRUE,width=10,height=4)
outfile=paste(plotdir,"Fig5_BOTTOM.svg", sep="")
ggsave(outfile,plot = pB,scale = 1,dpi = 600,limitsize = TRUE,width=10,height=4)
outfile=paste(plotdir,"Fig5_TOP.png", sep="")
ggsave(outfile,plot = top_panels,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=4)
outfile=paste(plotdir,"Fig5_BOTTOM.png", sep="")
ggsave(outfile,plot = pB,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=4)
# Combine and add sequences in Adobe Illustrator
