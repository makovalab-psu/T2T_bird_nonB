################################################################################
### R code for plotting ChIP-seq data for chicken
### written by Linnéa Smeds Sept 11, 2026.

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(viridis)
require(patchwork)
library(jsonlite)
library(ggplot2)

plotdir="plots/"

# Set colors 
vcolors=viridis(7)
type_col=c("#440154","#9460a1", "#e3cfe8")
g4col=vcolors[6]
################################################################################
# CODE FOR FIGURE S11 - ChIP-SEQ AROUND TSS
################################################################################

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

## Combined profile plot, one line per group ------------------
pd <- ggplot(profile_df, aes(x = position, y = mean_val, color = group_lab, fill = group_lab)) +
  geom_ribbon(aes(ymin = mean_val - se, ymax = mean_val + se), alpha = 0.3, color = NA) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values=type_col)+
  scale_color_manual(values=type_col)+
  labs(
    x = "Distance from TSS (bp)",
    y = "Mean ChIP / Control ratio (RPKM)",
    title = "D",
  ) +
  theme_minimal(base_size = 12)+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    plot.title=element_text(size=18, face="bold"), 
    strip.background = element_blank(),
    legend.key = element_rect(fill = NA, color = NA),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.length = unit(4, "pt"),
  )+
  guides(fill = guide_legend(nrow = 3))
pd

ggsave("plots/TSS_profile_by_group_ggplot.png", p_profile,
       width = 7, height = 4.5, dpi = 300)

## Separate heatmap for each group 
# Shared color scale (99th percentile across ALL data) so heatmaps are comparable
cap <- quantile(mat_long$value, 0.99, na.rm = TRUE)

# Order genes within each group by mean signal (highest first), like plotHeatmap
gene_order_df <- mat_long %>%
  group_by(group, gene_id) %>%
  summarise(mean_sig = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(group) %>%
  arrange(desc(mean_sig), .by_group = TRUE) %>%
  mutate(row_idx = row_number()) %>%
  ungroup() %>%
  select(group, gene_id, row_idx)

heat_df <- mat_long %>%
  left_join(gene_order_df, by = c("group", "gene_id")) %>%
  mutate(value_capped = pmin(value, cap))

# Heatmap plot function 
make_heatmap <- function(data, g, title, showleg=TRUE) {
  
    df_g <- data %>% filter(group == g)
    if (nrow(df_g) == 0) return(NULL)
    
    n_genes_g <- n_distinct(df_g$gene_id)
    
    p <- ggplot(df_g, aes(x = position, y = row_idx, fill = value_capped)) +
      geom_tile(show.legend = showleg) +
      scale_fill_gradientn(colors = c("white", "orange", "red"),
                           name = "ratio", limits = c(0, cap)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
      labs(x = "Distance from TSS (bp)",
           y = paste0("Genes (n=", n_genes_g, ")"), 
           title = title) +
      scale_x_continuous(lim=c(-3000,3000), expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme(
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(size=18, face="bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
      )
    p
}

pa<-make_heatmap(heat_df, "Macro", "A", FALSE)
pa
pb<-make_heatmap(heat_df, "Micro", "B", FALSE)
pb
pc<-make_heatmap(heat_df, "Dot", "C")
pc

p_final11<- wrap_elements(pa) + wrap_elements(pb) + wrap_elements(pc) + plot_layout(widths=c(1,1,1.3))
p_final11

outfile=paste(plotdir,"Fig21_chipseq_heatmap.png", sep="")
ggsave(outfile,plot = p_final11,scale = 1,dpi = 300,limitsize = TRUE,width=7,height=9)


################################################################################
# CODE FOR FIGURE S12 - ChIP-SEQ PEAKS IN OTHER FUNCTIONAL REGIONS 
################################################################################


enfile<-"experimental/peaks/functional.enrichment.tsv"

ENR<-enfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(PrintName=case_when(Class=="promoter" ~ "1 kb upstream\nfrom TSS\n(putative promoters)",
    Class=="introns" ~ "Intronic",
    Class=="intergenic" ~ "Intergenic",
    Class=="UTR5" ~ "5'UTRs",
    Class=="UTR3" ~ "3'UTRs",
    TRUE ~ Class)) %>%
  mutate(Group=case_when(Group=="macro" ~ "Macro",
                         Group=="micro" ~ "Micro",
                         Group=="dot" ~ "Dot",
                         TRUE ~ Group)) %>% filter(Group!="genome") %>% 
  pivot_longer(c(Enrichment_gw, Enrichment_grp), names_to="Baseline", values_to="Enrichment")

ENR$PrintName <- factor(ENR$PrintName, levels=c("1 kb upstream\nfrom TSS\n(putative promoters)","5'UTRs","CDS", "3'UTRs", "Intronic", "Intergenic")) #"lncRNA", 
ENR$Group <- factor(ENR$Group, levels=c("Macro","Micro","Dot"))



# Function to make consistent bar plots
make_enrich_plot <- function(data, title, yaxis) {
  
  p<-ggplot(data, aes(x=PrintName, y=Enrichment)) + 
    geom_bar(stat="identity", alpha=0.6, fill=g4col, color=g4col) +
  facet_wrap(~Group, nrow=4) +
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8))) +
    labs(x = NULL, y = yaxis, title=title) +  
    theme_minimal(base_size = 13) +
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      legend.position = "none", # We'll keep legend only on one
      axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title=element_text(size=18, face="bold"),
      axis.ticks.y=element_line())
    return(p)
}

pa<-make_enrich_plot(ENR%>%filter(Baseline=="Enrichment_gw"), "A", "Fold enrichment compared to genome-wide")
pa
pb<-make_enrich_plot(ENR%>%filter(Baseline=="Enrichment_grp"), "B", "Fold enrichment comapred to chromosome category")
pb

p_final12<- wrap_elements(pa) + wrap_elements(pb)
p_final12

outfile=paste(plotdir,"FigS22_chipseq_enrichment.png", sep="")
ggsave(outfile,plot = p_final12,scale = 1,dpi = 300,limitsize = TRUE,width=8,height=10)




