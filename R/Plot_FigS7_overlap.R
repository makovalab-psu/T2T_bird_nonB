################################################################################
### R code for plotting non-B motif overlaps in birds
### written by Linnéa Smeds 11-July-2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(ggupset)
require(ggtext)
require(patchwork)
require(viridis)
library(dplyr)
library(purrr)
options(scipen = 999, decimals=2)
plotdir="plots/"

# SPECIES LIST 
spfile="helpfiles/species_list.txt"
sptib<-spfile %>% read.table(header=FALSE) %>% as_tibble()

# Colorscheme
viridis_colors <- viridis(4)

################################################################################
#  FUNCTION TO MAKE PAIRWISE OVERLAP PLOTS
################################################################################

make_pairwise<-function(data, title) {
  
  ppw<-ggplot(data, aes(x=NB2, y=forcats::fct_rev(NB1), fill=Frac)) +
    geom_tile(color="white", show.legend = FALSE) +
    geom_text(aes(label=round(Frac*100, digits=1), color=textcol), size=4, show.legend=FALSE) +
    facet_wrap(~PrintName, nrow=6, scales='free', strip.position="right") +
    labs(x="", y="") +
    scale_fill_gradient(low = "#f2f6ff", high = viridis_colors[2], name="Overlap", na.value = 'white') +
    scale_color_manual(values=c("black","white")) +
    theme(panel.grid.major = element_line(colour = 'white'),
          panel.grid.minor = element_line(colour = 'white'),
          panel.background = element_rect(fill = 'white', colour = 'white'),
          panel.spacing = unit(1, "lines"),
          legend.position="bottom",
          plot.title = element_text(size=18, face="bold"),
          legend.title=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.text=element_text(size=11),
          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5),
          strip.background = element_rect(fill='white'),
          strip.text = element_text(size=14, face="bold"))+
    labs(title = title)
  return(ppw)
}
  
################################################################################
#  FUNCTION TO MAKE MULTI OVERLAP PLOTS
################################################################################
make_upset<-function(data, title) {
  
  # 1) Define color scale based on number of overlap types
  overlap_levels <- sort(unique(data$n_types))
  upset_colors <- viridis(length(overlap_levels))
  
  # 2) Compute how many combinations appear at each overlap level
  counts <- data |>
    dplyr::select(Comb, n_types) |>
    distinct() |>
    count(n_types, name = "n")
  
  # Expand colors to match combmatrix structure in ggupset
  repeats <- counts$n * seq_along(counts$n)
  color_vector <- rep(upset_colors, times = repeats)
  
  # Check that there are no lines without data 
  data <- data |>
    dplyr::mutate(Bp = as.numeric(Bp))
  
  pup<-ggplot(data, aes(x=Types)) +
    geom_col(aes(y=Bp/1000000, fill=vcolor), alpha=0.8, show.legend = FALSE) +
    facet_wrap(~PrintName, nrow=6, scales="free_y", strip.position="right")+
    scale_fill_identity() +
    scale_x_upset(order_by="degree")+
    theme_combmatrix(
      combmatrix.panel.point.color.fill = color_vector,
      combmatrix.panel.line.color = color_vector,
      combmatrix.panel.point.color.empty = "gray90",
      combmatrix.panel.striped_background.color.two = "grey95",
      combmatrix.panel.point.size = 2,
    ) +
    scale_y_continuous(name='Mb')+
    theme(panel.background = element_rect(fill = 'white', colour="gray"),
          panel.grid.major = element_line(colour = 'gray95'),
          strip.background = element_rect(fill = 'white', colour="white"),
          plot.title = element_text(size=18, face="bold"),
          strip.text = element_markdown(size=14, face="bold"),
          panel.spacing = unit(2, "lines"),
          axis.text.y=element_text(size=12),
          axis.title.y = element_text(size=14,vjust = 0.5), # Adjust 'r' to a smaller value
          axis.ticks.x=element_blank())+
    labs(title=title, x="")
  
  return(pup)
}

################################################################################
#  LOOP OVER SPECIES, READ IN THE DATA AND PLOT 
################################################################################

upplots <- vector("list", length(sptib$V3))
names(upplots) <- sptib$V1
pwplots <- vector("list", length(sptib$V3))
names(pwplots) <- sptib$V1
header<-c("A", "B", "C", "D", "E", "F", "G", "H")

for (i in 1:length(sptib$V1)){
  show(paste("Running pipeline for:",sptib$V1[i]))
  
  # READ IN DATA 
  fileMulti=paste("overlap/",sptib$V3[i],".merged.7types.summary.txt", sep="")
  filePair=paste("overlap/",sptib$V3[i],".merged.7types.pairwise.txt", sep="")
  
  # MAKE TIBBLES
  tibM<-fileMulti %>% read.table(header=TRUE) %>% as_tibble() %>% 
    rename(Comb=NonB, Bp=Overlap) %>% mutate(Types=strsplit(Comb, "-")) %>% 
    filter(Region!="autosomes") %>%
    mutate(PrintName=case_when(Region=="chrW" ~ "Chr W",
                               Region=="chrZ" ~ "Chr Z",
                               Region=="macro" ~ "Macro",
                               Region=="micro" ~ "Micro",
                               Region=="dot" ~ "Dot")) 
  tibM$PrintName=factor(tibM$PrintName, levels=c("Macro", "Micro", "Dot", "Chr Z", "Chr W"))
  
  tibP<-filePair %>% read.table(header=TRUE) %>% as_tibble() %>% 
    mutate(textcol=if_else(Frac>0.5,"white","black")) %>% 
    filter(Region!="autosomes") %>%
    mutate(PrintName=case_when(Region=="chrW_mat" ~ "Chr W",
                               Region=="chrZ_pat" ~ "Chr Z",
                               Region=="chrW" ~ "Chr W",
                               Region=="chrZ" ~ "Chr Z",
                               Region=="macro" ~ "Macro",
                               Region=="micro" ~ "Micro",
                               Region=="dot" ~ "Dot"))
  tibP$NB1=factor(tibP$NB1, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))
  tibP$NB2=factor(tibP$NB2, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))
  tibP$PrintName=factor(tibP$PrintName, levels=c("Macro", "Micro", "Dot","Chr Z", "Chr W"))
  
  # GENERATE PAIRWISE PLOTS 
  pp<-make_pairwise(tibP, title=header[i])
  
  # SUBSET AND GENERATE UPSET PLOTS
  SUBM <- tibM %>% filter(Bp>10000) %>% 
    arrange(desc(Bp)) %>% 
    mutate(vcolor=viridis_colors[lengths(Types)]) %>%
    mutate(n_types = map_int(Types, length))
  pu<-make_upset(SUBM, title=" ")

  upplots[[i]] <- pu
  pwplots[[i]] <- pp
}
###############################################################################
# PLOT COMBINED FOR EACH SPECIES 
for (nm in names(upplots)) {
  pcomb <- pwplots[[nm]] + upplots[[nm]] + plot_layout(widths = c(1, 2.5)) 
  ggsave(
    filename = paste0("plots/", nm, "_Upset.png"),
    plot = pcomb,
    scale = 1,dpi = 300,limitsize = TRUE,width=10,height=13
  )
}
