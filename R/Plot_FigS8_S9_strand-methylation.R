################################################################################
### R code for plotting methylation levels in G4s (strand dependent)
### written by Linnéa Smeds August 27, 2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(viridis)
library(ggh4x)
library(patchwork)
library(purrr)
plotdir="plots/"

setwd("/Users/lbs5874/OneDrive - The Pennsylvania State University/Documents/Projects/ZebraFinch/")

# Input files 
file="methylation/bTaeGut7v0.4_MT_rDNA.allCpG.group.txt"
intergenfile="methylation/bTaeGut7v0.4_MT_rDNA.allCpG.intergenic.group.txt"
sumfile<-"methylation/bTaeGut7v0.4_MT_rDNA.summaryCpG.group.txt"

# Reading in the data 
genetib<-file %>% read.table(header=TRUE) %>% as_tibble() 
intergentib<-intergenfile %>% read.table(header=TRUE) %>% as_tibble() 
sumtib<-sumfile %>% read.table(header=TRUE) %>% as_tibble() 


# Summarize gene regions with median and mean per transcript 
DATA1 <-genetib %>% group_by(Group,Class,Type,Trx) %>% 
  summarize(median=median(Score), mean=mean(Score)) %>% 
  mutate(PrintName=case_when(Class=="lncrna" ~ "lncRNA",  
                             Class=="promoter" ~ "1 kb upstream\nfrom TSS",
                             Class=="all_repeats" ~ "Repeats",
                             Class=="introns" ~ "Intronic",
                             Class=="intergenic" ~ "Intergenic",
                             Class=="UTR5" ~ "5'UTRs",
                             Class=="UTR3" ~ "3'UTRs",
                             TRUE ~ Class)) %>%
  mutate(Classification=case_when(Group=="macro" ~ "Macro",
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
DATA1$Classification <- factor(DATA1$Classification, levels=c("Macro","Micro","Dot"))
DATA1$PrintType <- factor(DATA1$PrintType, levels=c("Background","Coding","Template"))
DATA1$PrintName <- factor(DATA1$PrintName, levels=c("1 kb upstream\nfrom TSS","5'UTRs","CDS", "3'UTRs", "Intronic", "lncRNA", "Intergenic"))
DATA1$Region <- factor(DATA1$Region, levels=c("inside G4","outside G4"))

# For intergenic we use all sites 
DATA2<-intergentib %>%
  mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                  Group=="micro" ~ "Micro",
                                  Group=="dot" ~ "Dot",
                                  TRUE ~ Group)) %>% 
  mutate(Region=case_when(Type=="ignorant" ~ "inside G4",
                          Type=="background" ~ "outside G4",
                          TRUE ~ NA)) %>% mutate(PrintName="Intergenic", median=Score)
DATA2$Classification <- factor(DATA2$Classification, levels=c("Macro","Micro","Dot"))


# Summary per group, for Piechart 

DATA3 <- sumtib %>%
  mutate(PrintName=case_when(Class=="lncrna" ~ "lncRNA",  
                            Class=="promoter" ~ "1 kb upstream\nfrom TSS",
                            Class=="all_repeats" ~ "Repeats",
                            Class=="introns" ~ "Intronic",
                            Class=="intergenic" ~ "Intergenic",
                            Class=="UTR5" ~ "5'UTRs",
                            Class=="UTR3" ~ "3'UTRs",
                            TRUE ~ Class)) %>%
  mutate(Classification=case_when(Group=="macro" ~ "Macro",
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
DATA3$PrintName <- factor(DATA3$PrintName, levels=c("1 kb upstream\nfrom TSS","5'UTRs","CDS", "3'UTRs", "Intronic", "lncRNA", "Intergenic"))
DATA3$PrintType <- factor(DATA3$PrintType, levels=c("Background","Coding","Template"))
DATA3$Region <- factor(DATA3$Region, levels=c("inside G4","outside G4"))
DATA3$Classification <- factor(DATA3$Classification, levels=c("Macro","Micro","Dot"))
DATA3$Region <- factor(DATA3$Region, levels=c("inside G4","outside G4"))

# Define colors
vcolors=viridis(7)
################################################################################
# PLOT

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# separate coding and non coding strand 

SUB <- DATA1 %>% filter(Type!="background")
make_densityplot <- function(data, showleg = TRUE, showx = TRUE, showy=TRUE, ytext=FALSE) {
  p <- ggplot(data, aes(x = median/100, fill = PrintType, color = PrintType)) + 
    geom_density(alpha = 0.6, show.legend = showleg) +
    facet_nested(Classification ~ PrintName, scales = "free_y") +  # single-column facet
    scale_fill_manual(values = c(vcolors[6],"skyblue")) +
    scale_color_manual(values = c(vcolors[6],"skyblue3")) +
    #  theme_minimal(base_size = 12) +
    labs(x="Methylation level", y="Density")+
    scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      strip.text.x = element_text(size = 10),
      strip.text.y = if(ytext==TRUE) element_text(size = 12) else element_blank(),
      strip.background = element_blank(),
      legend.key = element_rect(fill = NA, color = NA),
      axis.title.x = if(showx==TRUE) element_text(size = 10) else element_blank(),
      axis.title.y = if(showy==TRUE) element_text(size = 10) else element_blank(),
      axis.text.y = element_text(size = 10),
      axis.ticks.length = unit(4, "pt"),
      #    panel.spacing.y = unit(1, "lines")
    ) 
  return(p)
}
# Plot each class panel and combine 
p1<-make_densityplot(SUB%>%filter(Class=="promoter"), FALSE, FALSE, TRUE)
p2<-make_densityplot(SUB%>%filter(Class=="UTR5"), FALSE, FALSE, FALSE)
p3<-make_densityplot(SUB%>%filter(Class=="CDS"), FALSE, FALSE, FALSE)
p4<-make_densityplot(SUB%>%filter(Class=="UTR3"), TRUE, TRUE, FALSE)
p5<-make_densityplot(SUB%>%filter(Class=="introns"), FALSE, FALSE, FALSE)
p6<-make_densityplot(SUB%>%filter(Class=="lncrna"), FALSE, FALSE, FALSE, TRUE)

pcomb2<-p1+p2+p3+p4+p5+p6 + plot_layout(widths=c(1,1,1,1,1,1))+plot_annotation(title = "")
pcomb2
outfile=paste(plotdir,"FigS8_Methylation_density_StrandSpec.png", sep="")
ggsave(outfile,plot = pcomb2,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=6)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot piechart with number of genes with methylation

SUB2 <- DATA3 %>% mutate(Meth=if_else(Number>0, "Yes", "No")) %>%
  group_by(PrintName, Classification, PrintType, Meth) %>% summarize(counts=n()) %>%
  mutate(total_counts = sum(counts)) %>%
  ungroup() %>%
  mutate(fraction = counts / total_counts) %>% 
  mutate(MethClass = ifelse(Meth == "No", "No", as.character(Classification)))
SUB2$MethClass <- factor(SUB2$MethClass, levels=c("No", "Macro","Micro","Dot"))

# To add significance, I first need to perform a pariwise z-test 
# Function to run prop.test for Macro vs one other group
test_vs_macro <- function(df, other) {
  tab <- df %>%
    filter(Classification %in% c("Macro", other)) %>%
    select(Classification, Meth, counts) %>%
    pivot_wider(names_from = Meth, values_from = counts)
  m <- as.matrix(tab[, c("No", "Yes")])
  rownames(m) <- tab$Classification
  # run prop.test (two-sample z-test for proportions)
  test <- prop.test(m[, "Yes"], rowSums(m))
  tibble(
    group1 = "Macro",
    group2 = other,
    p.value = test$p.value,
    prop_macro = m["Macro", "Yes"] / sum(m["Macro", ]),
    prop_other = m[other, "Yes"] / sum(m[other, ]),
    diff = prop_other - prop_macro
  )
}

# Run comparisons per PrintName × PrintType
results_pairwise <- SUB2 %>%
  group_by(PrintName, PrintType) %>%
  group_modify(~ bind_rows(
    test_vs_macro(.x, "Dot"),
    test_vs_macro(.x, "Micro")
  )) %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p.value, method = "fdr"),
    signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ "n.s."
    )
  )
results_pairwise


#Make a new tibble with significance 
SIGN<-results_pairwise %>% select(PrintName,PrintType,group2,signif) %>% 
  rename(Classification=group2) 

COMB<-SUB2 %>% full_join(SIGN)

COMB$MethClass <- factor(COMB$MethClass, levels=c("No", "Macro","Micro","Dot"))
COMB$PrintType <- factor(COMB$PrintType, levels=c("Background","Coding","Template"))
COMB$Classification<-factor(COMB$Classification, levels=c("Macro","Micro","Dot"))

# Set colors  
my_colors <- c(
  "No" = "white",
  "Macro" = "#440154",
  "Micro" = "#9460a1", 
  "Dot" = "#e3cfe8" 
)

classification_colors <- c(
  "Macro" = "#896192",
  "Micro" = "#B99AC1", 
  "Dot" = "#E8DCEB" 
)
#text_colors<-c(rep("white",24), rep("black",)

strip_bg_list <- c(
  lapply(classification_colors, function(col) element_rect(fill = col, colour = NA)),
  replicate(9, element_rect(fill = "white", colour = NA), simplify = FALSE)
)

p3<-ggplot(COMB, aes(x = 1, y = fraction, fill = MethClass, color=MethClass)) +
  geom_bar(stat = "identity", width = 1, alpha=0.6, show.legend = FALSE) +  # single bar per pie
  coord_polar(theta = "y", clip = "off") +  # pie chart
  geom_text(data=SUB2%>%filter(Meth=="Yes"), aes(x = 1, y = 0, label = total_counts),
            color = "black", size = 3, fontface = "italic", vjust=-3)+
  geom_text(data=COMB%>%filter(Classification=="Micro" | Classification=="Dot"), aes(x = 1, y = 0.5, label = signif),
            color = "black", size = 4, fontface = "plain", vjust=0)+
  facet_nested(
    Classification + PrintType ~ PrintName # nested facets: rows = Classification + PrintType, cols = PrintName
#    strip = strip_nested(
#      background_y = strip_bg_list
#    )
  )+
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text.y=element_text(size=12),
    strip.placement.y = "left",
    strip.text.y.left = element_text(angle = 0),
  ) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  labs(fill = "Has CpG site", color = "Has CpG site")
p3

outfile=paste(plotdir,"FigS9_Piechart_Regions_with_or_without_CpG_Signif.png", sep="")
ggsave(outfile,plot = p3,scale = 1,dpi = 300,limitsize = TRUE,width=7,height=10)

