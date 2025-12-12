################################################################################
### R code for plotting non-B motif enrichment in satellite and non-satellite 
### part of introns in the zebra finch T2T genome. 
### written by Linnéa Smeds Sept-10-2025
################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(viridis)
require(patchwork)
plotdir="plots/"

# Input files 
enfile1="repeats/group_intron_enrichment.tsv"
enfile2="repeats/compartment_intron_enrichment.tsv"
lenfile="repeats/intron_TR_lengthsummary.tsv"
repfile="repeats/intron_TR_repeatsummary.tsv"

# Reading in the data 
en_tib1<-enfile1 %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(Set="ChrType")
en_tib2<-enfile2 %>% read.table(header=TRUE) %>% as_tibble() %>%
  rename(Group=Compartment)  %>% mutate(Set="Compartment")
lentib<-lenfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Set=case_when(Group=="A" ~ "Compartment",
                       Group=="B" ~ "Compartment",
                       TRUE ~ "ChrType"))
reptib<-repfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Set=case_when(Group=="A" ~ "Compartment",
                       Group=="B" ~ "Compartment",
                       TRUE ~ "ChrType"))

# Category and Compartment data together 
DATA <- en_tib1 %>% bind_rows(en_tib2) %>% 
  mutate(PrintName=case_when(Subset=="TRF" ~ "TR",  
                             Subset=="noTRF" ~ "no TR",
                            TRUE ~ Subset)) %>%
  mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                  Group=="micro" ~ "Micro",
                                  Group=="microdot" ~ "Dot",
                                  TRUE ~ Group)) 
DATA$NonB <- factor(DATA$NonB, levels=c("ALL","APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))
DATA$PrintName <- factor(DATA$PrintName, levels=c("TR","no TR"))

# Length with fractions of Bp
LEN <- lentib %>% group_by(Group, Set) %>% 
  mutate(totbp=sum(Length), frac = Length/totbp) %>%
  ungroup() %>% mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                                Group=="micro" ~ "Micro",
                                                Group=="microdot" ~ "Dot",
                                                TRUE ~ Group)) %>%
  mutate(ClassCol = ifelse(Subset == "noTRF", "No", as.character(Classification)))
LEN$Classification<-factor(LEN$Classification, levels=c("Macro", "Micro", "Dot", "A", "B"))
LEN$ClassCol<-factor(LEN$ClassCol, levels=c("No", "Macro", "Micro", "Dot", "A", "B"))


# Repeats 
REP <-reptib %>% mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                                Group=="micro" ~ "Micro",
                                                Group=="microdot" ~ "Dot",
                                                TRUE ~ Group)) 
REP$Repeat<-factor(REP$Repeat, levels=c("1-4mer", "5-10mer", "11-50mer", "51-100mer", "100+mer"))
REP$Classification<-factor(REP$Classification, levels=c("Macro", "Micro", "Dot", "A", "B"))


# Set colors 
vcolors=viridis(8)
classification_fill <- c("Macro" = "#896192","Micro" = "#B99AC1", "Dot" = "#E8DCEB", "No" = "white", "A"="#E38AAA", B="#6F9DD0")
classification_colors <- c("Macro" = "#896192","Micro" = "#B99AC1", "Dot" = "#E8DCEB","No" = "lightgray", "A"="#E38AAA", B="#6F9DD0")
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
    labs(x = NULL, y = "Fold enrichment in introns") +  
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      legend.position = "none", # We'll keep legend only on one
      axis.text.x = element_text(size=10, angle=45, vjust=1, hjust=1),
      panel.grid = element_blank(),
      axis.ticks.y=element_line(),
      legend.title=element_blank(),
      plot.title=element_text(size=18, face="bold"),
    )+
    ggtitle("A")
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
    scale_color_manual(values = classification_colors)+
  ggtitle("B")
   return(p)
} 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO PLOT REPEAT CLASS NUMBER 
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
  ggtitle("C")
  if (!is.null(y_max)) {
    p <- p + 
      coord_cartesian(ylim = c(0, y_max)) +
      scale_y_continuous(expand = c(0.0, 0.0))
    
  }
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
  
# Repeats, for class plot 
data3_macro <- REP %>% filter(Classification == "Macro")
data3_micro <- REP %>% filter(Classification == "Micro")
data3_dot <- REP %>% filter(Classification == "Dot")

################################################################################
#  MAKE NINE SEPARATE PLOTS 
################################################################################
# Create the individual plots
p1 <- make_plot(data1_macro,25) + theme(axis.text.x = element_blank()) +
  theme(legend.position = "top",  legend.justification = "right", axis.title.y= element_blank()) + 
  guides(fill = guide_legend(nrow = 1))
p2 <- make_plot(data1_micro,25) + theme(axis.text.x = element_blank(), plot.title=element_blank())
p3 <- make_plot(data1_dot,25) + theme(axis.title.y= element_blank(), plot.title=element_blank())
p4<-make_pie(data2_macro)
p5<-make_pie(data2_micro)+ theme(plot.title=element_blank())
p6<-make_pie(data2_dot)+ theme(plot.title=element_blank())
p7<-make_bar(data3_macro,160000)+ theme(axis.text.x = element_blank(), axis.title.y= element_blank()) 
p8<-make_bar(data3_micro,160000)+ theme(axis.text.x = element_blank(), plot.title=element_blank()) 
p9<-make_bar(data3_dot,160000)+ theme(axis.title.y = element_blank(), plot.title=element_blank()) 


################################################################################
# COMBINE THE PLOTS INTO ONE 
################################################################################

# Combine plots with patchwork with different heights
pcomb<-((p1 / p2 / p3) | (p4 / plot_spacer() / p5 / plot_spacer() / p6 +plot_layout(heights=c(3,1,3,1,3))) | (p7 / p8 / p9)) + plot_layout(widths=c(3,1,1)) 
pcomb
outfile=paste(plotdir,"FigSX_Introns_and_TRs.png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 600,limitsize = TRUE,width=10,height=9)


################################################################################
# REPEAT ABOVE: EXTRACT SUBSETS OF DATA, COMPARTMENTS
# EXTRACT CHROMOSOME CATEGORIES 
# Enrichment, first panels 
data1_a <- DATA %>% filter(Classification == "A") %>%
  filter(NonB!="ALL")
data1_b <- DATA %>% filter(Classification == "B") %>%
  filter(NonB!="ALL")
# Length, for pie
data2_a <- LEN %>% filter(Classification == "A")
data2_b <- LEN %>% filter(Classification == "B")

# Repeats, for class plot 
data3_a <- REP %>% filter(Classification == "A")
data3_b <- REP %>% filter(Classification == "B")


################################################################################
#  MAKE SIX SEPARATE PLOTS 
# Create the individual plots
p10 <- make_plot(data1_a,22) + theme(axis.text.x = element_blank()) +
  theme(legend.position = "top",  legend.justification = "right") + 
  guides(fill = guide_legend(nrow = 1))
p11 <- make_plot(data1_b,22) + theme(plot.title=element_blank())
p12<-make_pie(data2_a)
p13<-make_pie(data2_b)+ theme(plot.title=element_blank())
p14<-make_bar(data3_a,130000)+ theme(axis.text.x = element_blank()) 
p15<-make_bar(data3_b,130000)+ theme(plot.title=element_blank()) 


################################################################################
# COMBINE THE PLOTS INTO ONE 
# Combine plots with patchwork with different heights
pcomb2<-((p10 / p11) | (p12 / plot_spacer() / p13 +plot_layout(heights=c(3,1,3))) | (p14 / p15)) + plot_layout(widths=c(3,1,1)) 
pcomb2
outfile2=paste(plotdir,"FigS15_Introns_and_TRs_ABCompartments.png", sep="")
ggsave(outfile2,plot = pcomb2,scale = 1,dpi = 600,limitsize = TRUE,width=10,height=7)

