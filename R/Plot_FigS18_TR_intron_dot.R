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
enfile1="repeats/bTaeGut7v0.4_MT_rDNA.intron_enrichment.group.tsv"
enfile2="repeats/bTaeGut7v0.4_MT_rDNA.intron_enrichment.compartment.tsv"
lenfile="repeats/bTaeGut7v0.4_MT_rDNA.introns_TR.lengthsummary.tsv"
allrepfile="repeats/bTaeGut7v0.4_MT_rDNA.introns_TR.dot.compart.lengths.tsv"

# Reading in the data 
en_tib1<-enfile1 %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(Set="ChrType")
en_tib2<-enfile2 %>% read.table(header=TRUE) %>% as_tibble() %>%
  rename(Group=Compartment, Coverage=Density)  %>% mutate(Set="Compartment")
lentib<-lenfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Set=case_when(Group=="A" ~ "Compartment",
                       Group=="B" ~ "Compartment",
                       TRUE ~ "ChrType"))
reptib<-allrepfile %>% read.table(header=TRUE) %>% as_tibble() 

# Category and Compartment data together 
DATA <- en_tib1 %>% bind_rows(en_tib2) %>% 
  mutate(PrintName=case_when(Subset=="TRF" ~ "tandem repeat",  
                             Subset=="noTRF" ~ "outside tandem repeat",
                            TRUE ~ Subset)) %>%
  mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                  Group=="micro" ~ "Micro",
                                  Group=="dot" ~ "Dot",
                                  TRUE ~ Group)) 
DATA$NonB <- factor(DATA$NonB, levels=c("Any","APR", "DR", "STR", "IR", "TRI", "G4", "Z"))
DATA$PrintName <- factor(DATA$PrintName, levels=c("tandem repeat","outside tandem repeat"))

# Length with fractions of Bp
LEN <- lentib %>% group_by(Group, Set) %>% 
  mutate(totbp=sum(Length), frac = Length/totbp) %>%
  ungroup() %>% mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                                Group=="micro" ~ "Micro",
                                                Group=="dot" ~ "Dot",
                                                TRUE ~ Group)) %>%
  mutate(ClassCol = ifelse(Subset == "noTRF", "No", as.character(Classification)))
LEN$Classification<-factor(LEN$Classification, levels=c("Macro", "Micro", "Dot", "A", "B"))
LEN$ClassCol<-factor(LEN$ClassCol, levels=c("No", "Macro", "Micro", "Dot", "A", "B"))


# Repeats 
REP <-reptib
REP$Comp<-factor(REP$Comp, levels=c("A", "B"))
REP$RepUnitLen<-as.numeric(REP$RepUnitLen)

# Set colors 
vcolors=viridis(7)
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
      axis.text.x = element_text(size=12),
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
    scale_color_manual(values = classification_colors)
   return(p)
} 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO PLOT REPEAT CLASS NUMBER 
make_bar <- function(data, y_max = NULL) {
  p<-ggplot(data, aes(x = RepUnitLen, fill = Comp, color = Comp)) + 
    geom_bar(position=position_dodge(), alpha=0.6, show.legend = FALSE) +
    scale_fill_manual(values=classification_colors) +
    scale_color_manual(values=classification_colors) +
    scale_x_continuous(expand = c(0.01, 0.01), limits=c(0,100)) +
    facet_wrap(Comp~., nrow=2, strip.position = "right")+
    labs(x ="Tandem repeat unit length (bp)", y="Annotated TR in introns") +  
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      strip.background = element_blank(),
      axis.text.x = element_text(size=10),
      axis.title.y=element_text(size=12),
      panel.grid = element_blank(),
      strip.text = element_text(size=14),
      axis.ticks.x=element_blank(),
      plot.title=element_text(size=18, face="bold"),
    )+
  ggtitle("B")
  if (!is.null(y_max)) {
    p <- p + 
      coord_cartesian(ylim = c(0, y_max)) +
      scale_y_continuous(expand = c(0.0, 0.0))
    
  }
  return(p)
}

################################################################################
# EXTRACT SUBSETS OF DATA FOR COMPARTMENTS
################################################################################
# Enrichment, first panels 
data1_a <- DATA %>% filter(Classification == "A") %>%
  filter(NonB!="Any")
data1_b <- DATA %>% filter(Classification == "B") %>%
  filter(NonB!="Any")
# Length, for pie
data2_a <- LEN %>% filter(Classification == "A")
data2_b <- LEN %>% filter(Classification == "B")
# Repeats, for class plot 
data3_a <- REP %>% filter(Comp == "A")
data3_b <- REP %>% filter(Comp == "B")


################################################################################
#  MAKE SIX SEPARATE PLOTS 
# Create the individual plots
p1 <- make_plot(data1_a,22) + theme(axis.text.x = element_blank()) +
  theme(legend.position = "top",  legend.justification = "right") + 
  guides(fill = guide_legend(nrow = 1))
p2 <- make_plot(data1_b,22) + theme(plot.title=element_blank())
p3<-make_pie(data2_a)
p4<-make_pie(data2_b)
p5<-make_bar(data3_a,20000)+ theme(axis.text.x = element_blank()) 
p6<-make_bar(data3_b,20000)+ theme(plot.title=element_blank()) 


################################################################################
# COMBINE THE PLOTS INTO ONE 

# Place pie charts inside enrichment plots
pA1<- p1 + inset_element(p3, left = 0.6, bottom = 0.3, right = 0.99, top = 0.99)
pA2<- p2 + inset_element(p4, left = 0.6, bottom = 0.3, right = 0.99, top = 0.99)
pB3<- p3 + inset_element(p6, left = 0.6, bottom = 0.3, right = 0.99, top = 0.99)
pA <- (pA1 / pA2)
pA
# Combine with repeat histograms 
pB <- (p5 / p6)

pcomb <- (pA | pB) + plot_layout(widths = c(2, 1))
pcomb

outfile=paste(plotdir,"FigS18_Introns_and_TRs_ABCompartments.png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 300,limitsize = TRUE,width=10,height=10)

