################################################################################
### R code for plotting non-B motif enrichment in repeats in Zebra finch.
### written by Linnéa Smeds Sept 10, 2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(tidyverse)
require(viridis)
require(patchwork)
plotdir="plots/"
setwd("/Users/lbs5874/Documents/Projects/ZebraFinch/")

################################################################################
################################################################################
################################################################################
# PART A, General repeat enrichment 

# Input files 
repfile="repeats/genome_repeat_enrichment.tsv"
lenfile="repeats/TE_TRF_SAT_lengths.ordered.txt"
repgroupfile="repeats/group_repeat_enrichment.tsv"
lenfile1="repeats/macro.TE_TRF_SAT_lengths.txt"
lenfile2="repeats/micro.TE_TRF_SAT_lengths.txt"
lenfile3="repeats/microdot.TE_TRF_SAT_lengths.txt"

# Reading in the data 
rep_tib1<-repfile %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(Group="GW")
rep_tib2<-repgroupfile %>% read.table(header=TRUE) %>% as_tibble()
len_tib<-lenfile %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Order=V1, Repeat=V2, GW=V3)
len_tib1<-lenfile1 %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Repeat=V1, macro=V2)
len_tib2<-lenfile2 %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Repeat=V1, micro=V2)
len_tib3<-lenfile3 %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Repeat=V1, microdot=V2)

LEN<-len_tib %>% full_join(len_tib1) %>% full_join(len_tib2) %>% 
  full_join(len_tib3) %>% mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>%
  mutate(Lenfilter=if_else(GW>=10000, "Yes", "No")) %>%
  pivot_longer(cols=c(GW,macro, micro, microdot), names_to="Group", values_to="Length")

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
                                  Group=="microdot" ~ "Dot",
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
  
  
DATA$nonB <- factor(DATA$nonB, levels=c("APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z","Any"))
DATA$Class <- factor(DATA$Class, levels=c("TEs", "TRs", "Satellites", "RNA", "Other"))
DATA$Classification<-factor(DATA$Classification, levels=c("Macro", "Micro", "Dot"))

################################################################################
# PLOT REPEAT ENRICHMENT; ONLY Repeat types with more than 10kb in total, and skip 
# "Satellite" as a joint class 
# Also tried plotting without logscale, but that doesn't work at all, the 
# highest enrichment value is way too high. 

SUB1<-DATA %>% filter(Repeat!="Satellite/Satellite" & Group!="GW" & Lenfilter=="Yes")
#outfile=paste(plotdir,"Fig3_NonB_enrichment_repeats_log2_5kb_",group,".png", sep="")
outfile=paste(plotdir,"Fig3_NonB_enrichment_repeats_log2_facet.png", sep="")

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
  ggtitle("A")
  scale_y_discrete(position="left")
pA

################################################################################
################################################################################
################################################################################
# PART B,C,D, Intron TR enrichment 


# Input files 
enfile="repeats/group_intron_enrichment.tsv"
lenfile="repeats/intron_TR_lengthsummary.tsv"
repfile="repeats/intron_TR_repeatsummary.tsv"

# Reading in the data 
en_tib<-enfile %>% read.table(header=TRUE) %>% as_tibble()
lentib<-lenfile %>% read.table(header=TRUE) %>% as_tibble() 
reptib<-repfile %>% read.table(header=TRUE) %>% as_tibble() 

# Category and Compartment data together 
DATA <- en_tib  %>% 
  mutate(PrintName=case_when(Subset=="TRF" ~ "TR",  
                             Subset=="noTRF" ~ "not TR",
                             TRUE ~ Subset)) %>%
  mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                  Group=="micro" ~ "Micro",
                                  Group=="microdot" ~ "Dot",
                                  TRUE ~ Group)) 
DATA$NonB <- factor(DATA$NonB, levels=c("ALL","APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))
DATA$PrintName <- factor(DATA$PrintName, levels=c("TR","not TR"))

# Length with fractions of Bp
LEN <- lentib %>% group_by(Group) %>% 
  mutate(totbp=sum(Length), frac = Length/totbp) %>%
  ungroup() %>% mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                                Group=="micro" ~ "Micro",
                                                Group=="microdot" ~ "Dot",
                                                TRUE ~ Group)) %>%
  mutate(ClassCol = ifelse(Subset == "noTRF", "No", as.character(Classification)))
LEN$Classification<-factor(LEN$Classification, levels=c("Macro", "Micro", "Dot"))
LEN$ClassCol<-factor(LEN$ClassCol, levels=c("No", "Macro", "Micro", "Dot"))


# Repeats 
REP <-reptib %>% mutate(Classification=case_when(Group=="macro" ~ "Macro",
                                                 Group=="micro" ~ "Micro",
                                                 Group=="microdot" ~ "Dot",
                                                 TRUE ~ Group)) 
REP$Repeat<-factor(REP$Repeat, levels=c("1-4mer", "5-10mer", "11-50mer", "51-100mer", "100+mer"))
REP$Classification<-factor(REP$Classification, levels=c("Macro", "Micro", "Dot"))


# Set colors 
vcolors=viridis(8)
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
    scale_color_manual(values = classification_colors)+
    ggtitle("C")
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
    ggtitle("D")
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
#  MAKE NINE SEPARATE PLOTS (THREE FOR EACH B, C, D)
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

# Combine plots with patchwork with different heights
pBCD<-((p1 / p2 / p3) | (p4 / plot_spacer() / p5 / plot_spacer() / p6 +plot_layout(heights=c(3,1,3,1,3))) | (p7 / p8 / p9)) + plot_layout(widths=c(3,1,1)) 
pBCD

################################################################################
# COMBINE THE TWO PARTS INTO ONE 
final=(pA/pBCD)+plot_layout(heights=c(3,3))
final
outfile=paste(plotdir,"Fig3.svg", sep="")
outfile=paste(plotdir,"Fig3.png", sep="")
ggsave(outfile,plot = final,scale = 1,dpi = 600,limitsize = TRUE,width=11,height=13)

