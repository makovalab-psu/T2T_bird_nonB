################################################################################
### R code for plotting non-B motif enrichment in bird centromeres.
### written by Linnéa Smeds 2-July-2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(ggtext)
plotdir="plots/"

setwd("/Users/lbs5874/Documents/Projects/ZebraFinch/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input files 
# PANEL A
infile="enrichment/bTaeGut7v0.4_MT_rDNA.centromere.txt"
chrfile="centromeres/chr_len.ordered.txt"
bgfile="centromeres/background_enrichment_merged.tsv"
cenfile="ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.bed"
# PANEL B
cinfile="enrichment/chicken.v23.centromere.new.txt"
cchrfile="ref/chicken.v23.CEN.new.bed"
cbgfile="centromeres/chicken.v23.background_enrichment_merged.new.tsv"
# SUPPL FIGURES
bgchrfile="centromeres/background_enrichmentCHR_merged.tsv"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the data 
# PANEL A & SUPPL FIGURES 13 AND 14, ZF
tib<-infile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Type="centro") %>% filter(nonB!="ALL") %>% 
  rename(GW=Enrichment_GW, CHR=Enrichment_Chr) %>% 
  pivot_longer(c(GW,CHR), names_to="Set", values_to="Enrichment")
bgtib<-bgfile %>% read.table(header=TRUE) %>% as_tibble() %>%
  pivot_longer(c(-Chr,-Window), names_to="nonB", values_to="Enrichment")  %>%
  mutate(Type="background")
bgchrtib<-bgchrfile %>% read.table(header=TRUE) %>% as_tibble() %>%
  pivot_longer(c(-Chr,-Window), names_to="nonB", values_to="Enrichment")  %>%
  mutate(Type="background")
centib<-cenfile %>% read.table(header=FALSE) %>% as_tibble() %>% mutate(CenLen=V3-V2)
chrtib<-chrfile %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Chr=V1, Len=V2)
chr_vect <- chrtib$Chr %>% unique()

# PANEL B, CHICKEN
ctib<-cinfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Type="centro") %>% filter(nonB!="ALL") %>% 
  rename(GW=Enrichment_GW, CHR=Enrichment_Chr) %>% 
  pivot_longer(c(GW,CHR), names_to="Set", values_to="Enrichment")
cbgtib<-cbgfile %>% read.table(header=TRUE) %>% as_tibble() %>%
  pivot_longer(c(-Chr,-Window), names_to="nonB", values_to="Enrichment")  %>%
  mutate(Type="background")

ccentib<-cchrfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  mutate(CenLen=V3-V2) %>% mutate(Chr=V1)
cchr_vect <- ccentib$Chr %>% unique()


################################################################################
### FUNCTIONS 

# Make a function that calculates the significance 
calc_signif <- function(cendata, bgdata) {
  
  nonB_types=cendata$nonB %>% unique()
  chr_types=cendata$Chr %>% unique()
  
  sign_matrix<- matrix(nrow = 0, ncol = 3)
  colnames(sign_matrix) <- c("Chr", "nonB", "Significance")
  # Go through each cell and check significance 
  for (c in 1:length(chr_types)) {
    for (i in 1:length(nonB_types)) {
      set.seed(1234+i)
      background <- bgdata %>% filter(nonB==nonB_types[i] & Chr==chr_types[c])  %>% 
        select(Enrichment)  %>% pull()
      observed <- cendata %>% filter(nonB==nonB_types[i] & Chr==chr_types[c])  %>% 
        select(Enrichment)  %>% pull()
      if(length(background>0) & length(observed>0)) {
        lower_bound <- quantile(background, 0.025)
        upper_bound <- quantile(background, 0.975)
        # Check if the observed value is outside the bounds
        is_significant <- observed < lower_bound | observed > upper_bound
        if(is_significant==TRUE){
          sign=TRUE
        }
        else {
          sign=FALSE
        }
      }
      else {
        sign=NA
      }
      # Add a new row to the matrix
      sign_matrix <- rbind(sign_matrix, c(as.character(chr_types[c]), nonB_types[i], sign) )
    }
  } 
  sign_tib <- sign_matrix %>% as_tibble()
  return(sign_tib)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make a function that merge the data table 
make_table<-function(data, sign, vect) {
  DATA<-data %>% inner_join(sign) %>% 
    mutate(Print=if_else(Significance==TRUE, paste(sprintf("%.1f",Enrichment), "*", sep=""), paste(sprintf("%.1f",Enrichment)))) %>%
    mutate(textcol=if_else(Enrichment>0.8*max(data$Enrichment) | Enrichment<0.3,"white","black"))
  DATA$Chr <- factor(DATA$Chr, levels=vect)
  DATA$nonB <- factor(DATA$nonB, levels=c("APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))
  
  return(DATA)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make a function that plots the data horizontally 
make_horizontal<-function(data, header) {
  p<-ggplot(data, aes(y=forcats::fct_rev(nonB), x=Chr, fill=Enrichment))+
    geom_tile(color="white", show.legend=FALSE) +
    geom_text(aes(label=Print, color=textcol), size=2, show.legend=FALSE) +
    labs(x="", y="") +
    scale_fill_gradientn(name = "Enrichment",
                         colors=c("darkblue", "white","red"),
                         values=c(0, 1/max(data$Enrichment), 1),
                         breaks=c(-10 ,0, 1),
                         labels = c("0","1",as.character(round(max(data$Enrichment)))),
                         na.value="darkblue") +
    scale_color_manual(values=c("black","white")) +
    theme(panel.grid = element_line(colour = 'white'),
          panel.background = element_rect(fill = 'white', colour = 'white'),
          panel.spacing = unit(0.5, "lines"),
          plot.title=element_text(size=18, face="bold"),
         # legend.position="bottom",
          #legend=element_blank(),
          axis.ticks=element_blank(), 
          #legend.margin=margin(c(0,0,0,0)),
          strip.background = element_rect(fill='white'),
          axis.text.x = element_markdown(vjust=1,hjust=1, size=9, angle=55),
    )+ggtitle(header)
  #Fake plot to make legend in log scale 
  pfake<-ggplot(data, aes(y=forcats::fct_rev(nonB), x=Chr, fill=Enrichment))+
    geom_tile(color="white", alpha=0) +
    scale_fill_gradientn(name = "Fold enrichment",
                         colors=c("darkblue", "white","red"),
                         breaks=c(0 ,max(data$Enrichment/2), max(data$Enrichment)),
                         labels = c("0","1",as.character(round(max(data$Enrichment)))))+
    theme(legend.position="right",
          legend.justification="top",
          legend.margin=margin(c(30,1,0,3)),
          legend.key.height = unit(0.85, 'cm'),
          legend.key.width = unit(0.2, 'cm'),
          legend.title.position = "left",
          legend.title = element_text(hjust=0.5, angle=90),
          strip.background = element_rect(fill='white'))

  leg <- ggpubr::get_legend(pfake, position="right")
  
  pcomb<-cowplot::plot_grid(p, leg, ncol=2,
                            rel_widths=c(1, 0.09))
  return(pcomb)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make a function that plots the data vertically 

make_vertical<-function(data) {
  p<-ggplot(data, aes(x=nonB, y=forcats::fct_rev(Chr), fill=Enrichment))+
    geom_tile(color="white", show.legend=FALSE) +
    geom_text(aes(label=Print, color=textcol), size=2.5, show.legend=FALSE) +
    labs(x="", y="") +
    scale_fill_gradientn(name = "Enrichment",
                         colors=c("darkblue", "white","red"),
                         values=c(0, 1/max(data$Enrichment), 1),
                         breaks=c(-10 ,0, 1),
                         #labels = c("0","1","10"),
                         na.value="darkblue") +
    scale_color_manual(values=c("black","white")) +
    theme(panel.grid = element_line(colour = 'white'),
          panel.background = element_rect(fill = 'white', colour = 'white'),
          panel.spacing = unit(0.5, "lines"),
        #  legend.position="bottom",
         # legend.title=element_blank(),
          axis.ticks=element_blank(), 
       #   legend.margin=margin(c(0,0,0,0)),
          strip.background = element_rect(fill='white'),
          strip.text = element_markdown(hjust=0, size=12))
  p
  #Fake plot to make legend in log scale 
  p2<-ggplot(data, aes(x=nonB, y=forcats::fct_rev(Chr), fill=Enrichment))+
    geom_tile(color="white", alpha=0) +
    scale_fill_gradientn(name = "Fold enrichment",
                         colors=c("darkblue", "white","red"),
                         breaks=c(0 ,max(data$Enrichment)/2, round(max(data$Enrichment))),
                         labels = c("0","1",as.character(round(max(data$Enrichment)))))+
    theme(legend.position="bottom",
          legend.margin=margin(c(3,1,0,3)),
          legend.key.height = unit(0.2, 'cm'),
          legend.key.width = unit(1.5, 'cm'),
          legend.title.position = "top",
          legend.title = element_text(hjust=0.5),
          strip.background = element_rect(fill='white'))
  p2
  leg <- ggpubr::get_legend(p2, position="bottom")
  # Combine 
  pcomb<-cowplot::plot_grid(p, leg, ncol=1,
                            rel_heights=c(1, 0.07))
  return(pcomb)
}



################################################################################
# Make significance matrices and add to the data 
# PANEL A 
zfgw_sign<-calc_signif(tib%>%filter(Set=="GW"),bgtib)
DATAA<-make_table(tib%>%filter(Set=="GW"),zfgw_sign,chr_vect) %>% 
  filter(str_detect(Chr, "_mat") | str_detect(Chr, "chrZ_pat")) 
                  
# PANEL B
chgw_sign<-calc_signif(ctib%>%filter(Set=="GW"),cbgtib)
DATAB<-make_table(ctib%>%filter(Set=="GW"),chgw_sign,cchr_vect)
  
# S13 
zfgw_sign<-calc_signif(tib%>%filter(Set=="GW"),bgtib)
DATAS13<-make_table(tib%>%filter(Set=="GW"),zfgw_sign,chr_vect)

# S14
zfchr_sign<-calc_signif(tib%>%filter(Set=="CHR"),bgchrtib)
DATAS14<-make_table(tib%>%filter(Set=="CHR"),zfchr_sign,chr_vect)

################################################################################
# PLOT 
# FIGURE 4
# PANEL A
pa<-make_horizontal(DATAA, "A")
# PANEL B
pb<-make_horizontal(DATAB, "B")
# Combine into one 
pcomb<-cowplot::plot_grid(pa, pb, ncol=1, rel_heights=c(1, 0.9))
outfile=paste(plotdir,"Fig4.png", sep="")
ggsave(outfile,plot = pcomb, bg = "white",scale = 1,dpi = 300,limitsize = TRUE,width=9,height=6)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIGURE S13 
ps13<-make_vertical(DATAS13)
outfile2=paste(plotdir,"FigS13.png", sep="")
ggsave(outfile2,plot = pcomb, bg = "white",scale = 1,dpi = 300,limitsize = TRUE,width=3.9,height=9)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIGURE S14
ps14<-make_vertical(DATAS14)

outfile=paste(plotdir,"FigS14.png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 300,limitsize = TRUE,width=3.9,height=9)

  
################################################################################
# 

# SOME STATS 
sign_tib %>% filter(nonB=="Z") %>% group_by(Significance) %>% summarize(n=n())
# A tibble: 2 × 2
#Significance     n
#<chr>        <int>
#  1 FALSE            5
#2 TRUE            75




################################################################################
