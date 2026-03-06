################################################################################
### R code for plotting non-B motif enrichment in bird centromeres.
### written by Linnéa Smeds 2-July-2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(ggtext)
plotdir="plots/"

setwd("/Users/lbs5874/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Documents/Projects/ZebraFinch")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input files 
# PANEL A
infile="centromeres/bTaeGut7v0.4_MT_rDNA.enrichment.tsv"
chrfile="centromeres/chr_len.ordered.txt"
bgfile="centromeres/bTaeGut7v0.4_MT_rDNA.background.enrichment.tsv"
cenfile="ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.bed"
# PANEL B
cinfile="centromeres/chicken.v23.enrichment.tsv"
cchrfile="ref/chicken.v23.CEN.bed"
cbgfile="centromeres/chicken.v23.background.enrichment.tsv"
# SUPPLEMENT 
zfile="centromeres/bTaeGut7v0.4_MT_rDNA.ZDNA.enrichment.tsv"
zbgfile="centromeres/bTaeGut7v0.4_MT_rDNA.background.ZDNA.enrichment.tsv"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the data 
# PANEL A & SUPPL FIGURES 13 AND 14, ZF
tib<-infile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Type="centro") %>% filter(NonB!="Any") %>% 
  rename(GW=Enrichment_GW, CHR=Enrichment_Chr) %>% 
  pivot_longer(c(GW,CHR), names_to="Set", values_to="Enrichment")
bgtib<-bgfile %>% read.table(header=TRUE) %>% as_tibble() %>%
  mutate(Type="background") %>% filter(NonB!="Any") %>% 
  rename(GW=Enrichment_GW, CHR=Enrichment_CHR) %>% 
  pivot_longer(c(GW,CHR), names_to="Set", values_to="Enrichment")

centib<-cenfile %>% read.table(header=FALSE) %>% as_tibble() %>% mutate(CenLen=V3-V2)
chrtib<-chrfile %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Chr=V1, Len=V2)
chr_vect <- chrtib$Chr %>% unique()
nonBvect<-c("APR", "DR", "STR", "IR", "TRI", "G4", "Z")

# PANEL B, CHICKEN
ctib<-cinfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Type="centro") %>% filter(NonB!="Any") %>% 
  rename(GW=Enrichment_GW, CHR=Enrichment_Chr) %>% 
  pivot_longer(c(GW,CHR), names_to="Set", values_to="Enrichment")
cbgtib<-cbgfile %>% read.table(header=TRUE) %>% as_tibble() %>%
  mutate(Type="background") %>% filter(NonB!="Any") %>% 
  rename(GW=Enrichment_GW, CHR=Enrichment_CHR) %>% 
  pivot_longer(c(GW,CHR), names_to="Set", values_to="Enrichment")

ccentib<-cchrfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  mutate(CenLen=V3-V2) %>% mutate(Chr=V1)
cchr_vect <- ccentib$Chr %>% unique()

# SUPPLEMENTARY: COMPARE Z-DNA ANNOTATION (NO BACKGROUND)
ztib<-zfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Type="centro") %>%
  rename(GW=Enrichment_GW, CHR=Enrichment_Chr) %>% 
  pivot_longer(c(GW,CHR), names_to="Set", values_to="Enrichment") %>%
  mutate(NonB=case_when(NonB=="Z" ~ "ZDNA Hunter M2",
                        NonB=="ZDNAm1" ~ "ZDNA Hunter M1",
                        NonB=="Zgfa" ~ "GFA ZDNA",
                        NonB=="Zseeker" ~ "ZSeeker",
                        TRUE ~ NonB))
zbgtib<-zbgfile %>% read.table(header=TRUE) %>% as_tibble() %>%
  mutate(Type="background") %>%
  rename(GW=Enrichment_GW, CHR=Enrichment_CHR) %>% 
  pivot_longer(c(GW,CHR), names_to="Set", values_to="Enrichment") %>% 
  mutate(NonB=case_when(NonB=="Z" ~ "ZDNA Hunter M2",
                        NonB=="ZDNAm1" ~ "ZDNA Hunter M1",
                        NonB=="Zgfa" ~ "GFA ZDNA",
                        NonB=="Zseeker" ~ "ZSeeker",
                        TRUE ~ NonB))

Zvect<-c("ZDNA Hunter M1", "ZDNA Hunter M2","GFA ZDNA", "ZSeeker")


################################################################################
### FUNCTIONS 

# Make a function that calculates the significance 
calc_signif <- function(cendata, bgdata) {
  nonB_types=cendata$NonB %>% unique()
  chr_types=cendata$Chr %>% unique()
  
  sign_matrix<- matrix(nrow = 0, ncol = 3)
  colnames(sign_matrix) <- c("Chr", "NonB", "Significance")
  # Go through each cell and check significance 
  for (c in 1:length(chr_types)) {
    for (i in 1:length(nonB_types)) {
      set.seed(1234+i)
      background <- bgdata %>% filter(NonB==nonB_types[i] & Chr==chr_types[c])  %>% 
        select(Enrichment)  %>% pull()
      observed <- cendata %>% filter(NonB==nonB_types[i] & Chr==chr_types[c])  %>% 
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
make_table<-function(data, sign, vect, NonB) {
  DATA<-data %>% inner_join(sign) %>% 
    mutate(Print=if_else(Significance==TRUE, paste(sprintf("%.1f",Enrichment), "*", sep=""), paste(sprintf("%.1f",Enrichment)))) %>%
    mutate(textcol=if_else(Enrichment>0.8*max(data$Enrichment) | Enrichment<0.3,"white","black"))
  DATA$Chr <- factor(DATA$Chr, levels=vect)
  DATA$NonB <- factor(DATA$NonB, levels=NonB)

  return(DATA)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make a function that plots the data horizontally 
make_horizontal<-function(data, header) {
  p<-ggplot(data, aes(y=forcats::fct_rev(NonB), x=Chr, fill=Enrichment))+
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
  pfake<-ggplot(data, aes(y=forcats::fct_rev(NonB), x=Chr, fill=Enrichment))+
    geom_tile(color="white", alpha=0) +
    scale_fill_gradientn(name = "Fold enrichment",
                         colors=c("darkblue", "white","red"),
                         breaks=c(0 ,max(data$Enrichment/2), max(data$Enrichment)),
                         labels = c("0","1",as.character(round(max(data$Enrichment)))))+
    theme(legend.position="right",
          legend.justification="center",
          #legend.margin=margin(c(30,1,0,3)),
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
  p<-ggplot(data, aes(x=NonB, y=forcats::fct_rev(Chr), fill=Enrichment))+
    geom_tile(color="white", show.legend=FALSE) +
    geom_text(aes(label=Print, color=textcol), size=2.5, show.legend=FALSE) +
    labs(x="", y="") +
    scale_fill_gradientn(name = "Enrichment",
                         colors=c("darkblue", "white","red"),
                         values=c(0, 1/max(data$Enrichment), 1),
                         breaks=c(-10 ,0, 1),
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
  p2<-ggplot(data, aes(x=NonB, y=forcats::fct_rev(Chr), fill=Enrichment))+
    geom_tile(color="white", alpha=0) +
    scale_fill_gradientn(name = "Fold enrichment",
                         colors=c("darkblue", "white","red"),
                         breaks=c(0 ,max(data$Enrichment)/2, floor(max(data$Enrichment))),
                         labels = c("0","1",as.character(floor(max(data$Enrichment)))))+
    theme(legend.position="bottom",
        #  legend.margin=margin(c(3,1,0,3)),
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
DATAA<-make_table(tib%>%filter(Set=="GW"),zfgw_sign,chr_vect, nonBvect) %>% 
  filter(str_detect(Chr, "_mat") | str_detect(Chr, "chrZ_pat")) %>% drop_na() 
                  
# PANEL B
chgw_sign<-calc_signif(ctib%>%filter(Set=="GW"),cbgtib)
DATAB<-make_table(ctib%>%filter(Set=="GW"),chgw_sign,cchr_vect, nonBvect)
  
# S15
zfgw_sign<-calc_signif(tib%>%filter(Set=="GW"),bgtib)
DATAS15<-make_table(tib%>%filter(Set=="GW"),zfgw_sign,chr_vect, nonBvect) %>% drop_na() 

# S16
zfchr_sign<-calc_signif(tib%>%filter(Set=="CHR"),bgtib)
DATAS16<-make_table(tib%>%filter(Set=="CHR"),zfchr_sign,chr_vect, nonBvect) %>% drop_na() 

# S17
zdna_zfgw_sign<-calc_signif(ztib%>%filter(Set=="GW"),zbgtib)
DATAS17<-make_table(ztib%>%filter(Set=="GW"),zdna_zfgw_sign,chr_vect,Zvect) %>% drop_na() 


################################################################################
# PLOT 
# FIGURE 4
# PANEL A
pa<-make_horizontal(DATAA, "A")
pa
outfilea=paste(plotdir,"Test_Fig5A.png", sep="")
ggsave(outfilea,plot = pa, bg = "white",scale = 1,dpi = 300,limitsize = TRUE,height=3.9,width=9)

# PANEL B
pb<-make_horizontal(DATAB, "B")
# Combine into one 
pcomb<-cowplot::plot_grid(pa, pb, ncol=1, rel_heights=c(1, 0.9))
pcomb
outfile=paste(plotdir,"Fig5.png", sep="")
ggsave(outfile,plot = pcomb, bg = "white",scale = 1,dpi = 300,limitsize = TRUE,width=9,height=6)
outfile=paste(plotdir,"Fig5.svg", sep="")
ggsave(outfile,plot = pcomb, bg = "white",scale = 1,dpi = 600,limitsize = TRUE,width=9,height=6)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIGURE S15 
ps15<-make_vertical(DATAS15)
outfile2=paste(plotdir,"FigS15_centromere_bothHaplo.png", sep="")
ggsave(outfile2,plot = ps15, bg = "white",scale = 1,dpi = 300,limitsize = TRUE,width=3.9,height=9)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIGURE S16
ps16<-make_vertical(DATAS16)
outfile=paste(plotdir,"FigS16_centromere_vs_chrLevels.png", sep="")
ggsave(outfile,plot = ps16,bg = "white", scale = 1,dpi = 300,limitsize = TRUE,width=3.9,height=9)

  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIGURE S17
psX<-make_vertical(DATASX)
psX
outfile=paste(plotdir,"FigS17_centromere.ZDNA.png", sep="")
ggsave(outfile,plot = ps17,bg = "white", scale = 1,dpi = 300,limitsize = TRUE,width=6,height=9)




################################################################################
# 

# SOME STATS 
sign_tib %>% filter(NonB=="Z") %>% group_by(Significance) %>% summarize(n=n())
# A tibble: 2 × 2
#Significance     n
#<chr>        <int>
#  1 FALSE            5
#2 TRUE            75




################################################################################
