################################################################################
### R code for plotting non-B motif enrichment in repeats in Zebra finch.
### written by Linnéa Smeds 1-July-2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
plotdir="plots/"

setwd("/Users/lbs5874/Documents/Projects/ZebraFinch/")

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

#DATA<-DATA %>% mutate(PrintName=fct_reorder(PrintName,Repeat))

################################################################################
# SATELLITES ONLY

SUB2<-DATA %>% filter(Class=="Satellites" & Repeat!="Satellite/Satellite" & Group!="GW") %>%
mutate(PrintName=paste(Repeat," (",sprintf("%.0f",Length/1000),"kb)", sep=""))

scalestart2=min(SUB2 %>% filter(Enrichment_gw>0) %>%select(logval))
scaleend2=max(SUB2 %>% filter(Enrichment_gw>0) %>%select(logval))
scalebreak2=-1*scalestart2/(scaleend2-scalestart2)

p2<-ggplot(SUB2, aes(reorder(Repeat, -Length), forcats::fct_rev(nonB))) +
  geom_tile(aes(fill=logval), color = "white",lwd = 0.3,linetype = 1) +
  coord_fixed(ratio = 1)+
  facet_wrap(Classification~., nrow=3, strip.position="right") +
  scale_fill_gradientn(name="log2(Enrichment)",
                       colors=c("darkblue", "white","red"),
                     #  colors=c("#075AFF", "white","#FF0000"),
                      #option="magma",
                       values=c(0, scalebreak2, 1), na.value="gray")+
  labs(x="", y="") +
  theme(axis.text.x = element_markdown(angle = 55, vjust = 1, hjust=1, size=8),
        #   plot.margin = margin(t=0,r=0,b=0,l=0, unit="pt"),
        legend.position="bottom",
        #  legend.justification="top",
        # legend.justification = c(0, 0.7),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10, hjust=0.5),
        legend.title.position="top",
        legend.key.width  = unit(3, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.margin = margin(t=0,r=0,b=0,l=0, unit="pt"),
        strip.background = element_rect(fill = 'white', colour="white"),
        plot.title = element_text(size=18),
        strip.text = element_text(colour = 'black', size=11),
        panel.background = element_rect(fill = 'white', colour="white"),
        panel.spacing.y=unit(1, "lines"),
        )+
  scale_y_discrete(position="left")
p2

outfile2=paste(plotdir,"FigS12_NonB_enrichment_Satellites_log2_facet.png", sep="")
ggsave(outfile2,plot = p2,scale = 1,dpi = 600,limitsize = TRUE,width=9,height=7)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CORRELATION BETWEEN TR SIZE AND NON-B ENRICHMENT
vcolors=viridis(8)
SUB3<-DATA %>% filter(Class=="TRs" & Group!="GW")

p3<-ggplot(SUB3, aes(x=reorder(Repeat, Order), y=Enrichment_gw, color=nonB, fill=nonB, shape=nonB, group=nonB)) +
  geom_point(show.legend=TRUE, alpha=0.6, size=3)+
  geom_line(show.legend=TRUE)+
  scale_fill_manual(values=c(vcolors,"gray"))+
  scale_color_manual(values=c(vcolors, "gray"))+
  scale_shape_manual(values=c(21,22,23,24,25,21,22,23,24))+
  labs(x="Repeat", y="Fold enrichment", scale="free_y") +
  facet_wrap(Classification~., nrow=7, strip.position="right", scales="free_y")+
  scale_x_discrete(expand=c(0.01,0.01))+
  theme(axis.text.x = element_markdown(angle = 55, vjust = 1, hjust=1, size=8),
      legend.title=element_blank(),
      legend.key=element_blank(),
      strip.background = element_rect(fill = 'white', colour="white"),
      strip.text = element_text(colour = 'black', size=11),
      panel.background = element_rect(fill = 'white', colour="black"),
      panel.spacing.y=unit(1, "lines"),
)
p3


outfile3=paste(plotdir,"FigS10_NonB_enrichment_vs_TandemRepeatUnit.png", sep="")
ggsave(outfile3,plot = p3,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=7)



p4<-ggplot(SUB3, aes(reorder(Repeat, Order), Density, color=nonB, fill=nonB)) +
  geom_bar(stat="identity", position=position_dodge(), alpha=0.6) +
  scale_fill_manual(values=c(vcolors,"gray"))+
  scale_color_manual(values=c(vcolors, "gray"))+
  labs(x="Repeat", y="Non-B Density", scale="free_y") +
  facet_wrap(Classification~., nrow=7, strip.position="right")+
  theme(axis.text.x = element_markdown(angle = 55, vjust = 1, hjust=1, size=8),
        legend.title=element_blank(),
        legend.key=element_blank(),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(colour = 'black', size=11),
        panel.background = element_rect(fill = 'white', colour="black"),
        panel.spacing.y=unit(1, "lines"),
  )
p4


