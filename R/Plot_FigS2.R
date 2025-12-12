################################################################################
### R code for plotting non-B motif content vs GC content
### written by Linnéa Smeds 30-June-2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(patchwork) 
library(ggpmisc)
plotdir="plots/"


# Input files 
densfile="densities/bTaeGut7v0.4_MT_rDNA.nonB_per_chrom.txt"
scaffile="ref/bTaeGut7v0.4_MT_rDNA.fa.fai"
gcfile="stats/bTaeGut7v0.4_MT_rDNA.GC_per_chrom.txt"
macrfile="macro.txt"
micrfile="micro.txt"
dotfile="microdot.txt"

# Reading in the data 
dens_tib<-densfile %>% read.table(header=TRUE) %>% as_tibble()
scaf_tib<- scaffile %>% read.table(header=FALSE) %>% as_tibble() %>% 
            select(V1,V2) %>% rename(Chr=V1, Length=V2)
gc_tib<- gcfile %>% read.table(header=TRUE) %>% as_tibble() 
macrtib<-macrfile %>% read.table(header=FALSE) %>% 
  as_tibble() %>% rename(Chr=V1) %>% mutate(Group="Macro")
micrtib<-micrfile %>% read.table(header=FALSE) %>% 
  as_tibble() %>% rename(Chr=V1) %>% mutate(Group="Micro")
dottib<- dotfile %>% read.table(header=FALSE) %>% 
  as_tibble() %>% rename(Chr=V1) %>% mutate(Group="Dot")

typetib<- macrtib %>% bind_rows(micrtib) %>% bind_rows(dottib)

DATA <- dens_tib %>% inner_join(scaf_tib) %>% inner_join(gc_tib) %>% 
    inner_join(typetib)
DATA$NonB <- factor(DATA$NonB, levels=c("ALL","APR", "DR", "STR", "IR", "MR","TRI", "G4", "Z"))
DATA$Group <- factor(DATA$Group, levels=c("Macro", "Micro", "Dot"))

nonBcol=c("darkgray","#AE75A2","#7BB0DF","#008A69","#E9DC6D","#F4A637","darkorange","#DB5829","#894B45")
type_col=c("#440154","#9460a1", "#e3cfe8")
################################################################################
# PLOT Non-B CCONTENT VS GC CONTENT  

p<-ggplot(DATA, aes(x=GC_cont*100, y=Density*100, fill=Group, color=Group))+
  geom_point(aes(shape=Group), show.legend=TRUE, alpha=0.6, size=3)+
  labs(x="GC content (%)", y="Non-B motif content (%)", scale="free_y") +
  scale_fill_manual(values=type_col)+
  scale_color_manual(values=type_col)+
  scale_shape_manual(values=c(21,24,22))+
  facet_wrap(NonB~., ncol=3, scales = "free_y")+
  theme(panel.background = element_rect(fill = 'white', colour="black"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        panel.grid = element_blank(),
        panel.spacing.y=unit(1, "lines"),
        legend.background = element_rect(color = NA),
        legend.position="bottom",
        #legend.justification = "right",
        legend.title = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(size=18),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))+
  guides(fill = guide_legend(nrow = 1))
p

outfile=paste(plotdir,"Fig_S2_NonB_density_vs_GCcontent.png", sep="")
ggsave(outfile,plot = p,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=7)


