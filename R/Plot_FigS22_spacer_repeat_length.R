################################################################################
### R code for plotting repeat length and spacer distribution for non-B motifs 
### written by Linnéa Smeds 09-December-2025

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(ggplot2)
require(patchwork)
require(viridis)
plotdir="plots/"

setwd("/Users/lbs5874/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Documents/Projects/ZebraFinch")


# Input files:
file1="stats/bTaeGut7v0.4_MT_rDNA_DR.rep-space.txt"
file2="stats/bTaeGut7v0.4_MT_rDNA_IR.rep-space.txt"
file3="stats/bTaeGut7v0.4_MT_rDNA_MR.rep-space.txt"

# read in data and make tibbles 
tib1<-file1 %>% read.table(header=TRUE) %>% as_tibble()
tib2<-file2 %>% read.table(header=TRUE) %>% as_tibble()
tib3<-file3 %>% read.table(header=TRUE) %>% as_tibble()

test2<-tib2 %>% filter(Spacer<11 & Repeat<31)


vcolors=viridis(8)

################################################################################
# Create function that plots a scatter plot with repeat length on one axis and 
# spacer length on the other 

make_plot <- function(data, header, c, verti=NULL) { 
 
  p_main <- ggplot(data, aes(Repeat, Spacer)) +
    geom_point(alpha = 0.1, size = 2, fill=c, color=c)+
    geom_hline(yintercept=10, color="darkred", linetype="dashed")+
    geom_vline(xintercept=verti, color="darkred", linetype="dashed")+
    theme_minimal()
  
  
  # Top histogram (Repeat)
  p_top <- ggplot(data, aes(Repeat)) +
    geom_histogram(binwidth = 1, alpha=0.6, fill=c, color=c) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank())+
    ggtitle(header)
  
  # Right histogram (Spacer)
  p_right <- ggplot(data, aes(Spacer)) +
    geom_histogram(binwidth = 1, alpha=0.6, fill=c, color=c) +
    coord_flip() +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank())
  
  # Combine using patchwork layout
  pcomb<-p_top +plot_spacer() + p_main + p_right +plot_layout(ncol=2, nrow=2, heights=c(1,3), widths=c(3,1))
  return(pcomb)
} 

p1<-make_plot(tib1, "Direct Repeats", vcolors[2])
p1
p2<-make_plot(tib2, "Inverted Repeats", vcolors[4], verti=30)
p2
p3<-make_plot(tib3, "Mirror Repeats", vcolors[5])
p3

p123<-(wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3)) + plot_layout(height=c(1,1,1))
p123

outfile=paste(plotdir,"FigS22_Spacer_vs_repeat.png", sep="")
ggsave(outfile,plot = p123,scale = 1,dpi = 300,limitsize = TRUE,width=5,height=12)




#outfile1=paste(plotdir,"Spacer_vs_repeat_DR.png", sep="")
#ggsave(outfile1,plot = p1,scale = 1,dpi = 300,limitsize = TRUE,width=5,height=5)
#outfile2=paste(plotdir,"Spacer_vs_repeat_IR.png", sep="")
#ggsave(outfile2,plot = p2,scale = 1,dpi = 300,limitsize = TRUE,width=5,height=5)
#outfile3=paste(plotdir,"Spacer_vs_repeat_MR.png", sep="")
#ggsave(outfile3,plot = p3,scale = 1,dpi = 300,limitsize = TRUE,width=5,height=5)








# Old version, using ggMarginal (works, but not so flexible)
require(ggExtra)
make_plot <- function(data, header) { 
  p <- ggplot(data, aes(Repeat, Spacer)) +
    geom_point(alpha = 0.1, size = 0.5) +  # lower opacity for big datasets
    theme_minimal()
  pcomb<- ggMarginal(
    p,
    type = "histogram",
    bins = 50,
    color="#6F9DD0",
    fill="#6F9DD0"
  )
  return(pcomb)
} 

