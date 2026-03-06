################################################################################
### R code for plotting a comparison between zebra finch and chicken 
### written by Linnéa Smeds 
### September 4, 2025
################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(patchwork) 
require(cowplot)
plotdir="plots/"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIES LIST 
spfile="helpfiles/species_list.txt"
sptib<-spfile %>% read.table(header=FALSE) %>% as_tibble()

# DEFINE COLORS 
type_col=c("#440154","#9460a1", "#e3cfe8")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function for plotting scatterplots

make_scatter_GC<-function(data,hd) {
  p<-ggplot(data, aes(x=GC_cont*100, y=Coverage*100, fill=Group, color=Group))+
    geom_point(aes(shape=Group), alpha=0.6, size=3)+
    labs(x="GC content (%)", y="Non-B motif content (%)") +
    scale_fill_manual(values=type_col)+
    scale_color_manual(values=type_col)+
    scale_shape_manual(values=c(21,24,22))+
    theme(panel.background = element_rect(fill = 'white', colour="black"),
          strip.background = element_rect(fill = 'white', colour="white"),
          strip.text = element_text(size=14),
          plot.margin = margin(b=2,l=2),
          plot.title = element_text(size=18, face="bold"),
    )+ggtitle(hd)
    guides(fill = guide_legend(nrow = 3))
  return(p)
}

make_scatter_len<-function(data,hd) {
  p<-ggplot(data, aes(x=Length/1e6, y=Coverage*100, fill=Group, color=Group))+
    geom_point(aes(shape=Group), alpha=0.6, size=3)+
    labs(x="Chromosome size Mb", y="Non-B motif content (%)") +
    scale_x_log10()+
    scale_fill_manual(values=type_col)+
    scale_color_manual(values=type_col)+
    scale_shape_manual(values=c(21,24,22))+
    theme(panel.background = element_rect(fill = 'white', colour="black"),
          strip.background = element_rect(fill = 'white', colour="white"),
          strip.text = element_text(size=14),
          plot.margin = margin(b=2,l=2),
          plot.title = element_text(size=18, face="bold"),
    )+ggtitle(hd)
  guides(fill = guide_legend(nrow = 3))
  return(p)
}


################################################################################
# PLOT FIG S7, NON-B VS %GC FOR ALL SPECIES 

### PLOT TITLES FOR COMBINED PLOTS
header=c("A", "B", "C", "D", "E", "F", "G", "H")

# LOOP OVER ALL SPECIES AND MAKE PLOTS
plots <- vector("list", length(sptib$V3))
plots_len <- vector("list", length(sptib$V3))
names(plots) <- sptib$V1
for (i in 1:length(sptib$V1)){
  show(paste("Running pipeline for:",sptib$V1[i]))
  
  # READ IN DATA 
  densfile=paste("coverage/",sptib$V3[i],".per_chrom.tsv", sep="")
  gcfile=paste("stats/",sptib$V3[i],".GC_per_chrom.txt", sep="")
  groupfile=paste("helpfiles/",sptib$V3[i],".groups.txt", sep="")
  dens_tib<-densfile %>% read.table(header=TRUE) %>% as_tibble()
  gc_tib<- gcfile %>% read.table(header=TRUE) %>% as_tibble() 
  grouptib<-groupfile %>% read.table(header=FALSE) %>% 
    rename(Chr=V1, Length=V2, Group=V3) %>% 
    mutate(Group=case_when(Group=="macro" ~ "Macro",
                           Group=="micro" ~ "Micro",
                           Group=="dot" ~ "Dot",
                           TRUE ~ Group))
  DATA<-dens_tib %>% inner_join(gc_tib) %>% 
    inner_join(grouptib) %>%  # mutate(Letter=header[i]) %>% 
    filter(NonB=="Any") %>%
    select(Group, GC_cont, Coverage,Length)
  DATA$Group <- factor(DATA$Group, levels=c("Macro", "Micro", "Dot"))
  p1<-make_scatter_GC(DATA, header[i])+
    theme(legend.position = "none",
          axis.title = element_blank())
  p2<-make_scatter_len(DATA, header[i])+
    theme(legend.position = "none",
          axis.title = element_blank())
  plots[[i]] <- p1
  plots_len[[i]]<- p2
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMBINE PLOTS WITH LEGEND 
# Get Legend
legend_plot <- get_legend(
  make_scatter_GC(DATA, header[1]) +
    theme(legend.position = "right",
          legend.background = element_rect(color = NA),
          legend.key = element_rect(color = NA),
          legend.text = element_text(size=12),
          legend.title = element_blank(),
          ))
legend_plot <- wrap_elements(legend_plot)

pgrid1<-(plots[[1]] + plots[[2]]+ plots[[3]]) / (plots[[4]] + plots[[5]]+ plots[[6]]) / (plots[[7]] + plots[[8]] + legend_plot)
pgrid2<-(plots_len[[1]] + plots_len[[2]]+ plots_len[[3]]) / (plots_len[[4]] + plots_len[[5]]+ plots_len[[6]]) / (plots_len[[7]] + plots_len[[8]] + legend_plot)

# Make Y axis 
y_axis <- ggplot() +
  labs(y = "Non-B motif content (%)") +
  theme_void() +
  theme(
    axis.title.y = element_text(angle = 90, size = 14, vjust=0),
    plot.margin = margin(r = 0)
  )

# Make X axes 
x_axis1 <- ggplot() +
  labs(x = "GC content (%)") +
  theme_void() +
  theme(
    axis.title.x = element_text(size = 14, vjust=0),
    plot.margin = margin(t = 0)
  )
x_axis2 <- ggplot() +
  labs(x = "Chromosome size (Mb)") +
  theme_void() +
  theme(
    axis.title.x = element_text(size = 14, vjust=0),
    plot.margin = margin(t = 0)
  )

# Add Y axis to plot, then X axis
# Non B vs GC
p_with_y <- (y_axis | pgrid1) +
  plot_layout(widths = c(0.01, 1))
pcomb1 <- p_with_y /
  x_axis1 +
  plot_layout(heights = c(1, 0.01))
pcomb1 
# NonB vs Len
p_with_y <- (y_axis | pgrid2) +
  plot_layout(widths = c(0.01, 1))
pcomb2 <- p_with_y /
  x_axis2 +
  plot_layout(heights = c(1, 0.01))
pcomb2 

# Save plots
outfile=paste(plotdir,"FigS7_NonB_GC_all_species.png", sep="")
ggsave(outfile,plot = pcomb1,scale = 1,dpi = 600,limitsize = TRUE,width=9,height=10)
outfile=paste(plotdir,"FigS5_NonB_Len_all_species.png", sep="")
ggsave(outfile,plot = pcomb2,scale = 1,dpi = 600,limitsize = TRUE,width=9,height=10)









################################################################################
################################################################################
# OLD FIG S5, COMPARING ONLY ZF AND CHICKEN 

# Extract those two and save in DATA
DATA<-tibble(Group=character(0), GC_cont=numeric(0), Coverage=numeric(0), Length=numeric(0), Species=character(0))
for (i in c(1,6)){
  show(paste("Running pipeline for:",sptib$V1[i]))
  
  # READ IN DATA 
  densfile=paste("coverage/",sptib$V3[i],".per_chrom.tsv", sep="")
  gcfile=paste("stats/",sptib$V3[i],".GC_per_chrom.txt", sep="")
  groupfile=paste("helpfiles/",sptib$V3[i],".groups.txt", sep="")
  dens_tib<-densfile %>% read.table(header=TRUE) %>% as_tibble()
  gc_tib<- gcfile %>% read.table(header=TRUE) %>% as_tibble() 
  grouptib<-groupfile %>% read.table(header=FALSE) %>% 
    rename(Chr=V1, Length=V2, Group=V3) %>% 
    mutate(Group=case_when(Group=="macro" ~ "Macro",
                           Group=="micro" ~ "Micro",
                           Group=="dot" ~ "Dot",
                           TRUE ~ Group))
  tmp<-dens_tib %>% inner_join(gc_tib) %>% 
    inner_join(grouptib) %>%
    filter(NonB=="Any") %>% mutate(Species=sptib$V1[i]) %>%
    select(Group, GC_cont, Coverage, Length, Species)
  DATA<-DATA %>% bind_rows(tmp)
} 
  
  
  
   DATA$Group <- factor(DATA$Group, levels=c("Macro", "Micro", "Dot"))
  
  
  
  p<-make_scatter(DATA, header[i])+
    theme(legend.position = "none",
          axis.title = element_blank())
  plots[[i]] <- p
}


#PANEL A, NON-B VS CHROMOSOME LENGTH  

# Just use "Any"
SUBSET<-DATA %>% filter(NonB=="Any")

pa<-ggplot(SUBSET, aes(x=Length/1000000, y=Coverage*100, fill=Group, color=Group))+
  geom_point(aes(shape=Group), show.legend=FALSE, alpha=0.6, size=3)+
  geom_segment(data=SUBSET%>%filter(Chr=="chrW"), aes(x = Length/1000000+15, y = Coverage*100+4, xend = Length/1000000+3, yend = Coverage*100+1), # Start and end coordinates of the arrow
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), # Adjust arrowhead length
               color = "red", # Set arrow color
               linewidth = 0.2) +
  labs(x="Chromosome Size (Mb)", y="Non-B motif content (%)") +
  scale_fill_manual(values=type_col)+
  scale_color_manual(values=type_col)+
  scale_shape_manual(values=c(21,24,22))+
  scale_x_log10()+
  facet_wrap(Species~.)+
  theme(panel.background = element_rect(fill = 'white', colour="black"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=14),
        panel.grid = element_blank(),
        panel.spacing=unit(1, "lines"),
        legend.background = element_rect(color = NA),
        legend.key = element_rect(color = NA),
        legend.position="right",
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        plot.title = element_text(size=18, face="bold"),
   )+
  ggtitle("A")+
  guides(fill = guide_legend(nrow = 3))
pa

#outfile=paste(plotdir,"FigS3A_NonB_density_vs_ChrLength.png", sep="")
#ggsave(outfile,plot = pc,scale = 1,dpi = 600,limitsize = TRUE,width=9,height=8)

################################################################################
# PLOT PANEL B, PLOT NON-B VS GC CONTENT  

pb<-ggplot(SUBSET, aes(x=GC_cont*100, y=Coverage*100, fill=Group, color=Group))+
  geom_point(aes(shape=Group), show.legend=TRUE, alpha=0.6, size=3)+
  geom_segment(data=SUBSET%>%filter(Chr=="chrW"), aes(x = GC_cont*100-4, y = Coverage*100+5, xend = GC_cont*100-1, yend = Coverage*100+1), # Start and end coordinates of the arrow
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), # Adjust arrowhead length
               color = "red", # Set arrow color
               linewidth = 0.2) +
  labs(x="GC content (%)", y="Non-B motif content (%)") +
  scale_fill_manual(values=type_col)+
  scale_color_manual(values=type_col)+ 
  scale_shape_manual(values=c(21,24,22))+
  facet_wrap(Species~.)+
  theme(panel.background = element_rect(fill = 'white', colour="black"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=14),
        panel.grid = element_blank(),
        panel.spacing=unit(1, "lines"),
        legend.background = element_rect(color = NA),
        legend.key = element_rect(color = NA),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.margin = margin(t=0,l=0,b=0,r=0),
        plot.title = element_text(size=18, face="bold"),
  )+
  guides(fill = guide_legend(nrow = 1))+
  ggtitle("B")
pb

combined_plot <- pa / pb
combined_plot
outfile=paste(plotdir,"FigS3.png", sep="")
ggsave(outfile,plot = combined_plot,scale = 1,dpi = 600,limitsize = TRUE,width=9,height=10)
