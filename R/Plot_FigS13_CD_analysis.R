# Script written by Jacob Sieg and modified by Linnéa Smeds.

library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
setwd("/Users/lbs5874/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Documents/Projects/ZebraFinch/")

####Define factor levels and palettes####

Buffer_levels = c("100 mM KCl 140 mM LiCl 20 mM MOPS pH 7.2", "140 mM LiCl 20 mM MOPS pH 7.2") #Intentionally hardcoded
Buffer_labels = c("100 mM KCl 140 mM LiCl", "140 mM LiCl")
Buffer_palette =c("#DB5829", "#F4A637")
Sample_levels = c("Blank", "Tgut368A", "CR1_LINE", "CR1_A", "CR1_B") #Intentionally hardcoded
Sample_labels = c("Blank", "Tgut368A", "CR1 LINE", "CR1 A", "CR1 B")

Titles= c("Blank", "Tgut368A", "LINE/CR1", "TR/Promoter", "LINE/CR1") 

####Read in data####

df_CD = read.csv("js3025_CD_of_Linneas_bird_satellite_pG4s/js3025_calculations_spreadsheet - Preprocessed_CD_Data.csv")
df_Ab = read.csv("js3025_CD_of_Linneas_bird_satellite_pG4s/js3025_calculations_spreadsheet - Preprocessed_UV-Absorbance_Data.csv")

####Check variables and set factor levels####

head(df_CD)
unique(df_CD$Buffer)
#[1] "100 mM KCl 140 mM LiCl 20 mM MOPS pH 7.2" "140 mM LiCl 20 mM MOPS pH 7.2"
unique(df_Ab$Buffer)
#[1] "100 mM KCl 140 mM LiCl 20 mM MOPS pH 7.2" "140 mM LiCl 20 mM MOPS pH 7.2"
unique(df_CD$Sample)
#[1] "Blank"    "CR1_A"    "CR1_B"    "CR1_LINE" "Tgut368A"
unique(df_Ab$Sample)
#[1] "Blank"    "CR1_A"    "CR1_B"    "CR1_LINE" "Tgut368A"
range(df_CD$Wavelength_nm)
#[1] 220 350
range(df_Ab$Wavelength)
#[1] 200 350

df_CD$Buffer = factor(df_CD$Buffer, levels = Buffer_levels, labels = Buffer_labels)
df_Ab$Buffer = factor(df_Ab$Buffer, levels = Buffer_levels, labels = Buffer_labels)

####Function that plots data for each Sample####

plot_data = function(x = "CR1_A", panelnames){
  df_CD_x = df_CD %>% filter(Sample == x)
  df_Ab_x = df_Ab %>% filter(Sample == x)
  
  P_CD = ggplot(df_CD_x, aes(x = Wavelength_nm, y = Molar_ellipticity, color = Buffer, linetype = Buffer)) +
    geom_line() +
    theme_classic() +
    scale_color_manual(values = Buffer_palette) +
    scale_x_continuous(limits = CD_Wavelength_limits, breaks = seq(CD_Wavelength_limits[1], CD_Wavelength_limits[2], by = 10)) +
    theme(axis.text = element_text(size = 6, color = "black"),
          axis.title = element_text(size = 6, color = "black"),
          legend.text = element_text(size = 6, color = "black"),
          legend.title = element_text(size = 6, color = "black"),
          legend.position = "bottom") +
    xlab("Wavelength (nm)") +
    ylab("Molar ellipticity")
  
  P_Tdiff = ggplot(df_Ab_x, aes(x = Wavelength, y = Tdiff, color = Temperature, group = Temperature)) +
    facet_wrap(~Buffer) +
    geom_vline(xintercept = 295) +
    geom_line() +
    theme_classic() +
    scale_color_viridis(option = "turbo") +
    scale_x_continuous(limits = Wavelength_limits, breaks = seq(Wavelength_limits[1], Wavelength_limits[2], by = 10)) +
    scale_y_continuous(limits = Delta_absorbance_limits) +
    theme(axis.text = element_text(size = 6, color = "black"),
          axis.title = element_text(size = 6, color = "black"),
          legend.text = element_text(size = 6, color = "black"),
          legend.title = element_text(size = 6, color = "black"),
          strip.text = element_text(size = 6, color = "black"),
          legend.position = "bottom") +
    xlab("Wavelength (nm)") +
    ylab("delta absorbance")
  
  P_A260 = ggplot(df_Ab_x %>% filter(Wavelength == 260), aes(x =Temperature, y = Aborbtivity, color = Buffer, linetype = Buffer)) +
    geom_line() +
    theme_classic() +
    scale_color_manual(values = Buffer_palette) +
    scale_x_continuous(limits = Temperature_limits, breaks = seq(Temperature_limits[1]+5, Temperature_limits[2]-5, by = 10)) +
    theme(axis.text = element_text(size = 6, color = "black"),
          axis.title = element_text(size = 6, color = "black"),
          legend.text = element_text(size = 6, color = "black"),
          legend.title = element_text(size = 6, color = "black"),
          legend.position = "bottom") +
    xlab("Temperature") +
    ylab("Molar absorbtivity (260 nm)")
  
  P_A294 = ggplot(df_Ab_x %>% filter(Wavelength == 294), aes(x =Temperature, y = Aborbtivity, color = Buffer, linetype = Buffer)) +
    geom_line() +
    theme_classic() +
    scale_color_manual(values = Buffer_palette) +
    scale_x_continuous(limits = Temperature_limits, breaks = seq(Temperature_limits[1]+5, Temperature_limits[2]-5, by = 10)) +
    theme(axis.text = element_text(size = 6, color = "black"),
          axis.title = element_text(size = 6, color = "black"),
          legend.text = element_text(size = 6, color = "black"),
          legend.title = element_text(size = 6, color = "black"),
          legend.position = "bottom") +
    xlab("Temperature") +
    ylab("Molar absorbtivity (294 nm)")
  
  
    
  P_A_B = plot_grid(P_CD, P_Tdiff, rel_widths = c(1, 1.5), labels = c(panelnames[1], panelnames[2]), label_size = 14)
  P_c_d = plot_grid(P_A260, P_A294, rel_widths = c(1, 1), labels = c(panelnames[3], panelnames[4]), label_size = 14)
  P_out = plot_grid(P_A_B, P_c_d, ncol = 1)
  return(P_out)
}

####Plot and save data####

CD_Wavelength_limits = c(220, 320)
Wavelength_limits = c(240, 320)
Temperature_limits = c(5, 95)
Delta_absorbance_limits = c(-0.2, 0.2)

#"Tgut368A" "CR1 LINE", "CR1 A", "CR1 B"
P_T = plot_data(x = "Tgut368A", c("A","B","C","D"))
P_L = plot_data(x = "CR1_LINE", c("E","F","G","H"))
P_A = plot_data(x = "CR1_A", c("I","J","K","L"))
P_B = plot_data(x = "CR1_B", c("M","N","O","P"))

ggsave("figures/Figure_js3025.2_Tgut368A_ABCD.png", P_T, width = 7, height =  5, units = "in")
ggsave("figures/Figure_js3025.3_CR1_LINE_EFGH.png", P_L, width = 7, height =  5, units = "in")
ggsave("figures/Figure_js3025.4_CR1_A_IJKL.png", P_A, width = 7, height =  5, units = "in")
ggsave("figures/Figure_js3025.5_CR1_B_MNOP.png", P_B, width = 7, height =  5, units = "in")
