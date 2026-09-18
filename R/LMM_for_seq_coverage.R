################################################################################
### Linear Mixed Model to investigate how different sequencing content affect
### PacBio and ONT coverage. 
# Written by Linnéa Smeds with input from Claude Sonnet 5 (Anthropic)

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(lme4)
require(lmerTest) 
library(performance)

# Seq method
#method="hifi"
method="ont"

# Input files 
infile=paste("coverage/bTaeGut7v0.4_MT_rDNA.mergedInfo.",method,".bed", sep="")
groupfile="helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt"   

# Reading in the data 
tabtib<-infile %>% read.table(header=TRUE) %>% as_tibble() 
grouptib<-groupfile %>% read.table(header=FALSE) %>% as_tibble() %>%
  select(V1,V3) %>% rename(Chr=V1, Group=V3)

DATA <- tabtib %>% inner_join(grouptib) %>% drop_na()

# -Potentially keep NAs (not used) 
# Keep NA compartment as its own explicit factor level (windows overlapping
# a compartment breakpoint) - if so remove drop.na() above
#DATA$Comp <- addNA(DATA$Comp)         # NA becomes a real factor level, not missing data
#levels(DATA$Comp)[is.na(levels(DATA$Comp))] <- "Unassigned"

DATA$Comp <- factor(DATA$Comp)
DATA$Group <- factor(DATA$Group)
DATA$Chr <- factor(DATA$Chr)


# Fit the mixed model, genome-wide data set
# Random intercept per chromosome accounts for chromosome-level clustering.
full_model <- lmer(
  SeqCov ~ NonBFrac + RepFrac + GCFrac + Comp + Group + (1 | Chr),
  data = DATA,
  REML = TRUE
)
r2_full <- r2(full_model)$R2_marginal

# Check Model summary
summary(full_model)

# Variance explained
# - Marginal R2 (fixed effects only) and conditional R2 (fixed + random)
r2(full_model)

# Type III ANOVA-style test for each fixed effect (with Satterthwaite df)
anova(full_model)

# Save model output 
#saveRDS(model, "full_genome_seqcov_lmm_dropNA.rds")
#capture.output(summary(model), anova(model), r2(model),
#               file = "full_genome_seqcov_lmm_summary_dropNA.txt")


# Get Delta R2 by dropping variables one at the time:
drop_nonb  <- lmer(SeqCov ~           RepFrac + GCFrac + Comp + Group + (1 | Chr), data = DATA, REML = TRUE)
drop_rep   <- lmer(SeqCov ~ NonBFrac +           GCFrac + Comp + Group + (1 | Chr), data = DATA, REML = TRUE)
drop_gc    <- lmer(SeqCov ~ NonBFrac + RepFrac +           Comp + Group + (1 | Chr), data = DATA, REML = TRUE)
drop_comp  <- lmer(SeqCov ~ NonBFrac + RepFrac + GCFrac         + Group + (1 | Chr), data = DATA, REML = TRUE)
drop_group  <- lmer(SeqCov ~ NonBFrac + RepFrac + GCFrac + Comp         + (1 | Chr), data = DATA, REML = TRUE)

delta_r2 <- c(
  NonBFrac = r2_full - r2(drop_nonb)$R2_marginal,
  RepFrac  = r2_full - r2(drop_rep)$R2_marginal,
  GCFrac   = r2_full - r2(drop_gc)$R2_marginal,
  Comp     = r2_full - r2(drop_comp)$R2_marginal,
  Group     = r2_full - r2(drop_group)$R2_marginal
)
delta_r2




################################################################################
# TEST PER GROUP INSTEAD 

# Function that runs the model and saves it
run_lmm_group <- function(data1, prefix) {
  model <- lmer(
    SeqCov ~ NonBFrac + RepFrac + GCFrac + Comp + (1 | Chr),
    data = data1,
    REML = TRUE
  )
  
  model_summary <- summary(model)
  model_r2      <- r2(model)
  model_anova   <- anova(model)
  
  # Save model output
#  saveRDS(model, paste0(prefix,"_seqcov_lmm_dropNA.rds"))
#  capture.output(model_summary, model_anova, model_r2,
#                 file = paste0(prefix,"_seqcov_lmm_summary_dropNA.txt"))
  
  # Drop one at the time 
  r2_full <- r2(model)$R2_marginal
  
  drop_nonb  <- lmer(SeqCov ~           RepFrac + GCFrac + Comp + (1 | Chr), data = data1, REML = TRUE)
  drop_rep   <- lmer(SeqCov ~ NonBFrac +           GCFrac + Comp + (1 | Chr), data = data1, REML = TRUE)
  drop_gc    <- lmer(SeqCov ~ NonBFrac + RepFrac +           Comp + (1 | Chr), data = data1, REML = TRUE)
  drop_comp  <- lmer(SeqCov ~ NonBFrac + RepFrac + GCFrac         + (1 | Chr), data = data1, REML = TRUE)
  
  delta_r2 <- c(
    NonBFrac = r2_full - r2(drop_nonb)$R2_marginal,
    RepFrac  = r2_full - r2(drop_rep)$R2_marginal,
    GCFrac   = r2_full - r2(drop_gc)$R2_marginal,
    Comp     = r2_full - r2(drop_comp)$R2_marginal
  )
  
  return(list(model = model,
              summary = model_summary,
              r2 = model_r2,
              anova = model_anova,
              delta = delta_r2))
}
  
# MAKE GROUP WISE DATA TIBBLES AND RUN MODEL

MACRO<-DATA %>% filter(Group=="macro")
hifi_macro_model<-run_lmm_group(MACRO, "macro")

MICRO<-DATA %>% filter(Group=="micro")
hifi_micro_model<-run_lmm_group(MICRO, "micro")

DOT<-DATA %>% filter(Group=="dot")
hifi_dot_model<-run_lmm_group(DOT, "dot")


################################################################################
# TEST With interaction

data_centred <- DATA %>% mutate(NonBFrac_c = NonBFrac - mean(NonBFrac),
                          GCFrac_c   = GCFrac - mean(GCFrac))
model_int <- lmer(SeqCov ~ NonBFrac_c * GCFrac_c + RepFrac + Comp + Group + (1 | Chr),
                  data = data_centred, REML = TRUE)

summary(model_int)
r2(model_int)
anova(model_int)

# delta (check that you use the correct r2_full)
r2(model_int)$R2_marginal-r2_full 

# Separate groups 
run_lmm_group_interact <- function(data1, prefix) {
  model <- lmer(
    SeqCov ~ NonBFrac * GCFrac + RepFrac + Comp + (1 | Chr),
    data = data1,
    REML = TRUE
  )
  model_summary <- summary(model)
  model_r2      <- r2(model)
  model_anova   <- anova(model)

  return(list(model = model,
              summary = model_summary,
              r2 = model_r2,
              anova = model_anova))
}

hifi_macro_model_interact<-run_lmm_group_interact(MACRO, "macro")
hifi_micro_model_interact<-run_lmm_group_interact(MICRO, "micro")
hifi_dot_model_interact<-run_lmm_group_interact(DOT, "dot")






