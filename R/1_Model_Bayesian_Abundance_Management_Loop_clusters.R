# ---------------------------------------------------------------------------------------------- #
# Bayesian SEM with brms: Beta mediators + hurdle lognormal abundance
# Loop for Functional Group (cluster)
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)

# -----------------------------
# Load and prepare data
# -----------------------------
# SINDEX_ZEROS contains counts and zeros, for the species x year not seen in a site when seen in another.
ubmsdata <- read.table("DATA/SINDEX_ZEROS.txt", header = TRUE, sep = "\t") %>%
  filter(YEAR > 2020) %>%
  mutate(
    SPECIES = case_when(
      SPECIES == "Glaucopsyche sp." ~ "Glaucopsyche melanops", 
      SPECIES == "Melanargia sp."   ~ "Melanargia lachesis",
      SPECIES == "Satyridae"   ~ "Pyronia sp",
      
      TRUE                          ~ SPECIES
    )
  )

space    <- read.delim("DATA/space_data.txt", sep = "\t")
clusters <- read.csv2("DATA/Species_clusters.csv")
mint  <- read.delim("DATA/MIntensity.txt", sep = "\t")
names(mint)[4] <- "Intensity"

mint <- mint %>%
  filter(Vegetation == "Herb") %>%
  .[ , -c(2,3)]
  
# Merge with space dataset
space <- space %>%
  left_join(mint, by = "Type") 

space <- space %>%
  mutate(
    C3.proportion   = C3ha / TotalA, 
    P.proportion    = PUha / TotalA, 
    Herb.proportion = pmin(HERbha / TotalA, 1), 
    Viv.proportion  = VIVha / TotalA
  )

ubmsdata <- ubmsdata %>%
  left_join(clusters, by = "SPECIES") %>%
  filter(!is.na(Cluster)) %>%
  left_join(space, by = "SITE_ID") %>%
  mutate(
    presence  = as.integer(SINDEX > 0),
    abundance = SINDEX,   # continuous proxy, >=0
    Cluster   = as.factor(Cluster)
  ) %>%
  rename(
    Area       = TotalA,
    Conn       = C500,
    C3_grass   = C3.proportion,
    P_grass    = P.proportion,
    Herb_grass = Herb.proportion,
    V_grass    = Viv.proportion
  )

# -----------------------------
# Standardize predictors + prepare mediators
# -----------------------------
ubmsdata <- ubmsdata %>%
  mutate(
    Intensity = factor(Intensity),
    # Smithson–Verkuilen transformation to keep proportions in ]0,1[ . This is needed for the Beta Family used afterwards
    C3_beta   = (C3_grass   * (n() - 1) + 0.5) / n(),
    P_beta    = (P_grass    * (n() - 1) + 0.5) / n(),
    Herb_beta = (Herb_grass * (n() - 1) + 0.5) / n(),
    Viv_beta  = (V_grass    * (n() - 1) + 0.5) / n(),
    # Standardize exogenous predictors for better improvement in the models and for better relative comparison
    z_Area = as.numeric(scale(Area)),
    z_Conn = as.numeric(scale(Conn)),
    # Logit-transform mediators for use as predictors in abundance model (improves model, better for comparison), then standardize (to improve relative comparison)
    zC3_logit   = as.numeric(scale(qlogis(C3_beta))),
    zP_logit    = as.numeric(scale(qlogis(P_beta))),
    zHerb_logit = as.numeric(scale(qlogis(Herb_beta))),
    zViv_logit  = as.numeric(scale(qlogis(Viv_beta)))
  )


# ------------------------------------------------------------------
# Function to test for multi-collinearity of explanatory variables
# ------------------------------------------------------------------
library(performance)

check_cluster_collinearity <- function(cluster_id, data) {
  df <- data %>% filter(Cluster == cluster_id)
  lm_abund <- lm(
    abundance ~ z_Area + z_Conn + Intensity + zC3_logit + zP_logit + zHerb_logit + zViv_logit,
    data = df
  )
  check_collinearity(lm_abund)
}

cluster <- levels(ubmsdata$Cluster)
collinearity_list <- lapply(cluster, check_cluster_collinearity, data = ubmsdata)
names(collinearity_list) <- cluster


# -----------------------------
# Function to fit SEM per cluster
# -----------------------------
fit_cluster_sem <- function(cluster_id, data) {
  
  message(">> Fitting SEM for cluster: ", cluster_id)
  
  df <- data %>% filter(Cluster == cluster_id)
  
  # --- Model formulas ---
  bf_C3   <- bf(C3_beta   ~ z_Area, family = Beta(link = "logit"))
  bf_P    <- bf(P_beta    ~ z_Area + Intensity, family = Beta(link = "logit"))
  bf_Herb <- bf(Herb_beta ~ z_Area + Intensity, family = Beta(link = "logit"))
  bf_Viv  <- bf(Viv_beta  ~ z_Area + Intensity, family = Beta(link = "logit"))
  
  bf_abund <- bf(
    abundance ~ z_Area + z_Conn + Intensity + zC3_logit + zP_logit + zHerb_logit + zViv_logit + (1 | YEAR) + (1 | SPECIES),
    hu       ~ z_Area + z_Conn + Intensity + zC3_logit + zP_logit + zHerb_logit + zViv_logit,
    family = hurdle_lognormal()
  )
  
  # --- Fit SEM ---
  fit <- brm(
    bf_C3 + bf_P + bf_Herb + bf_Viv + bf_abund,
    data = df,
    chains = 4, iter = 100000, warmup = 10000, cores = parallel::detectCores(),
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    seed = 20251004
  )
  
  # --- Save ---
  saveRDS(fit, file = paste0("Results/Abundance_brms_SEM_Cluster", cluster_id, ".rds"))
  message(">> Model for cluster ", cluster_id, " saved successfully.")
  return(fit)
}

# -----------------------------
# Run for all clusters
# -----------------------------
cluster_levels <- levels(ubmsdata$Cluster)

fits <- purrr::map(cluster_levels, ~ fit_cluster_sem(.x, ubmsdata))
names(fits) <- cluster_levels

# ---------------------------------------------------------------------------------------------- #
# Post-processing Bayesian SEM (hurdle lognormal + beta mediators) per cluster
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

# 1. Load models saved 
fit_c4 <- readRDS("Results/Abundance_brms_SEM_Cluster4.rds")
fit_c2 <- readRDS("Results/Abundance_brms_SEM_Cluster2.rds")
fit_c3 <- readRDS("Results/Abundance_brms_SEM_Cluster3.rds")

# 2. Basic diagnostics 
# Always check Rhat < 1.01 and reasonable neff_ratio,.pairs()!!!!!
summary(fit_c2)                          # detailed overview // check Rhat in summary for model fit Rhat ≈ 1.00 → chains have converged well.
summary(fit_c3)
summary(fit_c4)
# check in summary for R hat for each parameter SHOULD be < 1.05
range(neff_ratio(fit_c2), na.rm = TRUE)  # Effective sample size ratio, should be > 0.1 (better if > 0.2).

# 3. Posterior predictive checks 
# Includes both zeros and positive abundances (hurdle model)
pp_check(fit_c2, resp = "abundance") # includes de hu_ part

# 4. Marginal effects for plotting 
# Abundance conditional on presence
ce_abund_c2 <- conditional_effects(fit_c2, resp = "abundance")
plot(ce_abund_c2, points = TRUE)

# Probability of zero (absence process)
ce_hu_c2    <- conditional_effects(fit_c2, resp = "abundance", dpar = "hu")
plot(ce_hu_c2)

# create database file with estimates and sig /non-sig
fit_c2<- posterior::summarise_draws(fit_c2)  
fit_c2 <- fit_c2 %>%
  mutate(
    sig_95 = if_else(q5 > 0 | q95 < 0,"sig", "no_sig")
  ) 
fit_c2 <-fit_c2[,c(1,2, 4,6,7,11)]
write.csv(fit_c2, "Results/Results_C2_2026.csv", row.names = F)
