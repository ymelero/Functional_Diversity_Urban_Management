# ---------------------------------------------------------------------------------------------- #
# Bayesian SEM  for RaoQ, Shannon, Faith's PD
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)
library(car)   

# -----------------------------
# Load and prepare data
# -----------------------------
space <- read.delim("DATA/space_data.txt", sep = "\t")
mint  <- read.delim("DATA/MIntensity.txt", sep = "\t")
names(mint)[4] <- "Intensity"   # columna con low/medium/high

mint <- mint %>%
  filter(Vegetation == "Herb") %>%
  select(,-c(2,3))   # ajusta si los nombres difieren

space <- space %>%
  left_join(mint, by = "Type")

fd.df <- read.delim("DATA/fd_df.txt", sep = ";")

space <- space %>%
  mutate(
    C3.proportion   = C3ha  / TotalA,
    P.proportion    = PUha  / TotalA,
    Herb.proportion = pmin(HERbha / TotalA, 1),
    Viv.proportion  = VIVha / TotalA
  )

SEM_data <- space %>%
  select(SITE_ID, Area = TotalA, Conn = C500, Intensity = Intensity,
         C3_grass = C3.proportion, P_grass = P.proportion,
         Herb_grass = Herb.proportion, V_grass = Viv.proportion)

# -----------------------------
# Smithson–Verkuilen + estandarización (como hurdle)
# -----------------------------
n_sites <- nrow(SEM_data)
SEM_data <- SEM_data %>%
  mutate(
    C3_beta   = (pmin(pmax(C3_grass ,0),1) * (n_sites - 1) + 0.5) / n_sites,
    P_beta    = (pmin(pmax(P_grass  ,0),1) * (n_sites - 1) + 0.5) / n_sites,
    Herb_beta = (pmin(pmax(Herb_grass,0),1) * (n_sites - 1) + 0.5) / n_sites,
    Viv_beta  = (pmin(pmax(V_grass  ,0),1) * (n_sites - 1) + 0.5) / n_sites,
    z_Area = as.numeric(scale(Area)),
    z_Conn = as.numeric(scale(Conn)),
    zC3_logit   = as.numeric(scale(qlogis(C3_beta))),
    zP_logit    = as.numeric(scale(qlogis(P_beta))),
    zHerb_logit = as.numeric(scale(qlogis(Herb_beta))),
    zViv_logit  = as.numeric(scale(qlogis(Viv_beta)))
  )

# Add RaoQ, Shannon and PD (SITE_ID x YEAR)
SEM_data <- SEM_data %>%
  right_join(fd.df, by = "SITE_ID")

# -----------------------------
# Multicolinealidad (VIF)
# -----------------------------
library(performance)
lm_rao <- lm(RaoQ ~ z_Area + z_Conn +
               zC3_logit + zP_logit + zHerb_logit + zViv_logit +
               Intensity,data = SEM_data)
check_collinearity(lm_rao)

lm_rao <- lm(RaoQ ~ z_Area + z_Conn +
               zC3_logit + zP_logit + zHerb_logit + zViv_logit +
               Intensity,data = SEM_data)
check_collinearity(lm_rao)

lm_Shannon <- lm(Shannon ~ z_Area + z_Conn +
               zC3_logit + zP_logit + zHerb_logit + zViv_logit +
               Intensity,data = SEM_data)
check_collinearity(lm_Shannon)

lm_PD <- lm(PD ~ z_Area + z_Conn +
                   zC3_logit + zP_logit + zHerb_logit + zViv_logit +
                   Intensity,data = SEM_data)
check_collinearity(lm_PD)

# -----------------------------
# Model formulas: mediators
# -----------------------------
bf_C3   <- bf(C3_beta   ~ z_Area, family = Beta(link = "logit"))
bf_P    <- bf(P_beta    ~ z_Area + Intensity, family = Beta(link = "logit"))
bf_Herb <- bf(Herb_beta ~ z_Area + Intensity, family = Beta(link = "logit"))
bf_Viv  <- bf(Viv_beta  ~ z_Area + Intensity, family = Beta(link = "logit"))

# -----------------------------
# SEM 1: RaoQ
# -----------------------------
bf_Rao  <- bf(
  RaoQ ~ z_Area + z_Conn + Intensity +
    zC3_logit + zP_logit + zHerb_logit + zViv_logit +
    (1 | YEAR),
  family = gaussian()
)

fit_Rao <- brm(
  bf_C3 + bf_P + bf_Herb + bf_Viv + bf_Rao + set_rescor(FALSE),
  data = SEM_data,
  chains = 4, iter = 120000, warmup = 20000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = 0.999, max_treedepth = 15),
  seed = 20251004
)

saveRDS(fit_Rao, file = "Results/Rao_brms_SEM_2026.rds")

summary(fit_Rao)
range(neff_ratio(fit_Rao), na.rm = TRUE)
pp_check(fit_Rao, resp = "RaoQ")

ce_rao <- conditional_effects(fit_Rao, resp = "RaoQ")
plot(ce_rao, points = TRUE)

rm(fit_Rao,ce_rao)
gc()

# -----------------------------
# SEM 2: Shannon
# -----------------------------
bf_Shannon <- bf(
  Shannon ~ z_Area + z_Conn + Intensity +
    zC3_logit + zP_logit + zHerb_logit + zViv_logit +
    (1 | YEAR),
  family = gaussian()
)

fit_Shannon <- brm(
  bf_C3 + bf_P + bf_Herb + bf_Viv + bf_Shannon + set_rescor(FALSE),
  data = SEM_data,
  chains = 4, iter = 120000, warmup = 20000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = 0.999, max_treedepth = 15),
  seed = 20251005
)

saveRDS(fit_Shannon, "Results/Shannon_brms_SEM.rds")

summary(fit_Shannon)
range(neff_ratio(fit_Shannon), na.rm = TRUE)
pp_check(fit_Shannon, resp = "Shannon")

ce_shannon <- conditional_effects(fit_Shannon, resp = "Shannon")
plot(ce_shannon, points = TRUE)

rm(fit_Shannon,ce_shannon)
gc()

# -----------------------------
# SEM 3: Faith's PD
# -----------------------------
bf_PD <- bf(
  PD ~ z_Area + z_Conn + Intensity +
    zC3_logit + zP_logit + zHerb_logit + zViv_logit +
    (1 | YEAR),
  family = gaussian()
)

fit_PD <- brm(
  bf_C3 + bf_P + bf_Herb + bf_Viv + bf_PD + set_rescor(FALSE),
  data = SEM_data,
  chains = 4, iter = 120000, warmup = 20000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = 0.999, max_treedepth = 15),
  seed = 20251006
)

saveRDS(fit_PD, "Results/PD_brms_SEM.rds")

summary(fit_PD)
range(neff_ratio(fit_PD), na.rm = TRUE)
pp_check(fit_PD, resp = "PD")

ce_pd <- conditional_effects(fit_PD, resp = "PD")
plot(ce_pd, points = TRUE)


# Extract results

# RAOQ
library(posterior)
fit_Rao <- readRDS("Results/RaO_brms_SEM_2026.rds")
post.Rao <- posterior::summarise_draws(fit_Rao) 
post.Rao <- post.Rao %>%
  mutate(
    sig_95 = if_else(q5 > 0 | q95 < 0, "sig", "no_sig")
  ) 
post.Rao <-post.Rao[,c(1,2, 4,6,7,11)]
write.csv(post.Rao, "Results/Results_Rao_2026.csv", row.names = F)

# Shannon
fit_Shanon <- readRDS("Results/Shannon_brms_SEM.rds")
post.Shanon <- posterior::summarise_draws(fit_Shanon)  
post.Shanon <- post.Shanon %>%
  mutate(
    sig_95 = if_else(q5 > 0 | q95 < 0, "sig", "no_sig")
  ) 
post.Shanon <-post.Shanon[,c(1,2, 4,6,7,11)]
write.csv(post.Shanon, "Results/Results_Shanon_2026.csv", row.names = F)

# FPD
fit_PD <- readRDS("Results/PD_brms_SEM.rds")
post.PD <- posterior::summarise_draws(fit_PD)  
post.PD <- post.PD %>%
  mutate(
    sig_95 = if_else(q5 > 0 | q95 < 0,"sig", "no_sig")
  ) 
post.PD <-post.PD[,c(1,2, 4,6,7,11)]
write.csv(post.PD, "Results/Results_PD_2026.csv", row.names = F)

