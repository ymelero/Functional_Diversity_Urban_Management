# ---------------------------------------------------------------------------------------------- #
# Back-transform estimates (Estimate, SE, CI) from logit to proportion scale
# for the Beta submodels (C3, P, Herb, Viv)
# AND from log / logit scale for abundance and hurdle components
# Predict (marginal preditions) from SEM model for RaoQ 
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(tidyverse)
library(purrr)
library(brms)
detach("package:ggeffects")
detach("package:gridExtra")

fit_rao <- readRDS("Results/RaO_brms_SEM_2026.rds")
fx_rao <- as.data.frame(fixef(fit_rao)) %>% mutate(across(everything(), as.numeric))

# 1. Site variables
inv_logit <- function(x) plogis(as.numeric(x))

# Back-transform one row for Beta models (logit → proportion)
backtransform_row <- function(estimate, se, l95, u95) {
  tibble(
    Estimate_bt = inv_logit(estimate),
    SE_approx   = (inv_logit(estimate + se) - inv_logit(estimate - se)) / 2,
    CI95_low    = inv_logit(l95),
    CI95_high   = inv_logit(u95)
  )
}


extract_backtrans <- function(fx, mediator) {
  rows <- grep(paste0("^", mediator, "_"), rownames(fx), value = TRUE)
  
  if (length(rows) == 0) return(tibble())  # seguridad: si no hay coincidencias, devuelve vacío
  
  map_dfr(rows, function(rn) {
    vals <- as.numeric(fx[rn, c("Estimate", "Est.Error", "Q2.5", "Q97.5")])
    trans <- backtransform_row(vals[1], vals[2], vals[3], vals[4])
    tibble(
      Mediator = mediator,
      Parameter = rn,
      Estimate = vals[1],
      SE = vals[2],
      L95 = vals[3],
      U95 = vals[4],
      Estimate_bt = trans$Estimate_bt,
      SE_approx = trans$SE_approx,
      CI95_low = trans$CI95_low,
      CI95_high = trans$CI95_high
    )
  })
}

# Run for all Beta mediators
mediators_beta <- c("C3beta", "Pbeta", "Herbbeta", "Vivbeta")
bt_beta <- map_dfr(mediators_beta, ~ extract_backtrans(fx_rao, .x))

rao_rows <- grep("^RaoQ_", rownames(fx_rao), value = TRUE)

bt_rao <- map_dfr(rao_rows, function(rn) {
  vals <- as.numeric(fx_rao[rn, c("Estimate", "Est.Error", "Q2.5", "Q97.5")])
  tibble(
    Mediator = "RaoQ",
    Parameter = rn,
    Estimate = vals[1],
    SE = vals[2],
    L95 = vals[3],
    U95 = vals[4],
    Estimate_bt = vals[1],   # mismo valor (ya en escala original)
    SE_approx = vals[2],
    CI95_low = vals[3],
    CI95_high = vals[4]
  )
})

# Combine Beta + RaoQ result
bt_all <- bind_rows(bt_beta, bt_rao) %>%
  arrange(Mediator, Parameter)

# --- 4. Save ---
write.table(bt_all, "Results/backtransformed_RAO_2026.txt", sep = ";", row.names = FALSE)

# 3. Predict & Plot
library(ggeffects)
library(gridExtra)

# Pred Area -> C3
ce_C3_Area <- conditional_effects(fit_rao,
                                  effects = "z_Area",
                                  resp = "C3beta")

pred_C3_Area <- as.data.frame(ce_C3_Area[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = z_Area)

# Pred C3 -> RaoQ
ce_Rao_C3 <- conditional_effects(fit_rao,
                                 effects = "zC3_logit",
                                 resp = "RaoQ")

pred_Rao_C3 <- as.data.frame(ce_Rao_C3[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = zC3_logit)

# Plots
p_C3_Area <- ggplot() +
  geom_jitter(data = fit_rao$data, aes(x = z_Area, y = C3_beta),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_C3_Area, aes(x = x, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_C3_Area,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  labs(x = "Area", y = "C3-grass") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")
p_C3_Area

p_rao_C3 <- ggplot() +
  geom_jitter(data = fit_rao$data, aes(x = zC3_logit, y = RaoQ),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.002, width = 0.003) +
  geom_line(data = pred_Rao_C3, aes(x = x, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_Rao_C3,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  labs(x = "C3-grass", y = "RaoQ") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")
p_rao_C3

p_all<-grid.arrange(p_C3_Area, p_rao_C3, ncol =2, nrow = 1)


# =====================================================================
# SM plots: Area -> C3 proportion, C3 proportion -> RaoQ,
#           C3 area (ha) -> RaoQ
# =====================================================================

# SEM_data to be loaded from R file 4_Model_Baysian_SEM_Rao...
# SEM_data has: Area, C3_grass (proportion), C3_beta, z_Area, zC3_logit, RaoQ, etc.

# --- 1. Recover scaling (Area and C3_logit) --------------------------

mu_Area  <- attr(scale(SEM_data$Area),           "scaled:center")
sd_Area  <- attr(scale(SEM_data$Area),           "scaled:scale")
mu_C3log <- attr(scale(qlogis(SEM_data$C3_beta)),"scaled:center")
sd_C3log <- attr(scale(qlogis(SEM_data$C3_beta)),"scaled:scale")

z_to_Area   <- function(z) z * sd_Area  + mu_Area          # z_Area -> Area (ha)
zC3_to_prop <- function(z) plogis(z * sd_C3log + mu_C3log) # zC3_logit -> C3 proportion

# --- 2. Conditional effects (standardized scale) ---------------------

# Area -> C3_beta
ce_C3_Area <- conditional_effects(fit_rao,
                                  effects = "z_Area",
                                  resp    = "C3beta")

pred_C3_Area <- as.data.frame(ce_C3_Area[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = z_Area)

# C3_beta -> RaoQ
ce_Rao_C3 <- conditional_effects(fit_rao,
                                 effects = "zC3_logit",
                                 resp    = "RaoQ")

pred_Rao_C3 <- as.data.frame(ce_Rao_C3[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = zC3_logit)

# --- 3. Convert predictions and data to real scales ------------------

# Real-scale predictions
pred_C3_Area_real <- pred_C3_Area %>%
  mutate(Area_real = z_to_Area(x))   # Area in ha

pred_Rao_C3_prop <- pred_Rao_C3 %>%  # C3 proportion
  mutate(C3_prop = zC3_to_prop(x))

pred_Rao_C3_ha <- pred_Rao_C3_prop %>%   # C3 area in ha for a typical Area
  mutate(C3_ha = C3_prop * mean(SEM_data$Area, na.rm = TRUE))

# Real data with matching scales
data_real <- SEM_data %>%
  mutate(
    Area_real = Area,
    C3_prop   = C3_grass,
    C3_ha     = C3_grass * Area
  )

# --- 4. SM plot 1: Area (ha) -> C3 proportion ------------------------

p_C3_Area_real <- ggplot() +
  geom_jitter(data = data_real,
              aes(x = Area_real, y = C3_prop),
              size   = 2,
              alpha  = 0.5,
              color  = "orange",
              height = 0.002,
              width  = 0.003) +
  geom_line(data = pred_C3_Area_real,
            aes(x = Area_real, y = predicted),
            linewidth = 1.2,
            color     = "peachpuff") +
  geom_ribbon(data = pred_C3_Area_real,
              aes(x = Area_real, ymin = conf.low, ymax = conf.high),
              fill  = "peachpuff",
              alpha = 0.35) +
  theme_classic(base_size = 14) +
  labs(x = "Area (ha)", y = "C3-grass")

# --- 5. SM plot 2: C3 proportion -> RaoQ -----------------------------

p_Rao_C3_prop <- ggplot() +
  geom_jitter(data = data_real,
              aes(x = C3_prop, y = RaoQ),
              size   = 2,
              alpha  = 0.5,
              color  = "orange",
              height = 0.002,
              width  = 0.003) +
  geom_line(data = pred_Rao_C3_prop,
            aes(x = C3_prop, y = predicted),
            linewidth = 1.2,
            color     = "peachpuff") +
  geom_ribbon(data = pred_Rao_C3_prop,
              aes(x = C3_prop, ymin = conf.low, ymax = conf.high),
              fill  = "peachpuff",
              alpha = 0.35) +
  theme_classic(base_size = 14) +
  labs(x = "C3-grass", y = "RaoQ")

# --- 6. SM plot 3: C3 area (ha) -> RaoQ ------------------------------

p_Rao_C3_ha <- ggplot() +
  geom_jitter(data = data_real,
              aes(x = C3_ha, y = RaoQ),
              size   = 2,
              alpha  = 0.5,
              color  = "orange",
              height = 0.002,
              width  = 0.003) +
  geom_line(data = pred_Rao_C3_ha,
            aes(x = C3_ha, y = predicted),
            linewidth = 1.2,
            color     = "peachpuff") +
  geom_ribbon(data = pred_Rao_C3_ha,
              aes(x = C3_ha, ymin = conf.low, ymax = conf.high),
              fill  = "peachpuff",
              alpha = 0.35) +
  theme_classic(base_size = 14) +
  labs(x = "C3 area (ha)", y = "RaoQ")

# Example: show two plots side by side
grid.arrange(p_C3_Area_real, p_Rao_C3_prop, ncol = 2)
# And C3 area (ha) -> RaoQ separately:
p_Rao_C3_ha
grid.text("(a)", x = 0.114, y = 0.97, gp = gpar(fontsize = 12))
