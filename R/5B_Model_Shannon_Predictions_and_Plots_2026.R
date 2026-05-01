# ---------------------------------------------------------------------------------------------- #
# Back-transform estimates (Estimate, SE, CI) from logit to proportion scale
# for the Beta submodels (C3, P, Herb, Viv)
# Predict (marginal predictions) from SEM model for Shannon
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(tidyverse); library(purrr); library(brms)

fit_shannon <- readRDS("Results/Shannon_brms_SEM.rds")
fx_shannon  <- as.data.frame(fixef(fit_shannon)) %>% mutate(across(everything(), as.numeric))

inv_logit <- function(x) plogis(as.numeric(x))

backtransform_row <- function(estimate, se, l95, u95) {
  tibble(Estimate_bt = inv_logit(estimate),
         SE_approx   = (inv_logit(estimate + se) - inv_logit(estimate - se)) / 2,
         CI95_low    = inv_logit(l95),
         CI95_high   = inv_logit(u95))
}

extract_backtrans <- function(fx, mediator) {
  rows <- grep(paste0("^", mediator, "_"), rownames(fx), value = TRUE)
  if (length(rows) == 0) return(tibble())
  map_dfr(rows, function(rn) {
    vals  <- as.numeric(fx[rn, c("Estimate", "Est.Error", "Q2.5", "Q97.5")])
    trans <- backtransform_row(vals[1], vals[2], vals[3], vals[4])
    tibble(Mediator    = mediator,
           Parameter   = rn,
           Estimate    = vals[1],
           SE          = vals[2],
           L95         = vals[3],
           U95         = vals[4],
           Estimate_bt = trans$Estimate_bt,
           SE_approx   = trans$SE_approx,
           CI95_low    = trans$CI95_low,
           CI95_high   = trans$CI95_high)
  })
}

mediators_beta <- c("C3beta", "Pbeta", "Herbbeta", "Vivbeta")
bt_beta        <- map_dfr(mediators_beta, ~ extract_backtrans(fx_shannon, .x))

shannon_rows <- grep("^Shannon_", rownames(fx_shannon), value = TRUE)

bt_shannon <- map_dfr(shannon_rows, function(rn) {
  vals <- as.numeric(fx_shannon[rn, c("Estimate", "Est.Error", "Q2.5", "Q97.5")])
  tibble(Mediator    = "Shannon",
         Parameter   = rn,
         Estimate    = vals[1],
         SE          = vals[2],
         L95         = vals[3],
         U95         = vals[4],
         Estimate_bt = vals[1],
         SE_approx   = vals[2],
         CI95_low    = vals[3],
         CI95_high   = vals[4])
})

bt_all <- bind_rows(bt_beta, bt_shannon) %>% arrange(Mediator, Parameter)

write.table(bt_all, "Results/backtransformed_Shannon_2026.txt", sep = ";", row.names = FALSE)

# =====================================================================
# MS plots (standardized): Area -> Shannon, C3 -> Shannon, Viv -> Shannon, Area -> C3
# =====================================================================

library(ggplot2)
library(gridExtra)

# 1) Area (z_Area) -> Shannon
ce_Shan_Area <- conditional_effects(
  fit_shannon,
  effects = "z_Area",
  resp    = "Shannon"
)

pred_Shan_Area <- as.data.frame(ce_Shan_Area[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = z_Area)

p_Shan_Area <- ggplot() +
  geom_jitter(data = fit_shannon$data,
              aes(x = z_Area, y = Shannon),
              size = 2, alpha = 0.5, color = "#A1D99B",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_Shan_Area,
            aes(x = x, y = predicted),
            linewidth = 1.2, color = "#A1D99B") +
  geom_ribbon(data = pred_Shan_Area,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "#e0f3db", alpha = 0.35) +
  labs(x = "Area", y = "Shannon index") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# 2) C3 (zC3_logit) -> Shannon
ce_Shan_C3 <- conditional_effects(
  fit_shannon,
  effects = "zC3_logit",
  resp    = "Shannon"
)

pred_Shan_C3 <- as.data.frame(ce_Shan_C3[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = zC3_logit)

p_Shan_C3 <- ggplot() +
  geom_jitter(data = fit_shannon$data,
              aes(x = zC3_logit, y = Shannon),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_Shan_C3,
            aes(x = x, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_Shan_C3,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  labs(x = "C3-grass", y = "Shannon index") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# 3) Viv (zViv_logit) -> Shannon
ce_Shan_Viv <- conditional_effects(
  fit_shannon,
  effects = "zViv_logit",
  resp    = "Shannon"
)

pred_Shan_Viv <- as.data.frame(ce_Shan_Viv[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = zViv_logit)

p_Shan_Viv <- ggplot() +
  geom_jitter(data = fit_shannon$data,
              aes(x = zViv_logit, y = Shannon),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_Shan_Viv,
            aes(x = x, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_Shan_Viv,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  labs(x = "Viv-grass", y = "Shannon index") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# 4) Area (z_Area) -> C3 proportion (C3beta)
ce_C3_Area_shan <- conditional_effects(
  fit_shannon,
  effects = "z_Area",
  resp    = "C3beta"
)

pred_C3_Area_shan <- as.data.frame(ce_C3_Area_shan[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = z_Area)

p_C3_Area_shan <- ggplot() +
  geom_jitter(data = fit_shannon$data,
              aes(x = z_Area, y = C3_beta),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_C3_Area_shan,
            aes(x = x, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_C3_Area_shan,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  labs(x = "Area", y = "C3-grass") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

grid.arrange(p_Shan_Area, p_C3_Area_shan, p_Shan_C3, p_Shan_Viv, ncol = 2)

# =====================================================================
# SM plots (real scale): Area (ha) -> Shannon, C3 (ha) -> Shannon, Viv (ha) -> Shannon
# =====================================================================

# 1) Back-transform helpers -------------------------------------------

mu_Area   <- attr(scale(SEM_data$Area),             "scaled:center")
sd_Area   <- attr(scale(SEM_data$Area),             "scaled:scale")
mu_C3log  <- attr(scale(qlogis(SEM_data$C3_beta)),  "scaled:center")
sd_C3log  <- attr(scale(qlogis(SEM_data$C3_beta)),  "scaled:scale")
mu_Vivlog <- attr(scale(qlogis(SEM_data$Viv_beta)), "scaled:center")
sd_Vivlog <- attr(scale(qlogis(SEM_data$Viv_beta)), "scaled:scale")

z_to_Area    <- function(z) z * sd_Area   + mu_Area
zC3_to_prop  <- function(z) plogis(z * sd_C3log  + mu_C3log)
zViv_to_prop <- function(z) plogis(z * sd_Vivlog + mu_Vivlog)

# 2) Conditional effects (standardized scale) -------------------------

ce_Shan_Area <- conditional_effects(fit_shannon, effects = "z_Area",    resp = "Shannon")
ce_Shan_C3   <- conditional_effects(fit_shannon, effects = "zC3_logit", resp = "Shannon")
ce_Shan_Viv  <- conditional_effects(fit_shannon, effects = "zViv_logit",resp = "Shannon")

pred_Shan_Area <- as.data.frame(ce_Shan_Area[[1]]) %>%
  rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = z_Area)

pred_Shan_C3 <- as.data.frame(ce_Shan_C3[[1]]) %>%
  rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = zC3_logit)

pred_Shan_Viv <- as.data.frame(ce_Shan_Viv[[1]]) %>%
  rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = zViv_logit)

# 3) Convert predictions + data to real scales ------------------------

pred_Shan_Area_real <- pred_Shan_Area %>%
  mutate(Area_real = z_to_Area(x))

pred_Shan_C3_ha <- pred_Shan_C3 %>%
  mutate(C3_prop = zC3_to_prop(x),
         C3_ha   = C3_prop * mean(SEM_data$Area, na.rm = TRUE))

pred_Shan_Viv_ha <- pred_Shan_Viv %>%
  mutate(Viv_prop = zViv_to_prop(x),
         Viv_ha   = Viv_prop * mean(SEM_data$Area, na.rm = TRUE))

SEM_data <- SEM_data %>%
  mutate(
    Area_real = Area,
    C3_ha     = C3_grass * Area,
    Viv_ha    = V_grass  * Area
  )

# 4) SM plot 1: Area (ha) -> Shannon ---------------------------------

p_Shan_Area_real <- ggplot() +
  geom_jitter(data = SEM_data,
              aes(x = Area_real, y = Shannon),
              size = 2, alpha = 0.5, color = "#A1D99B",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_Shan_Area_real,
            aes(x = Area_real, y = predicted),
            linewidth = 1.2, color = "#A1D99B") +
  geom_ribbon(data = pred_Shan_Area_real,
              aes(x = Area_real, ymin = conf.low, ymax = conf.high),
              fill = "#e0f3db", alpha = 0.35) +
  theme_classic(base_size = 14) +
  labs(x = "Area (ha)", y = "Shannon index")

# 5) SM plot 2: C3 (ha) -> Shannon -----------------------------------

p_Shan_C3_ha <- ggplot() +
  geom_jitter(data = SEM_data,
              aes(x = C3_ha, y = Shannon),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_Shan_C3_ha,
            aes(x = C3_ha, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_Shan_C3_ha,
              aes(x = C3_ha, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  theme_classic(base_size = 14) +
  labs(x = "C3 area (ha)", y = "Shannon index")

# 6) SM plot 3: Viv (ha) -> Shannon (0 hasta max Viv_ha) --------------

max_Viv_ha_data <- max(SEM_data$Viv_ha, na.rm = TRUE)

min_row <- pred_Shan_Viv_ha[which.min(pred_Shan_Viv_ha$Viv_ha), ]
extra0  <- min_row; extra0$Viv_ha <- 0

pred_Shan_Viv_ha_plot <- bind_rows(extra0, pred_Shan_Viv_ha) %>%
  arrange(Viv_ha)

p_Shan_Viv_ha <- ggplot() +
  geom_jitter(data = SEM_data,
              aes(x = Viv_ha, y = Shannon),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_Shan_Viv_ha_plot,
            aes(x = Viv_ha, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_Shan_Viv_ha_plot,
              aes(x = Viv_ha, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  coord_cartesian(xlim = c(0, max_Viv_ha_data), expand = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "Viv area (ha)", y = "Shannon index")

grid.arrange(p_Shan_Area_real, p_Shan_C3_ha, p_Shan_Viv_ha, ncol = 2)
