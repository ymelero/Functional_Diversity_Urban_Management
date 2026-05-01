# ---------------------------------------------------------------------------------------------- #
# Back-transform estimates (Estimate, SE, CI) from logit to proportion scale
# for the Beta submodels (C3, P, Herb, Viv)
# Predict (marginal predictions) from SEM model for Faith's PD
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(tidyverse); library(purrr); library(brms)

fit_PD <- readRDS("Results/PD_brms_SEM.rds")
fx_PD  <- as.data.frame(fixef(fit_PD)) %>% mutate(across(everything(), as.numeric))

inv_logit <- function(x) plogis(as.numeric(x))

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
  if (length(rows) == 0) return(tibble())
  map_dfr(rows, function(rn) {
    vals  <- as.numeric(fx[rn, c("Estimate", "Est.Error", "Q2.5", "Q97.5")])
    trans <- backtransform_row(vals[1], vals[2], vals[3], vals[4])
    tibble(
      Mediator    = mediator,
      Parameter   = rn,
      Estimate    = vals[1],
      SE          = vals[2],
      L95         = vals[3],
      U95         = vals[4],
      Estimate_bt = trans$Estimate_bt,
      SE_approx   = trans$SE_approx,
      CI95_low    = trans$CI95_low,
      CI95_high   = trans$CI95_high
    )
  })
}

# Beta mediators
mediators_beta <- c("C3beta", "Pbeta", "Herbbeta", "Vivbeta")
bt_beta_PD     <- map_dfr(mediators_beta, ~ extract_backtrans(fx_PD, .x))

# PD rows (already on original scale)
PD_rows <- grep("^PD_", rownames(fx_PD), value = TRUE)

bt_PD <- map_dfr(PD_rows, function(rn) {
  vals <- as.numeric(fx_PD[rn, c("Estimate", "Est.Error", "Q2.5", "Q97.5")])
  tibble(
    Mediator    = "PD",
    Parameter   = rn,
    Estimate    = vals[1],
    SE          = vals[2],
    L95         = vals[3],
    U95         = vals[4],
    Estimate_bt = vals[1],   # original scale (gaussian, identity link)
    SE_approx   = vals[2],
    CI95_low    = vals[3],
    CI95_high   = vals[4]
  )
})

bt_all_PD <- bind_rows(bt_beta_PD, bt_PD) %>%
  arrange(Mediator, Parameter)

write.table(bt_all_PD, "Results/backtransformed_PD_2026.txt", sep = ";", row.names = FALSE)

# =====================================================================
# MS plots (standardized): Connectivity -> PD, Area -> C3, C3 -> PD, Viv -> PD
# =====================================================================

library(ggplot2)
library(gridExtra)

# 1) Connectivity (z_Conn) -> PD
ce_PD_Conn <- conditional_effects(
  fit_PD,
  effects = "z_Conn",
  resp    = "PD"
)

pred_PD_Conn <- as.data.frame(ce_PD_Conn[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = z_Conn)

p_PD_Conn <- ggplot() +
  geom_jitter(data = fit_PD$data,
              aes(x = z_Conn, y = PD),
              size = 2, alpha = 0.5, color = "#A1D99B",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_PD_Conn,
            aes(x = x, y = predicted),
            linewidth = 1.2, color = "#A1D99B") +
  geom_ribbon(data = pred_PD_Conn,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "#e0f3db", alpha = 0.35) +
  labs(x = "Connectivity", y = "Faith's PD") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# 2) Area (z_Area) -> C3 proportion (C3beta)
ce_C3_Area_PD <- conditional_effects(
  fit_PD,
  effects = "z_Area",
  resp    = "C3beta"
)

pred_C3_Area_PD <- as.data.frame(ce_C3_Area_PD[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = z_Area)

p_C3_Area_PD <- ggplot() +
  geom_jitter(data = fit_PD$data,
              aes(x = z_Area, y = C3_beta),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_C3_Area_PD,
            aes(x = x, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_C3_Area_PD,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  labs(x = "Area", y = "C3-grass") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# 3) C3 (zC3_logit) -> PD
ce_PD_C3 <- conditional_effects(
  fit_PD,
  effects = "zC3_logit",
  resp    = "PD"
)

pred_PD_C3 <- as.data.frame(ce_PD_C3[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = zC3_logit)

p_PD_C3 <- ggplot() +
  geom_jitter(data = fit_PD$data,
              aes(x = zC3_logit, y = PD),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_PD_C3,
            aes(x = x, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_PD_C3,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  labs(x = "C3-grass", y = "Faith's PD") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# 4) Viv (zViv_logit) -> PD
ce_PD_Viv <- conditional_effects(
  fit_PD,
  effects = "zViv_logit",
  resp    = "PD"
)

pred_PD_Viv <- as.data.frame(ce_PD_Viv[[1]]) %>%
  rename(predicted = estimate__,
         conf.low = lower__,
         conf.high = upper__,
         x = zViv_logit)

p_PD_Viv <- ggplot() +
  geom_jitter(data = fit_PD$data,
              aes(x = zViv_logit, y = PD),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_PD_Viv,
            aes(x = x, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_PD_Viv,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  labs(x = "Viv-grass", y = "Faith's PD") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# Panel MS en el orden pedido
grid.arrange(p_PD_Conn, p_C3_Area_PD, p_PD_C3, p_PD_Viv, ncol = 2)


# =====================================================================
# SM plots (real scale): Conn -> PD, Area (ha) -> PD, C3 (ha) -> PD, Viv (ha) -> PD
# =====================================================================

# 1) Back-transform helpers -------------------------------------------

mu_Area   <- attr(scale(SEM_data$Area),             "scaled:center")
sd_Area   <- attr(scale(SEM_data$Area),             "scaled:scale")
mu_Conn   <- attr(scale(SEM_data$Conn),             "scaled:center")
sd_Conn   <- attr(scale(SEM_data$Conn),             "scaled:scale")
mu_C3log  <- attr(scale(qlogis(SEM_data$C3_beta)),  "scaled:center")
sd_C3log  <- attr(scale(qlogis(SEM_data$C3_beta)),  "scaled:scale")
mu_Vivlog <- attr(scale(qlogis(SEM_data$Viv_beta)), "scaled:center")
sd_Vivlog <- attr(scale(qlogis(SEM_data$Viv_beta)), "scaled:scale")

z_to_Area   <- function(z) z * sd_Area  + mu_Area
z_to_Conn   <- function(z) z * sd_Conn  + mu_Conn
zC3_to_prop <- function(z) plogis(z * sd_C3log  + mu_C3log)
zViv_to_prop<- function(z) plogis(z * sd_Vivlog + mu_Vivlog)

# 2) Conditional effects (standardized scale) -------------------------

ce_PD_Conn <- conditional_effects(fit_PD, effects = "z_Conn",    resp = "PD")
ce_PD_Area <- conditional_effects(fit_PD, effects = "z_Area",    resp = "PD")
ce_PD_C3   <- conditional_effects(fit_PD, effects = "zC3_logit", resp = "PD")
ce_PD_Viv  <- conditional_effects(fit_PD, effects = "zViv_logit",resp = "PD")

pred_PD_Conn <- as.data.frame(ce_PD_Conn[[1]]) %>%
  rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = z_Conn)

pred_PD_Area <- as.data.frame(ce_PD_Area[[1]]) %>%
  rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = z_Area)

pred_PD_C3 <- as.data.frame(ce_PD_C3[[1]]) %>%
  rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = zC3_logit)

pred_PD_Viv <- as.data.frame(ce_PD_Viv[[1]]) %>%
  rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = zViv_logit)

# 3) Convert predictions + data to real scales ------------------------

pred_PD_Conn_real <- pred_PD_Conn %>%
  mutate(Conn_real = z_to_Conn(x))

pred_PD_Area_real <- pred_PD_Area %>%
  mutate(Area_real = z_to_Area(x))

pred_PD_C3_ha <- pred_PD_C3 %>%
  mutate(C3_prop = zC3_to_prop(x),
         C3_ha   = C3_prop * mean(SEM_data$Area, na.rm = TRUE))

pred_PD_Viv_ha <- pred_PD_Viv %>%
  mutate(Viv_prop = zViv_to_prop(x),
         Viv_ha   = Viv_prop * mean(SEM_data$Area, na.rm = TRUE))

SEM_data <- SEM_data %>%
  mutate(
    Conn_real = Conn,
    Area_real = Area,
    C3_ha     = C3_grass * Area,
    Viv_ha    = V_grass  * Area
  )

# 4) SM plot 1: Connectivity (real) -> PD -----------------------------

p_PD_Conn_real <- ggplot() +
  geom_jitter(data = SEM_data,
              aes(x = Conn_real, y = PD),
              size = 2, alpha = 0.5, color = "#A1D99B",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_PD_Conn_real,
            aes(x = Conn_real, y = predicted),
            linewidth = 1.2, color = "#A1D99B") +
  geom_ribbon(data = pred_PD_Conn_real,
              aes(x = Conn_real, ymin = conf.low, ymax = conf.high),
              fill = "#e0f3db", alpha = 0.35) +
  theme_classic(base_size = 14) +
  labs(x = "Connectivity (C500)", y = "Faith's PD")

# 5) SM plot 2: Area (ha) -> PD --------------------------------------

p_PD_Area_real <- ggplot() +
  geom_jitter(data = SEM_data,
              aes(x = Area_real, y = PD),
              size = 2, alpha = 0.5, color = "#A1D99B",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_PD_Area_real,
            aes(x = Area_real, y = predicted),
            linewidth = 1.2, color = "#A1D99B") +
  geom_ribbon(data = pred_PD_Area_real,
              aes(x = Area_real, ymin = conf.low, ymax = conf.high),
              fill = "#e0f3db", alpha = 0.35) +
  theme_classic(base_size = 14) +
  labs(x = "Area (ha)", y = "Faith's PD")

# 6) SM plot 3: C3 (ha) -> PD ----------------------------------------

p_PD_C3_ha <- ggplot() +
  geom_jitter(data = SEM_data,
              aes(x = C3_ha, y = PD),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_PD_C3_ha,
            aes(x = C3_ha, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_PD_C3_ha,
              aes(x = C3_ha, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  theme_classic(base_size = 14) +
  labs(x = "C3 area (ha)", y = "Faith's PD")

# 7) SM plot 4: Viv (ha) -> PD (0 hasta max Viv_ha) -------------------

max_Viv_ha_data <- max(SEM_data$Viv_ha, na.rm = TRUE)

min_row <- pred_PD_Viv_ha[which.min(pred_PD_Viv_ha$Viv_ha), ]
extra0  <- min_row; extra0$Viv_ha <- 0

pred_PD_Viv_ha_plot <- bind_rows(extra0, pred_PD_Viv_ha) %>%
  arrange(Viv_ha)

p_PD_Viv_ha <- ggplot() +
  geom_jitter(data = SEM_data,
              aes(x = Viv_ha, y = PD),
              size = 2, alpha = 0.5, color = "orange",
              height = 0.02, width = 0.2) +
  geom_line(data = pred_PD_Viv_ha_plot,
            aes(x = Viv_ha, y = predicted),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_PD_Viv_ha_plot,
              aes(x = Viv_ha, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff", alpha = 0.35) +
  coord_cartesian(xlim = c(0, max_Viv_ha_data), expand = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "Viv area (ha)", y = "Faith's PD")

# Panel (ajusta orden a gusto)
grid.arrange(p_PD_Conn_real, p_PD_Area_real, p_PD_C3_ha, p_PD_Viv_ha, ncol = 2)
