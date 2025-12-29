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

fit_rao <- readRDS("Results/RaO_brms_SEM.rds")
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
write.table(bt_all, "Results/backtransformed_RAO.txt", sep = ";", row.names = FALSE)

# 3. Predict & Plot
library(ggeffects)
library(gridExtra)

# Function for main responses & RAO
get_preds <- function(mediator) {
  ce <- conditional_effects(fit_rao, effects = "Type", resp = mediator)
  
  as.data.frame(ce[[1]]) %>%
    rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = Type) %>%
    mutate(Mediator = mediator)
}

mediators_main <- c("C3beta", "Pbeta", "Herbbeta", "Vivbeta")
preds_main <- purrr::map_dfr(mediators_main, get_preds)

preds_rao <- get_preds("RaoQ")


# Format for plotting
preds_all <- bind_rows(preds_main, preds_rao) %>%
  mutate(
    Type = factor(x, levels = c("BIO","PATRIMONIAL","USOS"),
                  labels = c("Bio.","Hist.","Soc.")),
    Mediator = factor(Mediator,
                      levels = c("C3beta","Pbeta","Herbbeta","Vivbeta","RaoQ"),
                      labels = c("C3","P","Herb","Viv","RaoQ"))
  )


# Plot (optional)
ggplot(preds_all, aes(x = Type, y = predicted, color = Mediator)) +
  geom_point(position = position_dodge(width = 0.8), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Space Type", y = "Proportion") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12)
  )


# separado
preds_props <- preds_all %>%
  filter(Mediator %in% c("C3", "P", "Herb", "Viv"))

p1 <- ggplot(preds_props, aes(x = Type, y = predicted, color = Mediator)) +
  geom_point(position = position_dodge(width = 0.8), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Space Type", y = "Vegetation") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12)
  )
p1

p_rao <- preds_all %>%
  filter(Mediator == "RaoQ")

p2<- ggplot(p_rao, aes(x = Type, y = predicted, color = Type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.25, position = position_dodge(width = 0.5)) +
  labs(x = "Space Type", y = "RaoQ") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12)
  )
p2

# for RaoQ the only actually significant was the C3 grass (negatively), which at the same time is - affected by Bio type
valid_YEAR <- if ("YEAR" %in% names(fit_rao$data)) {
  names(sort(table(fit_rao$data$YEAR), decreasing = TRUE))[1]
} else NULL

valid_SITE <- if ("SITE_ID" %in% names(fit_rao$data)) {
  names(sort(table(fit_rao$data$SITE_ID), decreasing = TRUE))[1]
} else NULL

seq_C3 <- seq(
  min(fit_rao$data$zC3_logit, na.rm = TRUE),
  max(fit_rao$data$zC3_logit, na.rm = TRUE),
  length.out = 200
)

predict_manual_RAO <- function(effect, seqx) {
  nd <- tibble(!!effect := seqx)
  for (v in names(fit_rao$data)) {
    if (!(v %in% names(nd)) && v != effect && is.numeric(fit_rao$data[[v]])) {
      nd[[v]] <- mean(fit_rao$data[[v]], na.rm = TRUE)
    }
  }
  if ("Type" %in% names(fit_rao$data)) {
    nd$Type <- factor("BIO", levels = levels(fit_rao$data$Type))
  }
  if (!is.null(valid_YEAR)) {
    nd$YEAR <- factor(valid_YEAR, levels = levels(fit_rao$data$YEAR))
  }
  if (!is.null(valid_SITE)) {
    nd$SITE_ID <- factor(valid_SITE, levels = levels(fit_rao$data$SITE_ID))
  }
  fr <- fitted(
    fit_rao,
    newdata = nd,
    resp = "RaoQ",
    allow_new_levels = TRUE,
    summary = TRUE
  )
  tibble(
    x = seqx,
    predicted = fr[, "Estimate"],
    conf.low  = fr[, "Q2.5"],
    conf.high = fr[, "Q97.5"]
  )
}

pred_Rao_C3 <- predict_manual_RAO("zC3_logit", seq_C3)

p_rao_C3 <- ggplot() +
  geom_jitter(data = fit_rao$data,
              aes(x = zC3_logit, y = RaoQ),
              size = 2,
              alpha = 0.25,
              color = "orange",
              height = 0.05,
              width = 0.3) +
  geom_line(data = pred_Rao_C3,
            aes(x = x, y = predicted),
            linewidth = 1.2,
            color = "peachpuff") +
  geom_ribbon(data = pred_Rao_C3,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "peachpuff",
              alpha = 0.15) +
  labs(x = "C3-grass", y = "RaoQ") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 11),
    axis.title   = element_text(size = 12)
  )

p_rao_C3

p_all<-grid.arrange(p1, p_rao_C3, ncol =3, nrow = 1)
