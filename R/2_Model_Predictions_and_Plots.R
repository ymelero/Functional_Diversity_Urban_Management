# ---------------------------------------------------------------------------------------------- #
# Back-transform estimates (Estimate, SE, CI) from logit to proportion scale
# for the Beta submodels (C3, P, Herb, Viv)
# AND from log / logit scale for abundance and hurdle components
# Predict (marginal predictions) from SEM models 
# Predict and plot effects
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(tidyverse)
library(purrr)
library(brms)
detach("package:ggeffects")
detach("package:gridExtra")

fit_c2 <- readRDS("Results/Abundance_brms_SEM_Cluster2.rds") #load the corresponding one
fx_all <- as.data.frame(fixef(fit_c2)) %>% mutate(across(everything(), as.numeric))

# 1. Site variables
inv_logit <- function(x) plogis(as.numeric(x))
inv_exp   <- function(x) exp(as.numeric(x))

# Back-transform one row for Beta models (logit → proportion)
backtransform_row <- function(estimate, se, l95, u95) {
  tibble(
    Estimate_bt = inv_logit(estimate),
    SE_approx   = (inv_logit(estimate + se) - inv_logit(estimate - se)) / 2,
    CI95_low    = inv_logit(l95),
    CI95_high   = inv_logit(u95)
  )
}

# Apply to one mediator (C3, P, Herb, Viv)
extract_backtrans <- function(fx, mediator) {
  rows <- grep(paste0("^", mediator, "_"), rownames(fx), value = TRUE)
  
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
bt_beta <- map_dfr(mediators_beta, ~ extract_backtrans(fx_all, .x))


# 2. Back-transform abundance (log) and hurdle (logit) parts
backtransform_abundance_row <- function(estimate, se, l95, u95, type) {
  if (type == "log") {
    tibble(
      Estimate_bt = inv_exp(estimate),
      SE_approx   = (inv_exp(estimate + se) - inv_exp(estimate - se)) / 2,
      CI95_low    = inv_exp(l95),
      CI95_high   = inv_exp(u95)
    )
  } else if (type == "logit") {
    tibble(
      Estimate_bt = inv_logit(estimate),
      SE_approx   = (inv_logit(estimate + se) - inv_logit(estimate - se)) / 2,
      CI95_low    = inv_logit(l95),
      CI95_high   = inv_logit(u95)
    )
  }
}

extract_backtrans_abundance <- function(fx, mediator, type) {
  rows <- grep(paste0("^", mediator, "_"), rownames(fx), value = TRUE)
  
  map_dfr(rows, function(rn) {
    vals <- as.numeric(fx[rn, c("Estimate", "Est.Error", "Q2.5", "Q97.5")])
    trans <- backtransform_abundance_row(vals[1], vals[2], vals[3], vals[4], type)
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

# Run for both abundance parts
bt_abund <- extract_backtrans_abundance(fx_all, "abundance", "log")
bt_hu    <- extract_backtrans_abundance(fx_all, "hu_abundance", "logit")

# Merge all results into one unified table
bt_all_full <- bind_rows(bt_beta, bt_abund, bt_hu) %>%
  arrange(Mediator, Parameter)
#write.table(bt_all_full, "Results/backtransdormed_C2.txt", sep=";", row.names = F)

# 3. Predict & Plot
library(ggeffects)
library(gridExtra)

# Function for main responses (Beta)
get_preds <- function(mediator) {
  ce <- conditional_effects(fit_c2, effects = "Type", resp = mediator)
  
  as.data.frame(ce[[1]]) %>%
    rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = Type) %>%
    mutate(Mediator = mediator)
}

# Function for hurdle component (distributional parameter 'hu')
get_preds_hurdle <- function() {
  ce <- conditional_effects(fit_c2, effects = "Type", resp = "abundance", dpar = "hu")
  
  as.data.frame(ce[[1]]) %>%
    rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = Type) %>%
    mutate(Mediator = "hu_abundance")
}

# Mediators present in the model
mediators_main <- c("C3beta", "Pbeta", "Herbbeta", "Vivbeta", "abundance")

# Run predictions for main responses
preds_main <- purrr::map_dfr(mediators_main, get_preds)

# Run predictions for hurdle (probability of zero)
preds_hu <- get_preds_hurdle()

# Combine all predictions
preds_all <- bind_rows(preds_main, preds_hu)

# Format for plotting
preds_all <- preds_all %>%
  mutate(
    Type = factor(x, levels = c("BIO", "PATRIMONIAL", "USOS"),
                  labels = c("Bio.", "Hist.", "Soc.")),
    Mediator = factor(Mediator,
                      levels = c("C3beta", "Pbeta", "Herbbeta", "Vivbeta",
                                 "abundance", "hu_abundance"),
                      labels = c("C3", "P", "Herb", "Viv", "Abundance", "Hurdle (P0)"))
  )



# Plot per variable
preds_props <- preds_all %>%
  filter(Mediator %in% c("C3", "P", "Herb", "Viv"))

p1 <- ggplot(preds_props, aes(x = Type, y = predicted, color = Mediator)) +
  geom_point(position = position_dodge(width = 0.8), size = 1.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Space Type", y = "Vegetation") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12)
  )
p1

# predict and plot area effect on vegetation:
get_preds_area <- function(mediator) {
  ce <- conditional_effects(fit_c2, effects = "z_Area", resp = mediator)
  
  as.data.frame(ce[[1]]) %>%
    rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = z_Area) %>%
    mutate(Mediator = mediator)
}
mediators_area <- c("C3beta", "Pbeta", "Herbbeta", "Vivbeta")
preds_area <- purrr::map_dfr(mediators_area, get_preds_area)

# Format for plotting
preds_area <- preds_area %>%
  mutate(Mediator = factor(Mediator,
                           levels = c("C3beta", "Pbeta", "Herbbeta", "Vivbeta"),
                           labels = c("C3", "P", "Herb", "Viv")))

p_area <- ggplot(preds_area, aes(x = x, y = predicted, color = Mediator)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Mediator), alpha = 0.15, color = NA) +
  labs(x = "Area", y = "Vegetation") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

p_area

# Plot PRESENCE and ABUNDANCE in relation to direct and indirect effects
# Direct effect of Space type on P  and N:
preds_hu_plot <- preds_all %>%
  filter(Mediator == "Hurdle (P0)")

p.hu.type <- ggplot(preds_hu_plot, aes(x = Type, y = 1-predicted, color = Type)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = 1- conf.low, ymax = 1-conf.high), width = 0.25) +
  labs(x = "Space Type", y = "Presence") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12)
  )
p.hu.type

preds_abund <- preds_all %>%
  filter(Mediator == "Abundance")

p.N.type <- ggplot(preds_abund, aes(x = Type, y = predicted, color = Type)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.25) +
  labs(x = "Space Type", y = "Abundance") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12)
  ) 
p.N.type


# Variable continuous
make_seq <- function(x) seq(min(x, na.rm = TRUE),
                            max(x, na.rm = TRUE),
                            length.out = 200)

seq_list <- list(
  Conn = make_seq(fit_c2$data$z_Conn),
  Area = make_seq(fit_c2$data$z_Area),
  C3   = make_seq(fit_c2$data$zC3_logit),
  P    = make_seq(fit_c2$data$zP_logit),
  Herb = make_seq(fit_c2$data$zHerb_logit),
  Viv  = make_seq(fit_c2$data$zViv_logit)
)

valid_YEAR <- if ("YEAR" %in% names(fit_c2$data)) {
  names(sort(table(fit_c2$data$YEAR), decreasing = TRUE))[1]
} else NULL

valid_SITE <- if ("SITE_ID" %in% names(fit_c2$data)) {
  names(sort(table(fit_c2$data$SITE_ID), decreasing = TRUE))[1]
} else NULL


predict_manual <- function(effect, seqx) {
  
  # build newdata
  nd <- tibble(!!effect := seqx)
  
  # numeric predictors fixed to mean
  for (v in names(fit_c2$data)) {
    if (!(v %in% names(nd)) && is.numeric(fit_c2$data[[v]])) {
      nd[[v]] <- mean(fit_c2$data[[v]], na.rm = TRUE)
    }
  }
  
  # fix factor Type
  if ("Type" %in% names(fit_c2$data)) {
    nd$Type <- factor("BIO", levels = levels(fit_c2$data$Type))
  }
  
  # random effects
  if (!is.null(valid_YEAR)) {
    nd$YEAR <- factor(valid_YEAR, levels = levels(fit_c2$data$YEAR))
  }
  
  if (!is.null(valid_SITE)) {
    nd$SITE_ID <- factor(valid_SITE, levels = levels(fit_c2$data$SITE_ID))
  }
  
  # abundance prediction
  ab <- fitted(
    fit_c2,
    newdata = nd,
    resp = "abundance",
    allow_new_levels = TRUE,
    summary = TRUE
  )
  
  # presence (1 - hurdle)
  hu <- fitted(
    fit_c2,
    newdata = nd,
    resp = "abundance",
    dpar = "hu",
    allow_new_levels = TRUE,
    summary = TRUE
  )
  
  tibble(
    x = seqx,
    predicted_abund = ab[, "Estimate"],
    abund_low       = ab[, "Q2.5"],
    abund_high      = ab[, "Q97.5"],
    predicted_pres  = 1 - hu[, "Estimate"],
    pres_low        = 1 - hu[, "Q97.5"],
    pres_high       = 1 - hu[, "Q2.5"]
  )
}

pred_conn <- predict_manual("z_Conn",      seq_list$Conn)
pred_area <- predict_manual("z_Area",      seq_list$Area)
pred_C3   <- predict_manual("zC3_logit",   seq_list$C3)
pred_P    <- predict_manual("zP_logit",    seq_list$P)
pred_Herb <- predict_manual("zHerb_logit", seq_list$Herb)
pred_Viv  <- predict_manual("zViv_logit",  seq_list$Viv)

#############################################
# Presence vs Connectivity
#############################################

pres.conn <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = z_Conn, y = ifelse(abundance > 0, 1, 0)),
              size = 2, alpha = 0.3,
              color = "#A1D99B",
              height = 0.05, width = 0.1) +
  geom_line(data = pred_conn,
            aes(x = x, y = predicted_pres),
            linewidth = 1.2, color = "#A1D99B") +
  geom_ribbon(data = pred_conn,
              aes(x = x, ymin = pres_low, ymax = pres_high),
              fill = "#e0f3db", alpha = 0.15) +
  labs(x = "Connectivity", y = "Presence") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


#############################################
# Abundance vs Connectivity
#############################################

abdn.conn <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = z_Conn, y = abundance),
              size = 2, alpha = 1,
              width = 0.05, height = 0.05,
              color = "#A1D99B",
              shape = 23) +
  geom_line(data = pred_conn,
            aes(x = x, y = predicted_abund),
            linewidth = 1.3, color = "#A1D99B") +
  geom_ribbon(data = pred_conn,
              aes(x = x, ymin = abund_low, ymax = abund_high),
              fill = "#e0f3db", alpha = 0.15) +
  labs(x = "Connectivity", y = "Abundance") +
  ylim(0, 50) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


#############################################
# Presence vs Area
#############################################

pres.area <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = z_Area, y = ifelse(abundance > 0, 1, 0)),
              size = 2, alpha = 0.3,
              color = "#A1D99B",
              height = 0.05, width = 0.09) +
  geom_line(data = pred_area,
            aes(x = x, y = predicted_pres),
            linewidth = 1.2, color = "#A1D99B") +
  geom_ribbon(data = pred_area,
              aes(x = x, ymin = pres_low, ymax = pres_high),
              fill = "#e0f3db", alpha = 0.15) +
  labs(x = "Area", y = "Presence") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


#############################################
# Abundance vs Area
#############################################

abdn.area <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = z_Area, y = abundance),
              size = 2, alpha = 0.35,
              width = 0.05, height = 0.05,
              color = "#FF8CC6") +
  geom_line(data = pred_area,
            aes(x = x, y = predicted_abund),
            linewidth = 1.3, color = "#FF8CC6") +
  geom_ribbon(data = pred_area,
              aes(x = x, ymin = abund_low, ymax = abund_high),
              fill = "#FFD6EA", alpha = 0.35) +
  labs(x = "Area", y = "Abundance") +
  ylim(0, 100) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


#############################################
# Presence vs C3
#############################################

pres.C3 <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = zC3_logit, y = ifelse(abundance > 0, 1, 0)),
              size = 2, alpha = 0.25,
              color = "orange",
              height = 0.05, width = 0.3) +
  geom_line(data = pred_C3,
            aes(x = x, y = predicted_pres),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_C3,
              aes(x = x, ymin = pres_low, ymax = pres_high),
              fill = "peachpuff", alpha = 0.15) +
  labs(x = "C3-grass", y = "Presence") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


#############################################
# Abundance vs C3
#############################################

abdn.C3 <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = zC3_logit, y = abundance),
              size = 2, alpha = 0.7,
              width = 0.05, height = 0.15,
              color = "orange",
              shape = 23) +
  geom_line(data = pred_C3,
            aes(x = x, y = predicted_abund),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_C3,
              aes(x = x, ymin = abund_low, ymax = abund_high),
              fill = "peachpuff", alpha = 0.25) +
  labs(x = "C3-grass", y = "Abundance") +
  ylim(0, 50) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


#############################################
# Presence vs P
#############################################

pres.P <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = zP_logit, y = ifelse(abundance > 0, 1, 0)),
              size = 2, alpha = 0.09,
              color = "#A1D99B",
              height = 0.05, width = 0.09) +
  geom_line(data = pred_P,
            aes(x = x, y = predicted_pres),
            linewidth = 1.2, color = "#A1D99B") +
  geom_ribbon(data = pred_P,
              aes(x = x, ymin = pres_low, ymax = pres_high),
              fill = "#e0f3db", alpha = 0.15) +
  labs(x = "P-grass", y = "Presence") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


#############################################
# Abundance vs P
#############################################

abdn.P <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = zP_logit, y = abundance),
              size = 2, alpha = 0.35,
              width = 0.05, height = 0.05,
              color = "#FF8CC6") +
  geom_line(data = pred_P,
            aes(x = x, y = predicted_abund),
            linewidth = 1.2, color = "#FF8CC6") +
  geom_ribbon(data = pred_P,
              aes(x = x, ymin = abund_low, ymax = abund_high),
              fill = "#FFD6EA", alpha = 0.35) +
  labs(x = "P-grass", y = "Abundance") +
  ylim(0, 100) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


#############################################
# Presence vs Herb
#############################################

pres.Herb <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = zHerb_logit, y = ifelse(abundance > 0, 1, 0)),
              size = 2, alpha = 0.3,
              color = "#A1D99B",
              height = 0.05, width = 0.09) +
  geom_line(data = pred_Herb,
            aes(x = x, y = predicted_pres),
            linewidth = 1.2, color = "#A1D99B") +
  geom_ribbon(data = pred_Herb,
              aes(x = x, ymin = pres_low, ymax = pres_high),
              fill = "#e0f3db", alpha = 0.15) +
  labs(x = "Herb-grass", y = "Presence") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


#############################################
# Abundance vs Herb
#############################################

abdn.Herb <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = zHerb_logit, y = abundance),
              size = 2, alpha = 0.7,
              width = 0.05, height = 0.15,
              color = "orange",
              shape = 23) +
  geom_line(data = pred_Herb,
            aes(x = x, y = predicted_abund),
            linewidth = 1.2, color = "peachpuff") +
  geom_ribbon(data = pred_Herb,
              aes(x = x, ymin = abund_low, ymax = abund_high),
              fill = "peachpuff", alpha = 0.25) +
  labs(x = "Herb-grass", y = "Abundance") +
  ylim(0, 50) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


#############################################
# Presence vs Viv
#############################################

pres.Viv <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = zViv_logit, y = ifelse(abundance > 0, 1, 0)),
              size = 2, alpha = 0.3,
              color = "#A1D99B",
              height = 0.05, width = 0.09) +
  geom_line(data = pred_Viv,
            aes(x = x, y = predicted_pres),
            linewidth = 1.2, color = "#A1D99B") +
  geom_ribbon(data = pred_Viv,
              aes(x = x, ymin = pres_low, ymax = pres_high),
              fill = "#e0f3db", alpha = 0.15) +
  labs(x = "Viv-grass", y = "Presence") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


#############################################
# Abundance vs Viv
#############################################

abdn.Viv <- ggplot() +
  geom_jitter(data = fit_c2$data,
              aes(x = zViv_logit, y = abundance),
              size = 2, alpha = 0.35,
              width = 0.05, height = 0.05,
              color = "#FF8CC6") +
  geom_line(data = pred_Viv,
            aes(x = x, y = predicted_abund),
            linewidth = 1.2, color = "#FF8CC6") +
  geom_ribbon(data = pred_Viv,
              aes(x = x, ymin = abund_low, ymax = abund_high),
              fill = "#FFD6EA", alpha = 0.35) +
  labs(x = "Viv-grass", y = "Abundance") +
  ylim(0, 100) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


# names are: 
#p.hu.type 
#p.N.type
#pres.conn
#abdn.conn
#pres.conn
#abdn.conn
#pres.area
#abdn.area
#pres.C3
#abdn.C3
#pres.P
#abdn.P
#pres.Herb
#abdn.Herb
#pres.Viv
#abdn.Viv

library(patchwork)
add_tag <- function(plot, tag) {
  plot +
    annotation_custom(
      grob = textGrob(tag,
                      x = 0.02, y = 0.98,       # esquina superior izquierda
                      hjust = 0, vjust = 1,
                      gp = gpar(fontsize = 14)),  # sin bold
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
}


# FG2
p_fg2 <-
  (
    add_tag(pres.conn, "(a)") |
      add_tag(p.hu.type, "(b)")
  ) /
  (
    add_tag(pres.C3, "(c)")   |
      add_tag(pres.P, "(d)")
  ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 14)
  )

# no legends    
p_fg2 <- p_fg2 + plot_layout(guides = "collect") &
  theme(legend.position = "none")
p_fg2

# FG3
p_fg3 <- (
    add_tag(pres.conn, "(a)") |
      add_tag(abdn.conn, "(b)") |
    add_tag(pres.area, "(c)") |
      add_tag(p.hu.type, "(d)")
  ) /
  (
    add_tag(pres.Viv, "(e)")   |
      add_tag(pres.C3, "(f)") |
      add_tag(abdn.C3, "(g)") 
  ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 14)
  )

# no legends    
p_fg3 <- p_fg3 + plot_layout(guides = "collect") &
  theme(legend.position = "none")

# FG4
p_fg4 <- (
  add_tag(pres.area , "(a)") |
    add_tag(pres.C3 , "(b)") |
    add_tag(abdn.Herb, "(c)")) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 14)
  )

# no legends    
p_fg4 <- p_fg4 + plot_layout(guides = "collect") &
  theme(legend.position = "none")

