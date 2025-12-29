# ---------------------------------------------------------------------------------------------- #
# Bayesian SEM with brms: Beta mediators + hurdle lognormal abundance for RAO
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)

# -----------------------------
# Load and prepare data
# -----------------------------
space <- read.delim("DATA/space_data.txt", sep = "\t")
fd.df <- read.delim("DATA/fd_df.txt", sep = ";")

space <- space %>%
  mutate(
    C3.proportion   = C3ha / TotalA, # ornamental veg., hugh hyric resources, highly managed (sega y riego)
    P.proportion    = PUha / TotalA, # native veg. with spontaneous herbs and flowers, low management(one sega por año??)
    Herb.proportion = pmin(HERbha / TotalA, 1), # native herbs and flowers managed (water) but no siega (only once?? - IN Biod types)
    Viv.proportion  = VIVha / TotalA # ??
  )

SEM_data <- space %>%
  select(SITE_ID, Area = TotalA, Conn = C500, Type,
         C3_grass = C3.proportion, P_grass = P.proportion,
         Herb_grass = Herb.proportion, V_grass = Viv.proportion) %>%
  mutate(
    z_Area = as.numeric(scale(Area)),
    z_Conn = as.numeric(scale(Conn)),
    # Smithson–Verkuilen para proporciones de vegetación
    n = n(),
    C3_beta   = (pmin(pmax(C3_grass,0),1)*(n-1)+0.5)/n,
    P_beta    = (pmin(pmax(P_grass ,0),1)*(n-1)+0.5)/n,
    Herb_beta = (pmin(pmax(Herb_grass,0),1)*(n-1)+0.5)/n,
    Viv_beta  = (pmin(pmax(V_grass ,0),1)*(n-1)+0.5)/n
  )

SEM_data <- SEM_data %>% right_join(fd.df, by = "SITE_ID") 

# -----------------------------
# Fit SEM 
# -----------------------------

# 1.Model formulas
bf_C3   <- bf(C3_beta   ~ z_Area + Type, family = Beta(link = "logit"))
bf_P    <- bf(P_beta    ~ z_Area + Type, family = Beta(link = "logit"))
bf_Herb <- bf(Herb_beta ~ z_Area + Type, family = Beta(link = "logit"))
bf_Viv  <- bf(Viv_beta  ~ z_Area + Type, family = Beta(link = "logit"))

bf_Rao  <- bf(RaoQ      ~ z_Area + z_Conn + Type + C3_beta + P_beta + Herb_beta + Viv_beta + (1 | YEAR),
              family = gaussian())

# 2.Fit Model
fit_comm_min <- brm(
  bf_C3 + bf_P + bf_Herb + bf_Viv + bf_Rao + set_rescor(FALSE),
  data = SEM_data,
  chains = 4, iter = 120000, warmup = 20000, cores = parallel::detectCores(),
    control = list(adapt_delta = 0.999, max_treedepth = 15),# updated to improve convergence
    seed = 20251004
  )
saveRDS(fit_comm_min, file = "Results/RaO_brms_SEM.rds")

# 3. Load models saved 
fit_rao <- readRDS("Results/RaO_brms_SEM.rds")
summary(fit_rao)

# check in summary for R hat for each parameter SHOULD be < 1.05
range(neff_ratio(fit_rao), na.rm = TRUE)  # Effective sample size ratio, should be > 0.1 (better if > 0.2).

# 3. Posterior predictive checks 
pp_check(fit_rao, resp = "RaoQ") 

# 4. Marginal effects for plotting 
ce_rao <- conditional_effects(fit_rao, resp = "RaoQ")
plot(ce_rao, points = TRUE)

