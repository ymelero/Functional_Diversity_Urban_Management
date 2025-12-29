# ---------------------------------------------------------------------------------------------- #
# Raw / Proportional Data Exploratory Plots
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #


library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lme4)
library(ggeffects)

space    <- read.delim("DATA/space_data.txt", sep = "\t")
space <- space %>%
  mutate(
    C3.proportion   = C3ha / TotalA, 
    P.proportion    = PUha / TotalA, 
    Herb.proportion = pmin(HERbha / TotalA, 1), 
    Viv.proportion  = VIVha / TotalA 
  )

# 1. Space composition versus Type
# Reshape to long format for easy plotting
space_long <- space %>%
  pivot_longer(cols = c(C3.proportion, P.proportion, Herb.proportion, Viv.proportion),
               names_to = "Group", values_to = "Proportion")

# Clean group names
space_long$Group <- factor(space_long$Group,
                           levels = c("C3.proportion", "P.proportion", "Herb.proportion", "Viv.proportion"),
                           labels = c("C3", "P", "Herb", "Viv"))

# Boxplot (overall)
ggplot(space_long, aes(x = Group, y = Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.6, width = 0.7, outlier.shape = 21) +
  scale_fill_brewer(palette = "Set2") + ylim(0, 0.5) +
  theme_classic(base_size = 12) +
  labs(x = "Functional group", y = "Proportion"
       ) +
  theme(legend.position = "none",
        plot.title = element_text(size = 13, face = "bold"))

# Boxplot by Type
p1<-ggplot(space_long, aes(x = Type, y = Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.6, width = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) + #ylim(0, 0.5) +
  labs(x = "", y = "Proportion") +
  theme(legend.title = element_blank(),
        legend.position = "top") + theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())+
  annotate("text", x = -Inf, y = Inf, label = "(a)", hjust = -0.5, vjust = 1.5, size= 3.5)

# remove C3 for better visulisation
space_long_red <- space_long %>% filter(Group !="C3")
cols <- RColorBrewer::brewer.pal(4, "Set2")[-1]
p2<-ggplot(space_long_red, aes(x = Type, y = Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.6, width = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = cols) +  
  theme_classic(base_size = 12) + ylim(0, 0.15) +
  labs(x = "Space type", y = "Proportion") +
  theme(legend.position = "none") +
  annotate("text", x = -Inf, y = Inf, label = "(b)", hjust = -0.5, vjust = 1.5, size= 3.5)

grid.arrange(p1, p2, ncol = 1, nrow = 2)

# 2. Area & Connectivity versus Type
p_area <- ggplot(space, aes(x = Type, y = TotalA, fill = Type)) +
  geom_boxplot(alpha = 0.6, width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  labs(x = "", y = "Area") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  annotate("text", x = -Inf, y = Inf, label = "(a)", hjust = -0.5, vjust = 1.5, size= 3.5)

p_con <- ggplot(space, aes(x = Type, y = C500, fill = Type)) +
  geom_boxplot(alpha = 0.6, width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  labs(x = "Space type", y = "Connectivity") +
  theme(legend.position = "none")+
  annotate("text", x = -Inf, y = Inf, label = "(b)", hjust = -0.5, vjust = 1.5, size= 3.5)

grid.arrange(p_area, p_con, ncol = 1, nrow = 2)

# 3. Plot variables correlations
library(corrplot) 
library(psych)
p.var <- space[, c(4:9,24:27)]
p.var <- p.var %>% rename("Open Area" = OpenA,
                          "Closed Area" = ClosedA,
                          "Total Area" = TotalA,
                          "Conn. 200m" = C200,
                          "Conn. 500m" = C500,
                          "Conn. 1km" = C1k,
                          "C3 Proportion" = C3.proportion,
                          "P Proportion" = P.proportion,
                          "Herb. Proportion" = Herb.proportion,
                          "Viv. Proportion" = Viv.proportion
                          )
M <- corr.test(p.var)
corrplot(M$r, method = "number", type = "upper", sig.level = c(0.001, 0.01, 0.05), 
         insig = "label_sig", pch.cex = 0.9, pch.col = "grey20", order = "original")
corrplot(M$r, p.mat = M$p, type = "upper", sig.level = c(0.001, 0.01, 0.05), 
         insig = "label_sig", pch.cex = 0.9, pch.col = "grey", order = "original", tl.col = "black")

# 4. Plot raw presence and abundance per cluster
# read ubms and cluster data from 1_Model_Bayesian....
# descriptive data:
# Presence
ubmsdata <- ubmsdata %>%
  mutate(presence = ifelse(SINDEX > 0, 1, 0))

presence_summary <- ubmsdata %>%
  group_by(Cluster) %>%
  summarise(
    n = n(),
    n_presence = sum(presence, na.rm = TRUE),
    prop_presence = n_presence / n,
    lower_CI = map2_dbl(n_presence, n, ~ prop.test(.x, .y)$conf.int[1]),
    upper_CI = map2_dbl(n_presence, n, ~ prop.test(.x, .y)$conf.int[2])
  )

glm_presence <- glmer(presence ~ Cluster + (1 | YEAR) + (1 | SITE_ID), data = ubmsdata, family = binomial)
summary(glm_presence)
pred_cluster <- ggeffect(glm_presence, terms = "Cluster")

p_presence <- ggplot() +
  geom_jitter(data = ubmsdata, 
              aes(x = Cluster, y = presence, color = Cluster),
              width = 0.5, height = 0.2,alpha = 0.1, size = 1.5) +
  geom_point(data = pred_cluster, 
             aes(x = x, y = predicted, color = x), 
             size = 3) +
  geom_errorbar(data = pred_cluster, 
                aes(x = x, ymin = conf.low, ymax = conf.high, color = x), 
                width = 0.15, linewidth = 0.6) +
  annotate("text", x = -Inf, y = Inf, label = "(c)", hjust = -0.5, vjust = 1.5, size= 3.5)+
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  labs(x = "", y = "Presence") +
  theme(legend.position = "none") + coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(labels = c("2" = "C2", "3" = "C3", "4" = "C4"))

# Abundance
abund_data <- ubmsdata %>%
  filter(SINDEX > 0)
abund_data %>%
  group_by(Cluster) %>%
  summarise(
    n = n(),
    mean = round(mean(SINDEX, na.rm = TRUE), 2),
    sd = round(sd(SINDEX, na.rm = TRUE), 2),
    min = round(min(SINDEX, na.rm = TRUE), 2),
    max = round(max(SINDEX, na.rm = TRUE), 2)
  )

glmer_abund <- glmer(SINDEX ~ Cluster + (1 | YEAR) + (1 | SITE_ID),
                         data = abund_data,
                         family = Gamma(link = "log"))
summary(glmer_abund)
pred_abund <- ggeffect(glmer_abund, terms = "Cluster")

p_abundance <- ggplot() +
  geom_jitter(data = abund_data,
              aes(x = Cluster, y = SINDEX, color = Cluster),
              width = 0.5, height = 0.2, alpha = 0.1, size = 1.5) +
  geom_point(data = pred_abund,
             aes(x = x, y = predicted, color = x),
             size = 3) +
  geom_errorbar(data = pred_abund,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = x),
                width = 0.15, linewidth = 0.6) +
  scale_color_brewer(
    palette = "Set2",
    name = "Functional cluster",
    labels = c("C2", "C3", "C4")
  ) +
  theme_classic(base_size = 12) +
  labs(x = "", y = "Abundance") +
  theme(
    legend.position = "none", axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

p_abdn <- p_abundance + annotate("text", x = -Inf, y = Inf, label = "(a)", hjust = -0.5, vjust = 1.5, size= 3.5)+
  theme(
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10))
p_abundance_zoom <- p_abundance + coord_cartesian(ylim = c(0, 100)) + 
  annotate("text", x = -Inf, y = Inf, label = "(b)", hjust = -0.5, vjust = 1.5, size= 3.5)

grid.arrange(p_abdn, p_abundance_zoom, p_presence, ncol = 1, nrow = 3)