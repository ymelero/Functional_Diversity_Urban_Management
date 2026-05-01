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

mint  <- read.delim("DATA/MIntensity.txt", sep = "\t")
names(mint)[4] <- "Intensity"   # columna con low/medium/high

mint <- mint %>%
  filter(Vegetation == "Herb") %>%
  select(,-c(2,3))   # ajusta si los nombres difieren

space <- space %>%
  left_join(mint, by = "Type")



# 1. Space composition versus Management
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
  labs(x = "Vegetatiob group", y = "Proportion"
  ) +
  theme(legend.position = "none",
        plot.title = element_text(size = 13, face = "bold"))

# Boxplot by Intensity
# remove C3 for better visualisation
space_long_red <- space_long %>% filter(Group !="C3")
cols <- RColorBrewer::brewer.pal(4, "Set2")[-1]
p2<-ggplot(space_long_red, aes(x = Intensity, y = Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.6, width = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = cols) +  
  theme_classic(base_size = 12) + ylim(0, 0.15) +
  labs(x = "M intensity", y = "Proportion") +
  theme(legend.title = element_blank(),
        legend.position = "top")


# 2. Area and Connectivity versus M. Intensity
p_area <- ggplot(space, aes(x = Intensity, y = TotalA, fill = Intensity)) +
  geom_boxplot(alpha = 0.6, width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  labs(x = "", y = "Area") +
  theme(legend.position = "top") + 
  theme(axis.ticks.x = element_blank()) +
  annotate("text", x = -Inf, y = Inf, label = "(a)", hjust = -0.5, vjust = 1.5, size= 3.5)

kruskal.test(TotalA ~ Intensity, data = space)
pairwise.wilcox.test(space$TotalA, space$Intensity, p.adjust.method = "holm")

p_con <- ggplot(space, aes(x = Intensity, y = C500, fill = Intensity)) +
  geom_boxplot(alpha = 0.6, width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  labs(x = "", y = "Connectivity") +
  theme(legend.position = "none") + 
  theme(axis.ticks.x = element_blank()) +
  annotate("text", x = -Inf, y = Inf, label = "(b)", hjust = -0.5, vjust = 1.5, size= 3.5)

kruskal.test(C500 ~ Intensity, data = space)
pairwise.wilcox.test(space$C500, space$Intensity, p.adjust.method = "holm")

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


