# ---------------------------------------------------------------------------------------------- #
# Rao’s Quadratic Entropy (RaoQ), Shannon y Faith's PD  per SITE × YEAR
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #


library(tidyverse)
library(FD)
library(vegan)
library(ape)
library(picante)

# 0. Load and prepare data
# SINDEX_ZEROS contains counts and zeros, for the species x year not seen in a site when seen in another.
ubmsdata <- read.table("DATA/SINDEX_ZERO_MS.txt", header = TRUE, sep = "\t") %>%
  filter(YEAR > 2020) %>%
  mutate(
    SPECIES = case_when(
      SPECIES == "Glaucopsyche sp." ~ "Glaucopsyche melanops",
      SPECIES == "Melanargia sp."   ~ "Melanargia lachesis",
      SPECIES == "Satyridae"   ~ "Pyronia bathseba",
      SPECIES == "Pyronia sp"   ~ "Pyronia bathseba",
      SPECIES == "Gonepteryx sp." ~ "Gonepteryx cleopatra",
      SPECIES == "Colias crocea" ~ "Colias croceus",
      SPECIES == "Pyrgus sp." ~ "Pyrgus malvoides",
      TRUE                          ~ SPECIES
    )
  )

list(unique(ubmsdata$SPECIES))

traits <- read.csv2("DATA/Traits.csv")[,1:9]
traits <- traits %>%
  rename("SPECIES" = Scientific.Name)

#ubmsdata <- ubmsdata %>%
#  left_join(traits, by = "SPECIES") 
#na.ubmsdata<-ubmsdata %>% filter(is.na(HSI))


# 1. Abundance matrix (SITE × SPECIES) 
abun <- ubmsdata %>%
  group_by(SITE_ID, YEAR, SPECIES) %>%
  summarise(abun = sum(SINDEX, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SPECIES, values_from = abun, values_fill = 0) %>%
  unite("site_year", SITE_ID, YEAR, sep = "_") %>%
  arrange(site_year) %>%                # opcional, solo para orden
  as.data.frame() %>%                   # <- convertir a data.frame
  column_to_rownames("site_year") %>%   # <- ahora sí quedan como rownames
  as.matrix()

abun <- abun[rowSums(abun) > 0, , drop = FALSE]

# 2. Trait Matrix (SPECIES × TRAITS) 
traits.m <- traits %>%
  filter(SPECIES %in% colnames(abun)) %>%
  distinct(SPECIES, .keep_all = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("SPECIES")

# order objects
common_sp <- intersect(rownames(traits.m), colnames(abun))
common_sp <- sort(common_sp)
traits.m <- traits.m[common_sp, , drop = FALSE]
abun   <- abun[, common_sp, drop = FALSE]

# clean DBs
abun.clean <- abun[rowSums(abun) > 0, , drop = FALSE]
abun.clean <- abun.clean[, colSums(abun.clean) > 0, drop = FALSE]
traits.m <- traits.m[rownames(traits.m) %in% colnames(abun.clean), , drop = FALSE]
traits.m  <- traits.m[,-1]
traits.m$SSI <- as.numeric(traits.m$SSI); traits.m$HSI <- as.factor(traits.m$HSI)
traits.m$Voltinism <- as.factor(traits.m$Voltinism); traits.m$Overwintering <- as.factor(traits.m$Overwintering)
traits.m$STI <- as.numeric(traits.m$STI); traits.m$Mobility <- as.factor(traits.m$Mobility)
traits.m$TAO <- as.numeric(traits.m$TAO)

# 3. Functional Distance (Gower: categorical and continous) & Rao's Q
stopifnot(identical(rownames(traits.m), colnames(abun.clean)))
gowerD <- gowdis(traits.m)   

rao <- dbFD(
  x = traits.m,
  a = abun.clean,
  calc.FRic = FALSE,
  calc.FDiv = FALSE,
  calc.CWM  = FALSE,
  stand.x   = TRUE,
  messages  = FALSE
)

# 4. Shannon diversity (taxonomic)
shannon <- diversity(abun.clean, index = "shannon")

# 5. Phylogenetic Diversity (Faith's PD)
tree <- read.tree("DATA/Butterflies_Europe_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)

sp_abun <- colnames(abun.clean)
sp_tree <- tree$tip.label

common_sp       <- intersect(sp_abun, sp_tree)
missing_in_tree <- setdiff(sp_abun, sp_tree)   # por si quieres inspeccionarlas

# Podar matriz de abundancias y árbol a las especies comunes
abun_phylo  <- abun.clean[, common_sp, drop = FALSE]
tree_pruned <- drop.tip(tree, setdiff(sp_tree, common_sp))

# Faith's PD por sitio-año
pd_res <- pd(
  samp = abun_phylo,
  tree = tree_pruned,
  include.root = TRUE
)


# 6. Final DB with functional diversity per site
fd.df <- tibble(
  site_year = names(rao$RaoQ),
  nbsp   = as.numeric(rao$nbsp),
  sing   = as.numeric(rao$sing.sp),
  RaoQ   = as.numeric(rao$RaoQ),
  Shannon = as.numeric(shannon[site_year]),
  PD      = as.numeric(pd_res$PD)
) %>%
  separate(site_year, into = c("SITE_ID", "YEAR"), sep = "_") %>%
  mutate(YEAR = as.numeric(YEAR)) %>%
  relocate(SITE_ID, YEAR)

write.table(fd.df, "DATA/fd_df.txt", sep = ";", row.names = FALSE)
