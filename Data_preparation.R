################################################################################

## Load packages ##
library(readr)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(factoextra)
library(vegan)
library(viridis)
library(ggfortify)
library(plotly)
library(dplyr)
library(tidyr)
library(janitor)
library(readxl)
library(forcats)
library(ggplot2)
library(ggdendro)
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)     
library(tidytext)
library(ggridges)

################################################################################
## Load data ##

# Set working directory
setwd("C:/Users/kleij057/OneDrive - Wageningen University & Research/DATA/SedaDNA clean code")

# Load sedaDNA data
nepaldat <- readr::read_tsv("260303_taxa_out.tsv")

# BUG: This needs to be moved from indexing to colnames
#nepaldat <- nepaldat[,c(1:4,7:70)]

#Load own sample labels
labels <- read.csv2("Lib_id2sample_id_new.csv", sep = ",")

# Add own labels to your data
# Rename the column
names(nepaldat)[names(nepaldat) == "library"] <- "Library_id"

# Reformat the IDs
nepaldat$Library_id <- gsub("_", ".", nepaldat$Library_id)
relab <- merge(labels, nepaldat, by = "Library_id")

# Now we remove negative controls
dat_exneg <- relab[!relab$Sample_id == "Extraction negative control",]
dat_exneg <- dat_exneg[!dat_exneg$Sample_id == "Library negative control",]

labels2 <- dat_exneg$Sample_id
labels2 <- as.data.frame(labels2)

labels2$loc <- sub("[0-9]+", "", labels2$labels2)
labels2$depth <- as.numeric(sub("[A-Za-z]+", "", labels2$labels2))

dat_exneg$Sample_id <- labels2$depth
names(dat_exneg)[names(dat_exneg) == "Sample_id"] <- "depth"
names(dat_exneg)[names(dat_exneg) == "Site_name"] <- "loc"
dat_exneg$loc <- labels2$loc

sedadna <- dat_exneg

## remove genus duplicates (these are cause because there is multiple species matches.
# so if you are working on species level, ignore this step)
sedadna <- sedadna[
  !duplicated(sedadna[, c("loc", "depth", "genus_name")]),
]

################################################################################
## Set your thresholds for what is included ##

# Only select your phyla of interest
sedadna <- subset(sedadna, phylum_name == "Streptophyta")

# Minimum amount of read counts
adna <- sedadna[sedadna$sum_genus_count >= 200,]

# Optional but often cleaner: only keep genera that occur in >X samples
adna <- adna %>%
  group_by(genus_name) %>%
  filter(n_distinct(depth) >= 3)


################################################################################
## Inspect and select bleedover taxa ##

#use unfiltered for readcounts here
#subset into columns you need
plotdna <- sedadna[,c("loc", "depth","phylum_name","order_name","family_name","genus_name",     
                      "taxa_name","sum_genus_count")]

## --- 1) Order genus by taxonomy: Phylum -> Order -> Family -> Genus ---

genus_levels <- plotdna %>%
  distinct(phylum_name, order_name, family_name, genus_name) %>%
  arrange(phylum_name, order_name, family_name, genus_name) %>%
  pull(genus_name)

plotdna <- plotdna %>%
  mutate(genus_name = factor(genus_name, levels = genus_levels, ordered = TRUE))


data <- plotdna[,c("loc", "depth", "order_name", "genus_name","sum_genus_count")]

# But do select only genera that you will have in your final data after filtering
full_taxa_list <- levels(as.factor(adna$genus_name))

data <- data %>%
filter(genus_name %in% full_taxa_list)


#One loc only (here MO)
data <- subset(data,  loc == "MO")
data$identifier<- paste(data$depth, data$loc, sep = "_")



#######
# Now only flag if in shared order
#######

# -------------------------
# 0) Genus → order lookup
# -------------------------
genus_order <- data %>%
  select(genus_name, order_name) %>%
  distinct()

# Safety check: ensure each genus maps to exactly one order
bad_map <- genus_order %>%
  count(genus_name) %>%
  filter(n > 1)

if (nrow(bad_map) > 0) {
  stop("Some genera map to multiple families. Fix genus-order mapping before proceeding.")
}

# -------------------------
# 1) Make wide community matrix (sum duplicates)
# -------------------------
comm_wide <- data %>%
  group_by(identifier, genus_name) %>%
  summarise(sum_genus_count = sum(sum_genus_count), .groups = "drop") %>%
  pivot_wider(
    names_from  = genus_name,
    values_from = sum_genus_count,
    values_fill = 0
  )

comm_mat <- comm_wide %>%
  select(-identifier) %>%
  as.matrix()

genera <- colnames(comm_mat)

# -------------------------
# 2) Pairwise Pearson tests (r, r2, p)
# -------------------------
results_list <- vector("list", length = (length(genera) * (length(genera) - 1)) / 2)
k <- 1

for (i in 1:(ncol(comm_mat) - 1)) {
  for (j in (i + 1):ncol(comm_mat)) {
    test <- cor.test(comm_mat[, i], comm_mat[, j], method = "pearson")
    
    r <- unname(test$estimate)
    
    results_list[[k]] <- data.frame(
      genus_1 = genera[i],
      genus_2 = genera[j],
      r       = r,
      r2      = r^2,
      p_value = test$p.value
    )
    k <- k + 1
  }
}

results <- bind_rows(results_list) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# -------------------------
# 3) Add order info + filter
# -------------------------
results_fam <- results %>%
  left_join(genus_order, by = c("genus_1" = "genus_name")) %>%
  rename(order_1 = order_name) %>%
  left_join(genus_order, by = c("genus_2" = "genus_name")) %>%
  rename(order_2 = order_name)

strong_same_order <- results_fam %>%
  filter(r2 > 0.8, order_1 == order_2) %>%
  arrange(desc(r2))

strong_same_order


#######
# Now plot the combo's
#######

library(ggplot2)
library(dplyr)

pairs_to_plot <- strong_same_order  # or strong_same_order_sig

plot_pair <- function(genus1, genus2, r2, p_value, p_adj, order_name = NA_character_) {
  dfp <- comm_wide %>%
    transmute(
      identifier = identifier,
      x = .data[[genus1]],
      y = .data[[genus2]]
    )
  
  ggplot(dfp, aes(x = x, y = y)) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_classic() +
    labs(
      x = genus1,
      y = genus2,
      title = paste0(genus1, " vs ", genus2,
                     if (!is.na(order_name)) paste0("  (", order_name, ")") else "")
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      hjust = 1.05, vjust = 1.2,
      label = sprintf("R² = %.3f\np = %.2g\np_adj = %.2g", r2, p_value, p_adj)
    )
}

pdf("flagged_pairs_scatterplots_0.8.pdf", width = 7, height = 6)

for (i in seq_len(nrow(pairs_to_plot))) {
  p <- plot_pair(
    genus1  = pairs_to_plot$genus_1[i],
    genus2  = pairs_to_plot$genus_2[i],
    r2      = pairs_to_plot$r2[i],
    p_value = pairs_to_plot$p_value[i],
    p_adj   = pairs_to_plot$p_adj[i],
    order_name  = pairs_to_plot$order_1[i]
  )
  print(p)
}

dev.off()


# Interpret and select Based on model fit & plots, bleed over taxa:

taxa_bleed <- c("Rubroshorea")

plotdna_clean <- adna %>%
  filter(!genus_name %in% taxa_bleed)

adna <- plotdna_clean


################################################################################
## Identify and remove outliers ##

# make label per observation
# taxa+site+depth

adna$identifier<- paste(adna$depth, adna$loc, adna$taxa_name, sep = "_")

# now subset to simple df
adna_sim <- adna[,c("identifier" , "depth", "loc", "order_name", "family_name", "genus_name", "genus_count_assembly", "sum_genus_count",     
                    "deam5p_frac_genus","deam3p_frac_genus" )]

# subset to per location if needed, rest of analyses will be done on MO only
adnaMO <- subset(adna_sim, adna$loc =="MO")

adna <- adnaMO


# 1. Reshape to long format: one row per depth × measure
plotdat <- adna %>%
  select(depth, identifier, deam3p_frac_genus, deam5p_frac_genus) %>%
  pivot_longer(
    cols = c(deam3p_frac_genus, deam5p_frac_genus),
    names_to = "end",
    values_to = "frac"
  )

# 2. Make interactive scatter plot to check if visually ok
figMO <- plot_ly(
  data = plotdat,
  x = ~depth,
  y = ~frac,
  type = "scatter",
  mode = "markers",
  color = ~end,  # one color for each column
  text = ~paste(
    "identifier:", identifier,
    "<br>depth:", depth,
    "<br>", end, ":", round(frac, 3)
  ),
  hoverinfo = "text"
) %>%
  layout(
    xaxis = list(title = "Depth"),
    yaxis = list(title = "Deamination fraction (genus)"),
    legend = list(title = list(text = "End"))
  )

figMO

### weighted, flags any outliers with an X. adapt outlier definition if needed ###
plot_dat <- adna %>%
  select(depth, identifier, genus_count_assembly,
         deam3p_frac_genus, deam5p_frac_genus) %>%
  pivot_longer(
    cols = c(deam3p_frac_genus, deam5p_frac_genus),
    names_to = "end",
    values_to = "frac"
  ) %>%
  group_by(depth, end) %>%
  mutate(
    # weighted mean
    w_mean = weighted.mean(frac, genus_count_assembly, na.rm = TRUE),
    
    # weighted variance → sd
    w_var  = sum(genus_count_assembly * (frac - w_mean)^2, na.rm = TRUE) /
      sum(genus_count_assembly, na.rm = TRUE),
    w_sd   = sqrt(w_var),
    
    # weighted Z-score
    z_w    = (frac - w_mean) / w_sd,
    
    # outlier definition
    outlier = abs(z_w) > 3
  ) %>%
  ungroup()

figMO <- plot_ly() %>%
  add_markers(
    data = plot_dat,
    x = ~depth,
    y = ~frac,
    color = ~end,
    symbols = c("circle", "x"),
    symbol = ~outlier,
    size = ~ifelse(outlier, 12, 8),
    text = ~paste(
      "identifier:", identifier,
      "<br>depth:", depth,
      "<br>", end, ":", round(frac, 3),
      "<br>Weighted Z:", round(z_w, 2),
      "<br>Outlier:", outlier
    ),
    hoverinfo = "text"
  ) %>%
  layout(
    xaxis = list(title = "Depth"),
    yaxis = list(title = "Deamination fraction (genus)"),
    legend = list(title = list(text = "Series"))
  )

figMO

# Write flagged outliers into file
outlier_table <- plot_dat %>%
  filter(outlier) %>%
  arrange(depth, end, desc(z_w)) %>%
  select(
    depth,
    identifier,
    end,
    frac,
    genus_count_assembly,
    w_mean,
    w_sd,
    z_w
  )

outlier_table <- subset(outlier_table, outlier_table$z_w < 0)
write.csv(outlier_table, "outlier_table.csv", row.names = FALSE)

# Now remove outliers from data set
outliers <- read.csv2("outlier_table.csv", sep = ",")

outliers <- outliers$identifier

plotdna <- adna %>%
  filter(!identifier %in% outliers)

################################################################################
# Add a relative abundance column and corresponding size bins

#transfer to relative abundances
plotdna <- plotdna %>%
  group_by(loc, depth) %>%
  mutate(rel_abundance = 100 * sum_genus_count / sum(sum_genus_count)) %>%
  ungroup()


#basic stats inspection to select binning
histdat <- plotdna[plotdna$rel_abundance < 20,]
hist(histdat$rel_abundance)

# make bins (make sure they are a factor)
# Define breaks
breaks <- c(0, 0.5, 1, 3, 5, 10, 20, 50, Inf)

# Define labels (in the desired order)
labels <- c("<0.5%","0.5-1%", "1-3%", "3-5%", "5-10%", "10-20%", "20-50%", ">50%")

# Add ordered factor column 'bin'
plotdna$bin <- cut(
  plotdna$rel_abundance,
  breaks = breaks,
  labels = labels,
  right = TRUE,
  include.lowest = TRUE,
  ordered_result = TRUE
)


################################################################################
# now write a .csv with cleaned data that can be loaded into new scripts
write.csv(plotdna, "cleandna.csv", row.names = FALSE)

