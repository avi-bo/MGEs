# Clear the environment
rm(list = ls())

# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(pheatmap)

# Read data
all_data <- read_tsv("path/to/file.tsv")
# Rename 'CasFinder' to 'CRISPR' in 'type' column
all_data$type[all_data$type == "CasFinder"] <- "CRISPR"
# Filter data to include only 'NPA' and 'PA' labels
all_data <- all_data %>% filter(label_AvP %in% c("NPA", "PA"))

# Define phylum data mapping for bacterial families
phylum_data <- tibble(
  family = c("Bacillaceae", "Bacillaceae_G", "Burkholderiaceae", "Enterobacteriaceae", 
             "Flavobacteriaceae", "Listeriaceae", "Microbacteriaceae", "Micrococcaceae", 
             "Micromonosporaceae", "Nocardioidaceae", "Paenibacillaceae", "Pseudomonadaceae", 
             "Pseudonocardiaceae", "Rhizobiaceae", "Rhodobacteraceae", "Sphingomonadaceae", 
             "Streptomycetaceae", "Weeksellaceae", "Xanthomonadaceae"),
  phylum = c("Firmicutes", "Firmicutes", "Proteobacteria", "Proteobacteria", "Bacteroidota", 
             "Firmicutes", "Actinobacteriota", "Actinobacteriota", "Actinobacteriota", "Actinobacteriota", 
             "Firmicutes", "Proteobacteria", "Actinobacteriota", "Proteobacteria", "Proteobacteria", 
             "Proteobacteria", "Actinobacteriota", "Bacteroidota", "Proteobacteria")
)

# Merge phylum information with the main dataset
all_data <- left_join(all_data, phylum_data, by = "family")


# Calculate the number of unique genome IDs per family and label_AvP category
family_lengths <- all_data %>%
  group_by(family, label_AvP) %>%
  summarise(family_len = n_distinct(genome_id), .groups = "drop")

# Merge family length data with the main dataset
all_data <- left_join(all_data, family_lengths, by = c("family", "label_AvP"))

# Count occurrences of each type per genome ID
df_counted <- all_data %>%
  count(genome_id, type, name = "type_len")

# Pivot data to wide format to get counts per type
df_pivoted <- df_counted %>%
  pivot_wider(names_from = type, values_from = type_len, values_fill = 0)

# Transform back to long format
df_long <- df_pivoted %>% 
  pivot_longer(
    cols = -c(genome_id), # exclude these columns
    names_to = "types",
    values_to = "type_len"
  )

# Select specific columns for merging
all_data_selected <- all_data %>%
  select(genome_id, family, label_AvP, phylum, family_len)

# Merge df_long with the selected columns from all_data and remove duplicate rows
df_merged <- left_join(df_long, all_data_selected, by = "genome_id")
df_merged <- df_merged %>% distinct()

# Perform Wilcoxon test
stat_data <- df_merged %>%
  group_by(types, family) %>%
  summarize(
    p_value = wilcox.test((type_len / family_len) ~ label_AvP, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>% 
  mutate(adjusted_p_value = p.adjust(p_value, method = "BH"))

# Summarize type occurrences
df_tmp <- df_merged %>%
  group_by(types,family,label_AvP,phylum,family_len) %>%
  summarise(
    family = first(family),
    label_AvP = first(label_AvP),
    phylum = first(phylum),
    family_len = first(family_len),
    type_len = sum(type_len)
  )

# Compute fold change between PA and NPA
df_fold_change <- df_tmp %>%
  group_by(family, types, label_AvP, phylum) %>%
  summarise(mean_len = median(type_len / family_len, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = label_AvP, values_from = mean_len, values_fill = list(NPA = NA, PA = NA)) %>%
  mutate(fold_change = case_when(
    NPA == 0 & PA != 0 ~ 30,
    PA == 0 & NPA != 0 ~ 0.001,
    TRUE ~ PA / NPA
  ))

# Merge fold change data with statistical results
df_fold_change <- df_fold_change %>% 
  left_join(stat_data, by = c("types", "family"))

# Mask non-significant results
df_fold_change$fold_change <- ifelse(
  is.na(df_fold_change$adjusted_p_value) | df_fold_change$adjusted_p_value > 0.05,
  NaN,
  df_fold_change$fold_change
)

# Pivot data for heatmap
df_wide <- df_fold_change %>%
  select(family, phylum, types, fold_change) %>%
  pivot_wider(names_from = types, values_from = fold_change, values_fill = NA)

df_wide <- df_wide[order(df_wide$phylum),]
df_wide <- df_wide %>% select_if(~ !all(is.na(.)))

heatmap_data <- df_wide[, -c(1,2)]  # Exclude family_label for matrix

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_data)
heatmap_matrix[is.na(heatmap_matrix)] <- 1  # Replace NA with 0

# Set row names
row.names(heatmap_matrix) <- df_wide$family

# Convert hex colors to RGB for colorRampPalette
start_color_1 <- "#8b0000"
end_color_1 <- "#ff6347"
start_color_2 <- "lightblue"
end_color_2 <- "darkblue"

# Calculate the midpoint color for dividing the lightblue to darkblue scale
rgb_start <- col2rgb(start_color_2)
rgb_end <- col2rgb(end_color_2)
rgb_mid <- (rgb_start + rgb_end) / 2
end_color_21 <- rgb(rgb_mid[1], rgb_mid[2], rgb_mid[3], maxColorValue = 255)

bk1 <- c(0.001)
bk2 <- seq(0.08, 0.99, by = 0.01)
bk3 <- c(1)
bk4 <- seq(1.01, 10, by = 0.5)
bk5 <- seq(11, 27, by = 5)
bk6 <- c(30)

bk <- c(bk1, bk2, bk3, bk4, bk5, bk6)

# Create color ramps
colors_1 <- c("#3f0000", colorRampPalette(c(start_color_1, end_color_1))(n = length(bk2) - 1)) # From 0.001 to 0.99
colors_2 <- colorRampPalette(c(start_color_2, end_color_21))(19) # From 1.01 to 10 by 0.5
colors_3 <- colorRampPalette(c(end_color_21, end_color_2))(3) # From 11 to 27 by 5

# Add white for value 1
colors <- c(colors_1, "white", colors_2, colors_3, "#000067")


# Now use these colors and breaks in your heatmap
pheatmap(heatmap_matrix,
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk,
         scale = "none",
         #main = "Fold Change between PA and NPA",
         color = colors,
         show_rownames = TRUE,
         display_numbers = F,
         #border_color = "black",
         cellheight = 20,
         legend_breaks = c(0.001, 0.08, 0.99, 1, 10, 27, 30),
         fontsize = 15
)
