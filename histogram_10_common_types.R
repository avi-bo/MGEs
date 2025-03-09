# Clear the environment
rm(list = ls())

# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# Read data
all_data <- read_tsv("path/to/file.tsv")
all_data$type[all_data$type == "CasFinder"] <- "CRISPR"
all_data <- all_data %>% filter(label_AvP %in% c("NPA", "PA"))

# Define phylum data
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

# Merge phylum data
all_data <- left_join(all_data, phylum_data, by = "family")

# Calculate family_len efficiently
family_lengths <- all_data %>%
  group_by(family, label_AvP) %>%
  summarise(family_len = n_distinct(genome_id), .groups = "drop")

all_data <- left_join(all_data, family_lengths, by = c("family", "label_AvP"))

# Calculate type count per family and label_AvP
df <- all_data %>%
  group_by(family, label_AvP, phylum, family_len, type) %>%
  summarise(type_len = n(), .groups = "drop")

# Ensure all type-family-label_AvP combinations exist
df_complete <- df %>%
  select(family, phylum, type, label_AvP, type_len, family_len) %>%
  complete(family, type, label_AvP, fill = list(type_len = 0, family_len = 0))

# Get unique family-label_AvP-phylum combinations
t_m <- df %>% select(family, label_AvP, phylum) %>% distinct()

# Merge the data and filter
merged_df <- left_join(df_complete, t_m, by = c("family", "label_AvP"))
  

# Compute fold change
df_fold_change <- merged_df %>%
  group_by(family, type, label_AvP, phylum.x) %>%
  summarise(mean_len = median(type_len / family_len, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = label_AvP, values_from = mean_len, values_fill = list(NPA = NA, PA = NA)) %>%
  mutate(fold_change = case_when(
    is.na(NPA) & !is.na(PA) ~ 30,
    is.na(PA) & !is.na(NPA) ~ 0.01,
    TRUE ~ PA / NPA
  ))



# List to store plots
plot_list <- list()

# Loop through each family and generate plots
for (fam in unique(df_fold_change$family)) {
  
  # Select top 10 most common types in PA or NPA within the family
  top_10_types <- df_fold_change %>%
    filter(family == fam) %>%
    pivot_longer(cols = c(PA, NPA), names_to = "category", values_to = "value") %>%
    group_by(type) %>%
    summarise(total_value = sum(value, na.rm = TRUE)) %>%
    top_n(10, total_value) %>%
    pull(type)
  
  # Filter the dataframe for only the top 10 types in the current family
  top_types_df <- df_fold_change %>%
    filter(family == fam, type %in% top_10_types) %>%
    pivot_longer(cols = c(PA, NPA), names_to = "category", values_to = "value")
  
  # Skip families with no data
  if (nrow(top_types_df) == 0) next
  
  # Generate the plot
  p <- ggplot(top_types_df, aes(x = reorder(type, -value), y = value, fill = category)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste(fam), x = "Type", y = "occurrences/genome", fill = "Category") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Store the plot in the list
  plot_list[[fam]] <- p
}

# Combine all plots in a grid
final_plot <- wrap_plots(plot_list) #+ plot_annotation(title = "Top 10 Most Common Types in PA & NPA Across Families")

# Print the final plot
print(final_plot)
