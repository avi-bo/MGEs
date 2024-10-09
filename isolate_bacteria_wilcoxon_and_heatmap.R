# Clear the environment
rm(list = ls())

# Load required libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(pheatmap)
library(RColorBrewer)

# Function to perform habitat comparison analysis
perform_habitat_comparison <- function(data, habitat1, habitat2) {
  result_data <- data %>%
    filter(habitat %in% c(habitat1, habitat2)) %>%
    group_by(category, family) %>%
    summarize(
      p_value = wilcox.test(norm ~ habitat, exact = FALSE)$p.value,
      fold_change = median(norm[habitat == habitat1]) / median(norm[habitat == habitat2]),
      .groups = "drop"
    )
  
  return(result_data)
}

# Updated function to prepare transposed heatmap data
prepare_heatmap_data <- function(result_data, phylum_data, p_value_threshold = 0.05) {
  significant_results <- result_data %>% filter(p_value < p_value_threshold)
  
  heatmap_data <- result_data %>%
    left_join(significant_results, by = c("category", "family")) %>%
    select(category, family, fold_change.x, fold_change.y) %>%
    rename(fold_change = fold_change.x, significant = fold_change.y) %>%
    select(category, family, significant) %>%
    pivot_wider(names_from = family, values_from = significant) %>%
    as.data.frame()
  
  heatmap_data[is.na(heatmap_data)] <- 1
  
  # Set row names and remove category column
  row.names(heatmap_data) <- heatmap_data$category
  heatmap_data <- heatmap_data[, -1]
  
  # Order columns (families) by phylum
  phylum_order <- phylum_data[order(phylum_data$phylum), "family"]
  heatmap_data <- heatmap_data[, match(phylum_order, names(heatmap_data))]
  
  # Reorder rows (categories) as specified
  category_order <- c("Prophages - pfam", "Plasmids - pfam", "Transposons - pfam", "Defense systems - pfam",
                      "Prophages - intact", "Defense systems - intact")
  
  # Reorder the rows
  heatmap_data <- heatmap_data[match(category_order, rownames(heatmap_data)), ]
  
  return(heatmap_data)
}


# Updated function to create phylum annotation (now for columns)
create_phylum_annotation <- function(heatmap_data, phylum_data) {
  phylum_annotation <- data.frame(Phyla = phylum_data$phylum[match(colnames(heatmap_data), phylum_data$family)])
  rownames(phylum_annotation) <- colnames(heatmap_data)
  return(phylum_annotation)
}

# Function to create color palette for phyla
create_phylum_colors <- function(phylum_annotation) {
  unique_phyla <- unique(phylum_annotation$Phyla)
  palette <- colorRampPalette(brewer.pal(n = length(unique_phyla), name = "Set3"))
  colors <- palette(length(unique_phyla))
  annotation_colors <- setNames(colors, unique_phyla)
  return(list(Phyla = annotation_colors))
}

# Updated function to create heatmap
create_heatmap <- function(heatmap_data, phylum_annotation, annotation_colors, comparison_name, habitat1, habitat2) {
  bk1 <- c(seq(0, 0.9, by = 0.2), 0.999)
  bk2 <- c(1.01, seq(1.1, 2.2, by = 0.2))
  bk <- c(bk1, bk2)
  
  # Choose color palette based on comparison
  if (comparison_name == "PA_vs_NPA") {
    my_palette <- c(colorRampPalette(colors = c("#370000", "tomato1"))(n = length(bk1) - 1),
                    "white", "white",
                    c(colorRampPalette(colors = c("lightblue", "darkblue"))(n = length(bk2) - 1)))
  } else if (comparison_name == "PA_vs_SOIL") {
    my_palette <- c(colorRampPalette(colors = c("#26200f", "#BFB599"))(n = length(bk1) - 1),
                    "white", "white",
                    c(colorRampPalette(colors = c("lightblue", "darkblue"))(n = length(bk2) - 1)))
  } else if (comparison_name == "SOIL_vs_NPA") {
    my_palette <- c(colorRampPalette(colors = c("#370000", "tomato1"))(n = length(bk1) - 1),
                    "white", "white",
                    c(colorRampPalette(colors = c("#BFB599", "#26200f"))(n = length(bk2) - 1)))
  } else {
    stop("Invalid comparison name")
  }
  
  # Create a gaps_row vector to add space after the "pfam" categories
  gaps_row <- 4  
  
  p <- pheatmap(heatmap_data,
                scale = "none",
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                color = my_palette,
                breaks = bk,
                ylab = "category",
                xlab = "family",
                fontsize_row = 10,
                fontsize_col = 10,
                cellwidth = 20,
                cellheight = 20,
                border_color = "black",
                legend_breaks = c(0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.22, 1.42, 1.62, 1.82, 2.02),
                legend_labels = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.6-0.99", "", "1.01 - 1.2", "1.2-1.4", "1.4-1.6", "1.6-1.8", ">1.8"),
                display_numbers = FALSE,
                fontsize = 12,
                annotation_col = phylum_annotation,
                annotation_colors = annotation_colors,
                angle_col = 45,
                gaps_row = gaps_row,
                silent = F) 
  
  # Save the heatmap as a TIFF file
  jpeg(filename = paste0("path/to/output_folder/", habitat1, "_", habitat2, ".jpeg"),
       width = 280, height = 80, units = "mm", res = 600)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  dev.off()
  
  # Return the plot object in case it's needed for further manipulation
  return(p)
}

# Main analysis function
perform_analysis <- function(data, habitat1, habitat2, comparison_name, phylum_data) {
  # Perform comparison
  result_data <- perform_habitat_comparison(data, habitat1, habitat2)
  
  # Prepare heatmap data
  heatmap_data <- prepare_heatmap_data(result_data, phylum_data)
  
  # Create phylum annotation
  phylum_annotation <- create_phylum_annotation(heatmap_data, phylum_data)
  
  # Create phylum colors
  annotation_colors <- create_phylum_colors(phylum_annotation)
  
  # Create and save heatmap
  create_heatmap(heatmap_data, phylum_annotation, annotation_colors, comparison_name, habitat1, habitat2)
}


# Read and preprocess data
all_data <- read_tsv("path/to/file.tsv")
colnames(all_data)[4] <- "habitat"

# Update category names
all_data$category <- case_when(
  all_data$category == "defense_system-PF" ~ "Defense systems - pfam",
  all_data$category == "plasmid-PF" ~ "Plasmids - pfam",
  all_data$category == "prophage-PF" ~ "Prophages - pfam",
  all_data$category == "transposon-PF" ~ "Transposons - pfam",
  all_data$category == "defense_systems" ~ "Defense systems - intact",
  all_data$category == "prophages" ~ "Prophages - intact",
  TRUE ~ all_data$category
)

# Define phylum data
phylum_data <- data.frame(
  family = c("Bacillaceae", "Bacillaceae_G", "Burkholderiaceae", "Enterobacteriaceae", "Flavobacteriaceae", "Listeriaceae", "Microbacteriaceae", "Micrococcaceae", "Micromonosporaceae", "Nocardioidaceae", "Paenibacillaceae", "Pseudomonadaceae", "Pseudonocardiaceae", "Rhizobiaceae", "Rhodobacteraceae", "Sphingomonadaceae", "Streptomycetaceae", "Weeksellaceae", "Xanthomonadaceae"),
  phylum = c("Firmicutes", "Firmicutes", "Proteobacteria", "Proteobacteria", "Bacteroidota", "Firmicutes", "Actinobacteriota", "Actinobacteriota", "Actinobacteriota", "Actinobacteriota", "Firmicutes", "Proteobacteria", "Actinobacteriota", "Proteobacteria", "Proteobacteria", "Proteobacteria", "Actinobacteriota", "Bacteroidota", "Proteobacteria")
)

# Perform analyses
perform_analysis(all_data, "PA", "NPA", "PA_vs_NPA", phylum_data)
perform_analysis(all_data, "PA", "SOIL", "PA_vs_SOIL", phylum_data)
perform_analysis(all_data, "SOIL", "NPA", "SOIL_vs_NPA", phylum_data)
