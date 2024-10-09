rm(list=ls())

library(readr)
library(ggplot2)
library(dplyr)
library(dunn.test)
library(ggsignif)
library(ggpubr)
library(multcompView)
library(rcompanion)
library(FSA)

# Function to clean habitat names
clean_habitat_names <- function(data) {
  data$ext_hab <- data$ext_hab %>%
    tolower() %>%
    gsub(" ", "_", .) %>%
    gsub("[^a-z0-9_]", "", .)
  return(data)
}

# Function to generate label dataframe
generate_label_df <- function(data, category) {
  # Check if there are at least two groups with sufficient data
  group_counts <- table(data$ext_hab)
  if (sum(group_counts >= 2) < 2) {
    return(NULL)  
  }
  
  # Perform Kruskal-Wallis test
  kruskal <- tryCatch({
    kruskal.test(norm ~ ext_hab, data = data)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(kruskal) || kruskal$p.value >= 0.05) {
    # If Kruskal-Wallis test fails or is not significant, assign same letter to all groups
    return(data.frame(groups = unique(data$ext_hab),
                      letters = rep("a", length(unique(data$ext_hab)))))
  }
  
  # Perform Dunn's test
  dunn_result <- tryCatch({
    dunn.test(data$norm, data$ext_hab, method = "bh", kw = FALSE, table = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(dunn_result) || length(dunn_result$comparisons) == 0) {
    return(data.frame(groups = unique(data$ext_hab),
                      letters = rep("a", length(unique(data$ext_hab)))))
  }
  
  # Create a data frame of comparisons and p-values
  comparison_df <- data.frame(
    comparison = dunn_result$comparisons,
    p.value = dunn_result$P.adjusted
  )
  
  # Generate compact letter display
  cld <- tryCatch({
    cldList(comparison = comparison_df$comparison,
            p.value = comparison_df$p.value,
            threshold = 0.01)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(cld)) {
    return(data.frame(groups = unique(data$ext_hab),
                      letters = rep("a", length(unique(data$ext_hab)))))
  }
  
  #print(cld)
  
  # Convert to data frame
  cld_df <- data.frame(groups = cld$Group,
                       letters = cld$Letter)
  
  return(cld_df)
}

# Analysis function
analyze <- function(data) {
  print(unique(data$pf_category))
  c <- unique(data$pf_category)
  
  # Clean habitat names
  data <- clean_habitat_names(data)
  
  # Remove NA values
  data <- data %>% filter(!is.na(norm) & !is.na(ext_hab))
  
  cols <- c("phyllo"= "#9DCA9D",
            "roots" = "#9DCA9D",
            "green_algae" = "#D5EFD5",
            "red_algae" = "#D5EFD5",
            "soil" = "#D9B2B2",
            "deep_subsurface" = "#D9B2B2",
            "skin" = "#F8C8DC",
            "nasal_cavity" = "#F8C8DC",
            "oral_cavity" = "#F8C8DC",
            "stomach" = "#F8C8DC",
            "large_intestine" = "#F8C8DC",
            "vagina" = "#F8C8DC",
            "fungi" = "#FFF9BF",
            "insects" = "#f8a19a",
            "birds" = "#ABBDCA",
            "fish" = "#EFBB83",
            "wetlands" = "#B2DFFF",
            "river" = "#B2DFFF",
            "lake" = "#B2DFFF",
            "groundwater" = "#B2DFFF",
            "hot_water" = "#B2DFFF",
            "coastal" = "#B2DFFF",
            "intertidal_zone" = "#B2DFFF",
            "oceanic" = "#B2DFFF",
            "hydrothermal_vents" = "#B2DFFF",
            "saline" = "#B2DFFF"
  )
  
  data$ext_hab <- factor(data$ext_hab, levels = names(cols))
  
  # Generate letter labels
  label_df <- generate_label_df(data, c)
  
  if (!is.null(label_df)) {
    # Merge letter labels with original data
    data_with_letters <- merge(data, label_df, by.x = "ext_hab", by.y = "groups", all.x = TRUE)
    # Calculate mean norm for each group to position letters
    data_summary <- data %>%
      group_by(ext_hab) %>%
      summarise(mean_norm = mean(norm, na.rm = TRUE),
                max_norm = max(norm, na.rm = TRUE))
    
    data_with_letters <- merge(data_with_letters, data_summary, by = "ext_hab")
    data_with_letters <- data_with_letters %>%
      mutate(text_y = case_when(
        pf_category == "defense" ~ 0.058,
        pf_category == "prophage" ~ 0.038,
        pf_category == "plasmid" ~ 0.0056,
        pf_category == "transposon" ~ 0.037,
        TRUE ~ NA_real_ 
      ))
    
    # Create the plot
    p <- ggplot(data_with_letters, aes(x = ext_hab, y = norm, fill = ext_hab)) +
      geom_boxplot(outlier.size=-1, show.legend = FALSE) +
      scale_fill_manual(values = cols) +
      xlab("Habitat") + ylab("normalized pfam values") +
      theme_minimal() +
      theme(text = element_text(size=20), 
            axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
            axis.text.y = element_text(size = 20),
            axis.title.y = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
      geom_text(aes(y = text_y, label = letters), 
            vjust = -0.5, size = 10, position = position_nudge(x = -0.1)) 
  } else {
    # Create plot without letters if statistical tests fail
    p <- ggplot(data, aes(x = ext_hab, y = norm, fill = ext_hab)) +
      geom_boxplot(outlier.size=-1) +
      scale_fill_manual(values = cols) +
      labs(title = unique(data$pf_category)) +
      xlab("Habitat") + ylab("normalized pfam values") +
      theme_minimal() +
      theme(text = element_text(size=20), 
            axis.text.x = element_text(angle = 45, hjust = 1, size = 20))
  }
  
  if(unique(data$pf_category) == "defense") {
    p <- p + coord_cartesian(ylim = c(0, 0.065))
  } else if(unique(data$pf_category) == "prophage") {
    p <- p + coord_cartesian(ylim = c(0, 0.04)) 
  } else if(unique(data$pf_category) == "plasmid") {
    p <- p + coord_cartesian(ylim = c(0, 0.006))
  } else if(unique(data$pf_category) == "transposon") {
    p <- p + coord_cartesian(ylim = c(0, 0.04))
  }
  
  list(plot = p)  
}

# Read in data 
df <- read_tsv("path/to/file.tsv")
df["ext_hab"][df["ext_hab"] == "Hot (42-90C)"] <- "hot_water"

# Clean habitat names in the main dataframe
df <- clean_habitat_names(df)

# filter out metagenomes with less then 10000 pfam domains
filter_df <- df[df$total_pfam > 10000, ]

# Split data
split_df <- split(filter_df, filter_df$pf_category)

# Apply function
results <- lapply(split_df, analyze)

# access to the plots:
results$defense$plot
results$prophage$plot
results$plasmid$plot
results$transposon$plot
