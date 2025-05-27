colors <- c("E13.5" = "#a9d0f5", 
            "E15.5" = "#e58b29",  
            "E13.5_ref" = "gray",
          "Control" = "#04a29d", 
            "Vcl cKO" = "#b2012e", 
           "E15.5_Neuroblast" = "#F97FBD",  
            "E13.5_BP" = "#6b853e",
            "E15.5_BranchA" = "#bca9f5", 
             "E13.5_BranchA" = "#be81f7" 

   )
results$name <- factor(results$name, levels = c("E13.5_BP", "E15.5_Neuroblast", "E13.5_BranchA", "E15.5_BranchA"))
results <- results %>%
  mutate(abs_log2meanfc = abs(log2meanfc)) %>%
  arrange(name)  # Sort by name to group similar categories

#################
terms <- list("E13.5_BP"= c("developmental growth involved in morphogenesis", "developmental cell growth"),
              "E15.5_Neuroblast"= c("developmental cell growth", "regulation of neurogenesis"),
              "E13.5_BranchA" = c("axon guidance", "neuron projection guidance"), 
              "E15.5_BranchA" = c("synaptic vesicle cycle","cell junction assembly",  "regulation of synapse structure or activity"))

plot.data <- data.frame()
for(i in names(terms)){ 
    temp <- subset(results, results$name == i & results$term %in% terms[[i]])
    plot.data <- rbind(plot.data, temp)    
}

###################
plot.data <- results

# Set a number of 'empty bars' to add at the end of each name group
empty_bar <- 1
to_add <- data.frame(matrix(NA, empty_bar * nlevels(as.factor(results$name)), ncol(results)))
colnames(to_add) <- colnames(results)
to_add$name <- rep(levels(as.factor(results$name)), each = empty_bar)
data <- rbind(plot.data, to_add)
data <- data %>% arrange(name)
data$id <- seq(1, nrow(data))

# Get the name and y position of each term label (above bars)
label_data_term <- data
number_of_bar <- nrow(label_data_term)
angle <- 90 - 360 * (label_data_term$id - 0.5) / number_of_bar  # Center the label in the bar
label_data_term$hjust <- ifelse(angle < -90, 1, 0)
label_data_term$angle <- ifelse(angle < -90, angle + 180, angle)

# Prepare data for name labels (below bars) and base lines
name_groups <- data %>%
  group_by(name) %>%
  summarise(start = min(id), end = max(id) - empty_bar) %>%  # Exclude empty bars for base lines
  rowwise() %>%
  mutate(title = mean(c(start, end)))  # Center of the group for name label

# Prepare name label data
label_data_name <- name_groups %>%
  mutate(
    label = name,
    y = -max(data$abs_log2meanfc, na.rm = TRUE) * 0.3,  # Position below bars
    angle = 90 - 360 * (title - 0.5) / number_of_bar,   # Use title (central id) for angle
    hjust = ifelse(angle < -90, 1, 0),
    angle = ifelse(angle < -90, angle + 180, angle)
  )

library(ggplot2)

max(data$pval,na.rm = T)

data$alph <- ifelse(data$pval < 0.001, 1, 
                    ifelse(data$pval < 0.01, 0.7, 0.4))

p <- ggplot(data, aes(x = id, y = abs_log2meanfc, fill = name, alpha = abs_log2meanfc*2)) +
  # Add bars
  geom_bar(stat = "identity", na.rm = TRUE) +
  # Set y-axis limits (small central circle)
  # ylim(-max(data$abs_log2meanfc, na.rm = TRUE) * 0.2, max(data$abs_log2meanfc, na.rm = TRUE) * 1) +
  # ylim(-max(data$abs_log2meanfc, na.rm = TRUE) * 0.2, max(data$abs_log2meanfc, na.rm = TRUE) * 1) +

  # Add term labels (inside bars) with curved alignment
  geom_text(data = label_data_term, 
            aes(x = id, 
                y = abs_log2meanfc * 0.2,  # Place labels at half the bar height
                label = term, 
                hjust = hjust, 
                angle = angle), 
            color = "black", fontface = "bold", alpha = 0.8, size = 2.5, 
            na.rm = TRUE, inherit.aes = FALSE) +
  # Add name labels (below bars) with curved alignment
  geom_text(data = label_data_name, 
            aes(x = title, 
                y = y, 
                label = label, 
                hjust = hjust, 
                angle = angle), 
            color = "black", fontface = "bold", alpha = 0.7, size = 0, 
            na.rm = TRUE, inherit.aes = FALSE) +
  # Add base lines for name groups
  geom_segment(data = name_groups, 
               aes(x = start, color = name, y = -0.1 * max(data$abs_log2meanfc, na.rm = TRUE), 
                   xend = end, yend = -0.1 * max(data$abs_log2meanfc, na.rm = TRUE)), alpha = 0.8, linewidth = 0.6, inherit.aes = FALSE) +
  # Apply polar coordinates
  coord_polar() +
  # Customize theme
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(0, 4), "cm"), 
      legend.text = element_text(size=6)
  ) + scale_fill_manual(values = colors) + scale_color_manual(values = colors)


# Set plot dimensions
options(repr.plot.width = 8, repr.plot.height = 9)
p
