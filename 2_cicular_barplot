library(ggplot2)
library(dplyr)
library(gghalves)

####### prepare data
results <- results %>%
  mutate(abs_log2meanfc = abs(log2meanfc)) %>%
  arrange(name)

plot.data <- results
plot.data$name <- factor(plot.data$name, levels = c("E13.5_BP", 
                                                    "E15.5_Neuroblast", 
                                                    "E13.5_BranchA", 
                                                    "E15.5_BranchA"))
wid.bar <- 0.5
plot.data$width <- wid.bar

######## add the empty bar and label data
# empty_bar <- 1
# to_add <- data.frame(matrix(NA, empty_bar * nlevels(as.factor(results$name)), ncol(results)))
# colnames(to_add) <- colnames(results)
# to_add$name <- rep(levels(as.factor(results$name)), each = empty_bar)
# to_add$width <- wid.bar/5 

# data <- rbind(plot.data, to_add)
data <- data %>% arrange(name)
data$id <- seq(1, nrow(data))

label_data_term <- data
number_of_bar <- nrow(data)
angle <- 90 - 360 * (label_data_term$id - 0.5) / number_of_bar
label_data_term$hjust <- ifelse(angle < -90, 1, 0)
label_data_term$angle <- ifelse(angle < -90, angle + 180, angle)

name_groups <- data %>%
  group_by(name) %>%
  summarise(start = min(id), end = max(id) - empty_bar) %>%
  rowwise() %>%
  mutate(title = mean(c(start, end)))

label_data_name <- name_groups %>%
  mutate(
    label = name,
    y = -max(data$abs_log2meanfc, na.rm = TRUE) * 0.3,
    angle = 90 - 360 * (title - 0.5) / number_of_bar,
    hjust = ifelse(angle < -90, 1, 0),
    angle = ifelse(angle < -90, angle + 180, angle)
  )


box.data <- plot.mat  ###################### using the data in net step.
data$label <- paste0(data$name, "_", data$term)
box.data$label <- paste0(box.data$anno, "_", box.data$term)

data.temp <- na.omit(data)
box.data <- subset(box.data, box.data$label %in% unique(data.temp$label))
box.data$id <- data.temp$id[match(box.data$label, data.temp$label)]
box.data$id <- factor(box.data$id)
data$id <- as.character(data$id)
box.data$id <- as.character(box.data$id)

# Plot with boxplots
p <- ggplot() +
  # Bar layer
  geom_bar(data = data, 
           aes(x = id, 
               y = abs_log2meanfc, 
               fill = name, 
               alpha = abs_log2meanfc * 0.5), 
           stat = "identity") +

  geom_half_boxplot(data = subset(box.data, group == "Control"),
                    aes(x = id,
                        y = path_sore,
                        outlier.stroke = 0,
                        fill = group),
                    alpha = 0.3,
                    width = 0.7,
                    outlier.shape = NA,  
                    outlier.size = 0,
                    outlier.stroke = 0,
                    errorbar.draw = FALSE,
                    side = 'r',
                    position = position_nudge(y = max(data$abs_log2meanfc) * -0.5)
                   ) + 

    geom_half_boxplot(data = subset(box.data, group == "Vcl cKO"), 
                      aes(x = id,
                        y = path_sore,
                        outlier.stroke = 0,
                        fill = group),
                      alpha = 0.3,
                      width = 0.7,
                      outlier.shape = NA,
                      outlier.size = 0,
                      outlier.stroke = 0,
                      errorbar.draw = FALSE,
                      side = 'l',
                      position = position_nudge(y = max(data$abs_log2meanfc) * -0.5)
                     ) + 
    geom_segment(data = name_groups, 
                   aes(x = start, 
                       y = -0.1 * max(data$abs_log2meanfc, na.rm = TRUE), 
                       xend = end, 
                       yend = -0.1 * max(data$abs_log2meanfc, na.rm = TRUE), 
                       color = name),
                 linewidth = 3,
                 alpha = 0.8) + 
    geom_text(data = label_data_term, 
                aes(x = id, 
                    y = abs_log2meanfc * 0.15,
                    label = term, 
                    hjust = hjust, 
                    angle = angle), 
                color = "black", 
                fontface = "bold", 
                alpha = 0.8, 
                size = 4.5, 
                na.rm = TRUE) +

    geom_text(data = label_data_term, 
            aes(x = id, 
                y = abs_log2meanfc * 0.15,
                label = term, 
                hjust = hjust, 
                angle = angle),
            color = "black",
            fontface = "bold", 
            alpha = 0.8, 
            size = 4.5, 
            na.rm = TRUE) +

    coord_polar() +
    theme_minimal() +
    theme(
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(0, 4), "cm"), 
        legend.text = element_text(size = 10)
    ) + 
    ylim(-max(data$abs_log2meanfc, na.rm = TRUE) * 0.5,
         max(data$abs_log2meanfc, na.rm = TRUE) * 1.2) +  
    scale_fill_manual(values = colors) + 
    scale_color_manual(values = colors)

options(repr.plot.width = 14, repr.plot.height = 12)
p
