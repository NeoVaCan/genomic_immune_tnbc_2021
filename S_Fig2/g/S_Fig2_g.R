library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)

res_quantiseq <- read.csv(file = "combined_dfComplete.txt", sep = "\t", stringsAsFactors = TRUE)

res_quantiseq %>%
  gather(sample, fraction, -cell_type) %>%
    mutate(sample = factor(sample, levels=c("M1", "M15", "M16", "M18", "M5A", "M5B", "M5C", "M6A", "M6B", "M7", "M8A", "M8B", "M9", "M10C", "M10D", "M10E", "M11A", "M11B", "M12A"))) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_quantiseq)))

pdf('302_On-therapy_Pie_Charts_040521.pdf', width = 10, height = 10)

slices <- c(0.013222768, 0.017999577, 0.05175178, 0.024165596, 0.045865485, 0.019887824, 0.010997091, 0.014163427, 0.02215491, 0.006877781, 0.772913761)
lbls <- c("B cells", "Macrophages M1", "Macrophages M2",  "Monocytes",   "Neutrophils", "NK cells", "T cells CD4", "T cells CD8", "Tregs", "Dendritic cells", "Other")
pie(slices, labels = lbls, col=c("#451C87", "#b3b300", "#CE0648", "#2363C5", "#efd22b", "#AB4CA1", "green", "#DD8C24", "#24A9DD", "#ED6D42", "#808A8D"), main="On-Therapy")

dev.off()

pdf('302_Post-mortem_Pie_Charts_040521.pdf', width = 10, height = 10)

slices <- c(0.011312563, 0.009628436, 0.024419402, 0.001700684, 0.035586177, 0.014798608, 0.000258648, 0.013351561, 0.015422219, 0.010336971, 0.863184733)
lbls <- c("B cells", "Macrophages M1", "Macrophages M2",  "Monocytes",   "Neutrophils", "NK cells", "T cells CD4", "T cells CD8", "Tregs", "Dendritic cells", "Other")
pie(slices, labels = lbls, col=c("#451C87", "#b3b300", "#CE0648", "#2363C5", "#efd22b", "#AB4CA1", "green", "#DD8C24", "#24A9DD", "#ED6D42", "#808A8D"), main="Post-mortem")

dev.off()