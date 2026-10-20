library(ggplot2)

# Read the saved expression summary
gene_summary <- read.csv(
  "data/processed/gene_expression_summary.csv"
)

# Select the 10 largest decreases and increases
ranked_genes <- gene_summary[
  order(gene_summary$mean_difference_rpkm),
]

plot_data <- rbind(
  head(ranked_genes, 10),
  tail(ranked_genes, 10)
)

plot_data$direction <- ifelse(
  plot_data$mean_difference_rpkm > 0,
  "Higher in psoriasis",
  "Lower in psoriasis"
)

# Plot differences in group averages
difference_plot <- ggplot(
  plot_data,
  aes(
    x = reorder(gene_symbol, mean_difference_rpkm),
    y = mean_difference_rpkm,
    fill = direction
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(
    "Higher in psoriasis" = "#B44939",
    "Lower in psoriasis" = "#32749B"
  )) +
  labs(
    title = "Largest differences in average gene expression",
    subtitle = "Psoriasis minus healthy control; ranked by absolute RPKM difference",
    x = "Gene symbol",
    y = "Difference in mean RPKM",
    fill = NULL,
    caption = "Exploratory comparison; statistical significance has not been tested."
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(difference_plot)

# Create a folder for the finished figures
figure_folder <- "results/figures"

dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)

# Save the chart
ggsave(
  filename = file.path(figure_folder, "mean_expression_differences.png"),
  plot = difference_plot,
  width = 11,
  height = 8,
  units = "in",
  dpi = 300
)



# Connect to the database
con <- DBI::dbConnect(
  RSQLite::SQLite(),
  dbname = "data/processed/psoriasis.sqlite"
)

# Retrieve sample-level expression for three example genes
selected_expression <- DBI::dbGetQuery(con, "
  SELECT
    e.gene_symbol,
    e.sample_id,
    e.rpkm,
    s.disease_group
  FROM expression AS e
  JOIN samples AS s
    ON e.sample_id = s.sample_id
  WHERE e.gene_symbol IN ('S100A8', 'PI3', 'DCD');
")

DBI::dbDisconnect(con)

# Set group order
selected_expression$disease_group <- factor(
  selected_expression$disease_group,
  levels = c("healthy_control", "psoriasis")
)

# Plot each sample, with boxplots summarizing the distributions
sample_plot <- ggplot(
  selected_expression,
  aes(x = disease_group, y = log2(rpkm + 1), fill = disease_group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0, seed = 42),
    alpha = 0.4,
    size = 1
  ) +
  facet_wrap(~ gene_symbol, scales = "free_y") +
  scale_fill_manual(values = c(
    "healthy_control" = "#32749B",
    "psoriasis" = "#B44939"
  )) +
  scale_x_discrete(labels = c(
    "healthy_control" = "Control",
    "psoriasis" = "Psoriasis"
  )) +
  labs(
    title = "Sample-level expression of selected genes",
    subtitle = "Examples selected from the exploratory group comparison",
    x = NULL,
    y = "log2(RPKM + 1)",
    caption = "Each point represents one sample. Panels use separate y-axis scales."
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(sample_plot)

# Save the figure
ggsave(
  file.path(figure_folder, "selected_gene_distributions.png"),
  plot = sample_plot,
  width = 10,
  height = 5,
  dpi = 300
)

# Load the saved model results
de_results <- read.csv(
  "results/tables/exploratory_differential_expression.csv"
)

# Plot model effects against adjusted p-values
volcano_plot <- ggplot(
  de_results,
  aes(
    x = logFC,
    y = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
    color = result_group
  )
) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50"
  ) +
  scale_color_manual(values = c(
    "Does not meet both thresholds" = "grey75",
    "Higher in psoriasis" = "#B44939",
    "Lower in psoriasis" = "#32749B"
  )) +
  labs(
    title = "Exploratory differential expression in psoriasis",
    subtitle = "489 higher and 642 lower at the selected thresholds",
    x = "Difference in mean log2(RPKM + 1): psoriasis minus control",
    y = "-log10(BH-adjusted p-value)",
    color = NULL,
    caption = "Group-only limma model; results are not validated biomarkers."
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(volcano_plot)

# Save the plot
ggsave(
  file.path(figure_folder, "exploratory_volcano_plot.png"),
  plot = volcano_plot,
  width = 11,
  height = 7,
  dpi = 300
)


# Load saved pathway results
enrichment_results <- read.csv(
  "results/tables/pathway_enrichment_results.csv"
)

# Select the 10 significant pathways with smallest adjusted p-values
top_pathways <- enrichment_results[
  !is.na(enrichment_results$adjusted_p_value) &
    enrichment_results$adjusted_p_value < 0.05,
]

top_pathways <- head(
  top_pathways[order(top_pathways$adjusted_p_value), ],
  10
)

# Make pathway names readable and wrap long labels
top_pathways$label <- vapply(
  gsub(
    "_", " ",
    sub("^REACTOME_", "", top_pathways$pathway_name)
  ),
  function(x) paste(strwrap(x, width = 40), collapse = "\n"),
  character(1)
)

# Put the strongest result at the top
top_pathways$label <- factor(
  top_pathways$label,
  levels = rev(top_pathways$label)
)

pathway_plot <- ggplot(
  top_pathways,
  aes(x = fold_enrichment, y = label)
) +
  geom_point(
    aes(
      size = selected_pathway_genes,
      color = -log10(adjusted_p_value)
    )
  ) +
  scale_color_gradient(low = "#32749B", high = "#B44939") +
  labs(
    title = "Top enriched pathways in the psoriasis comparison",
    subtitle = "Top 10 by BH-adjusted p-value",
    x = "Fold enrichment",
    y = NULL,
    size = "Selected genes",
    color = "-log10\nadjusted p-value",
    caption = paste(
      "Background: 8,124 tested, annotated genes;",
      "selected: 745. Higher and lower expression combined."
    )
  ) +
  theme_minimal()

print(pathway_plot)

# Save the figure
ggsave(
  file.path(figure_folder, "pathway_enrichment_plot.png"),
  plot = pathway_plot,
  width = 12,
  height = 8,
  dpi = 300
)
