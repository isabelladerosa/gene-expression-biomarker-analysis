# Connect to the existing psoriasis database
con <- DBI::dbConnect(
  RSQLite::SQLite(),
  dbname = "data/processed/psoriasis.sqlite"
)

# Calculate average expression for each gene in each group
gene_summary <- DBI::dbGetQuery(con, "
  SELECT
    e.gene_symbol,
    AVG(CASE WHEN s.disease_group = 'healthy_control'
        THEN e.rpkm END) AS mean_control_rpkm,
    AVG(CASE WHEN s.disease_group = 'psoriasis'
        THEN e.rpkm END) AS mean_psoriasis_rpkm
  FROM expression AS e
  JOIN samples AS s
    ON e.sample_id = s.sample_id
  GROUP BY e.gene_symbol
  ORDER BY e.gene_symbol;
")

# Inspect the results
head(gene_summary)
nrow(gene_summary)

# Close this script's database connection
DBI::dbDisconnect(con)

# Save the SQL-generated expression summary
write.csv(
  gene_summary,
  "data/processed/gene_expression_summary.csv",
  row.names = FALSE
)



# Calculate the difference in average RPKM
gene_summary$mean_difference_rpkm <-
  gene_summary$mean_psoriasis_rpkm -
  gene_summary$mean_control_rpkm

# Calculate the ratio where control expression is above zero
# Leave other ratios missing to avoid dividing by zero
gene_summary$expression_ratio <- NA_real_

has_control_expression <- gene_summary$mean_control_rpkm > 0

gene_summary$expression_ratio[has_control_expression] <-
  gene_summary$mean_psoriasis_rpkm[has_control_expression] /
  gene_summary$mean_control_rpkm[has_control_expression]

# Preview genes with the largest increases in average RPKM
top_increases <- gene_summary[
  order(gene_summary$mean_difference_rpkm, decreasing = TRUE),
]

head(top_increases, 10)



# Rank genes by the largest decreases in average RPKM
top_decreases <- gene_summary[
  order(gene_summary$mean_difference_rpkm),
]

# Preview the top 10
head(top_decreases, 10)

# Save the updated summary, including differences and ratios
write.csv(
  gene_summary,
  "data/processed/gene_expression_summary.csv",
  row.names = FALSE
)


# Create a folder for analysis results
table_folder <- "results/tables"

dir.create(table_folder, recursive = TRUE, showWarnings = FALSE)

# Rank genes by differences in average expression
ranked_genes <- gene_summary[
  order(gene_summary$mean_difference_rpkm),
]

# Save the 10 largest decreases
write.csv(
  head(ranked_genes, 10),
  file.path(table_folder, "top_10_mean_expression_decreases.csv"),
  row.names = FALSE
)

# Save the 10 largest increases, largest first
write.csv(
  head(
    gene_summary[
      order(gene_summary$mean_difference_rpkm, decreasing = TRUE),
    ],
    10
  ),
  file.path(table_folder, "top_10_mean_expression_increases.csv"),
  row.names = FALSE
)


