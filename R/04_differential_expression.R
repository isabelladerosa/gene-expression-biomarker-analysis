library(limma)

# Set your project folder
project_folder <- "."

# Load expression values
expression_data <- read.delim(
  gzfile(file.path(
    project_folder,
    "data/raw/GSE54456_RPKM_samples.txt"
  )),
  check.names = FALSE
)

names(expression_data)[1] <- "gene_symbol"

# Load sample information
samples <- read.csv(
  file.path(project_folder, "data/processed/sample_metadata.csv"),
  stringsAsFactors = FALSE
)

# Match sample information to expression-column order
sample_ids <- names(expression_data)[-1]

stopifnot(
  !anyDuplicated(samples$sample_id),
  setequal(sample_ids, samples$sample_id)
)

samples <- samples[match(sample_ids, samples$sample_id), ]

stopifnot(identical(sample_ids, samples$sample_id))

# Create a numeric matrix with genes as row names
expression_matrix <- as.matrix(expression_data[-1])
rownames(expression_matrix) <- expression_data$gene_symbol

# Define the comparison groups
group <- factor(
  samples$disease_group,
  levels = c("healthy_control", "psoriasis")
)

# Exploratory filter: retain genes with RPKM >= 1
# in at least 20% of samples in either group
keep <- (
  rowMeans(expression_matrix[, group == "healthy_control"] >= 1) >= 0.20
) | (
  rowMeans(expression_matrix[, group == "psoriasis"] >= 1) >= 0.20
)

filtered_expression <- expression_matrix[keep, ]

# Inspect the filtering results
c(
  total_genes = nrow(expression_matrix),
  retained_genes = sum(keep),
  removed_genes = sum(!keep)
)

table(group)


# Transform RPKM values; adding 1 accommodates zeros
log_expression <- log2(filtered_expression + 1)

# Model the group difference
# Healthy control is the reference group
design <- model.matrix(~ group)

# Fit a model for each gene
fit <- limma::lmFit(log_expression, design)

# Allow residual variance to depend on average expression
fit <- limma::eBayes(fit, trend = TRUE)

# Extract results and adjust for multiple testing
de_results <- limma::topTable(
  fit,
  coef = "grouppsoriasis",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Preserve gene symbols as a column
de_results <- data.frame(
  gene_symbol = rownames(de_results),
  de_results,
  row.names = NULL
)

# Preview results
head(de_results)
nrow(de_results)



# Exploratory thresholds:
# adjusted p-value < 0.05 and absolute model effect >= 1
meets_threshold <- with(
  de_results,
  !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) >= 1
)

de_results$result_group <- "Does not meet both thresholds"

de_results$result_group[meets_threshold & de_results$logFC > 0] <-
  "Higher in psoriasis"

de_results$result_group[meets_threshold & de_results$logFC < 0] <-
  "Lower in psoriasis"

# Count genes in each category
table(de_results$result_group)

# Save all tested genes, including those below the thresholds
results_folder <- file.path(project_folder, "results/tables")
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

write.csv(
  de_results,
  file.path(results_folder, "exploratory_differential_expression.csv"),
  row.names = FALSE
)


# Connect to the existing database
con <- DBI::dbConnect(
  RSQLite::SQLite(),
  dbname = file.path(
    project_folder,
    "data/processed/psoriasis.sqlite"
  )
)

# Save model results as a database table
# Replaces this results table if you rerun the analysis
DBI::dbWriteTable(
  con,
  "differential_expression",
  de_results,
  overwrite = TRUE,
  row.names = FALSE
)

# Verify the result counts using SQL
DBI::dbGetQuery(con, "
  SELECT
    result_group,
    COUNT(*) AS gene_count
  FROM differential_expression
  GROUP BY result_group;
")

DBI::dbDisconnect(con)
