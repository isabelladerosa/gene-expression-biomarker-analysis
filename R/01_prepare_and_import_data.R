# Import psoriasis gene-expression data
expression_data <- read.delim(
  gzfile("data/raw/GSE54456_RPKM_samples.txt"),
  check.names = FALSE
)

# Name the first column
names(expression_data)[1] <- "gene_symbol"

# Inspect dimensions and column names
dim(expression_data)
head(names(expression_data))

# Check for missing or blank gene symbols
sum(is.na(expression_data$gene_symbol) |
      trimws(expression_data$gene_symbol) == "")

# Check for duplicate gene symbols
sum(duplicated(expression_data$gene_symbol))

# Download sample metadata from GEO
geo_data <- GEOquery::getGEO(
  "GSE54456",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

# Extract the first platform's sample metadata
sample_metadata <- Biobase::pData(geo_data[[1]])

# Inspect its size and available columns
dim(sample_metadata)
names(sample_metadata)

# Inspect sample identifiers and tissue labels
head(sample_metadata[, c(
  "title",
  "geo_accession",
  "source_name_ch1",
  "tissue type:ch1"
)])

# Count samples in each tissue category
table(sample_metadata[["tissue type:ch1"]], useNA = "ifany")


# Get sample IDs from the expression table
expression_sample_ids <- names(expression_data)[-1]

# Check for duplicate metadata IDs
sum(duplicated(sample_metadata$title))

# Find expression samples missing from the metadata
setdiff(expression_sample_ids, sample_metadata$title)

# Find metadata samples missing from the expression table
setdiff(sample_metadata$title, expression_sample_ids)

# Align metadata rows with expression columns
sample_metadata <- sample_metadata[
  match(expression_sample_ids, sample_metadata$title),
]

# Confirm that the order matches exactly
identical(expression_sample_ids, sample_metadata$title)

# Create a clean sample table
samples <- data.frame(
  sample_id = sample_metadata$title,
  geo_accession = sample_metadata$geo_accession,
  tissue_type = sample_metadata[["tissue type:ch1"]],
  stringsAsFactors = FALSE
)

# Assign readable group labels
samples$disease_group <- ifelse(
  samples$tissue_type == "lesional psoriatic skin",
  "psoriasis",
  "healthy_control"
)

# Preview and check group counts
head(samples)
table(samples$disease_group)

# Save the table to your existing folder
write.csv(
  samples,
  "data/processed/sample_metadata.csv",
  row.names = FALSE
)

# Check that every sample column is numeric
all(vapply(expression_data[-1], is.numeric, logical(1)))

# Count missing expression values
sum(is.na(expression_data[-1]))

# Find the smallest and largest expression values
range(as.matrix(expression_data[-1]), na.rm = TRUE)

# Create a table containing one row per gene
genes <- data.frame(
  gene_symbol = expression_data$gene_symbol,
  stringsAsFactors = FALSE
)

# Check the number of genes and preview the table
nrow(genes)
head(genes)

# Save the gene table
write.csv(
  genes,
  "data/processed/genes.csv",
  row.names = FALSE
)

# Reshape expression values into a long table
expression_long <- data.frame(
  gene_symbol = rep(
    expression_data$gene_symbol,
    times = ncol(expression_data) - 1
  ),
  sample_id = rep(
    names(expression_data)[-1],
    each = nrow(expression_data)
  ),
  rpkm = as.vector(as.matrix(expression_data[-1])),
  stringsAsFactors = FALSE
)

# Preview the table and check its size
head(expression_long)
nrow(expression_long)

# Confirm that all genes and samples match our reference tables
all(expression_long$gene_symbol %in% genes$gene_symbol)
all(expression_long$sample_id %in% samples$sample_id)


# Connect to a SQLite database saved in your processed folder
con <- DBI::dbConnect(
  RSQLite::SQLite(),
  dbname = "data/processed/psoriasis.sqlite"
)

# Enforce relationships between tables
DBI::dbExecute(con, "PRAGMA foreign_keys = ON;")

# Create the genes table
DBI::dbExecute(con, "
  CREATE TABLE IF NOT EXISTS genes (
    gene_symbol TEXT PRIMARY KEY NOT NULL
  );
")

# Create the samples table
DBI::dbExecute(con, "
  CREATE TABLE IF NOT EXISTS samples (
    sample_id TEXT PRIMARY KEY NOT NULL,
    geo_accession TEXT NOT NULL UNIQUE,
    tissue_type TEXT NOT NULL,
    disease_group TEXT NOT NULL
      CHECK (disease_group IN ('psoriasis', 'healthy_control'))
  );
")

# Create the expression table
DBI::dbExecute(con, "
  CREATE TABLE IF NOT EXISTS expression (
    gene_symbol TEXT NOT NULL,
    sample_id TEXT NOT NULL,
    rpkm REAL NOT NULL CHECK (rpkm >= 0),
    PRIMARY KEY (gene_symbol, sample_id),
    FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
  );
")

# Confirm the tables exist
DBI::dbListTables(con)

# Import all three tables together
# If an import fails, undo this block's database changes
# Import only when all three tables are empty
existing_counts <- c(
  genes = DBI::dbGetQuery(
    con, "SELECT COUNT(*) AS n FROM genes;"
  )$n,
  samples = DBI::dbGetQuery(
    con, "SELECT COUNT(*) AS n FROM samples;"
  )$n,
  expression = DBI::dbGetQuery(
    con, "SELECT COUNT(*) AS n FROM expression;"
  )$n
)

if (all(existing_counts == 0)) {
  
  DBI::dbWithTransaction(con, {
    DBI::dbAppendTable(con, "genes", genes)
    DBI::dbAppendTable(con, "samples", samples)
    DBI::dbAppendTable(con, "expression", expression_long)
  })
  
  message("Data imported successfully.")
  
} else {
  stop(
    "Database already contains data. Import skipped to prevent duplicates. ",
    "Use the existing database for analysis."
  )
}

# Check the number of rows in each database table
DBI::dbGetQuery(con, "
  SELECT 'genes' AS table_name, COUNT(*) AS row_count FROM genes
  UNION ALL
  SELECT 'samples', COUNT(*) FROM samples
  UNION ALL
  SELECT 'expression', COUNT(*) FROM expression;
")



