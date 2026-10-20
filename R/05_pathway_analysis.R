# Set the project folder
project_folder <- "."

# Download human Reactome pathway annotations
reactome_annotations <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
)

# Inspect the available fields and database version
names(reactome_annotations)
unique(reactome_annotations$db_version)

# Count pathways and annotated genes
length(unique(reactome_annotations$gs_name))
length(unique(reactome_annotations$gene_symbol))

# Preserve the downloaded annotations for reproducibility
saveRDS(
  reactome_annotations,
  file.path(
    project_folder,
    "data/processed/reactome_annotations.rds"
  )
)



# Create one row per pathway
pathways <- unique(data.frame(
  pathway_id = reactome_annotations$gs_id,
  pathway_name = reactome_annotations$gs_name,
  source = "MSigDB Reactome",
  db_version = reactome_annotations$db_version,
  stringsAsFactors = FALSE
))

# Create unique gene-to-pathway links
gene_pathways <- unique(data.frame(
  gene_symbol = reactome_annotations$gene_symbol,
  pathway_id = reactome_annotations$gs_id,
  stringsAsFactors = FALSE
))

# Load the genes tested in our model
de_results <- read.csv(
  file.path(
    project_folder,
    "results/tables/exploratory_differential_expression.csv"
  ),
  stringsAsFactors = FALSE
)

# Check exact gene-symbol matches
has_pathway <- de_results$gene_symbol %in% gene_pathways$gene_symbol

c(
  tested_genes = nrow(de_results),
  matched_to_pathways = sum(has_pathway),
  without_pathway_match = sum(!has_pathway)
)

# Retain links for tested genes only
gene_pathways <- gene_pathways[
  gene_pathways$gene_symbol %in% de_results$gene_symbol,
]

# Preview the tables
head(pathways)
head(gene_pathways)



# Connect to the existing database
con <- DBI::dbConnect(
  RSQLite::SQLite(),
  dbname = file.path(
    project_folder,
    "data/processed/psoriasis.sqlite"
  )
)

# Check identifiers and relationships before importing
stopifnot(
  !anyDuplicated(pathways$pathway_id),
  !anyDuplicated(gene_pathways),
  all(gene_pathways$pathway_id %in% pathways$pathway_id)
)

# Save both tables together
DBI::dbWithTransaction(con, {
  DBI::dbWriteTable(
    con, "pathways", pathways,
    overwrite = TRUE, row.names = FALSE
  )
  
  DBI::dbWriteTable(
    con, "gene_pathways", gene_pathways,
    overwrite = TRUE, row.names = FALSE
  )
})

# Verify the saved pathway links
DBI::dbGetQuery(con, "
  SELECT
    COUNT(*) AS gene_pathway_links,
    COUNT(DISTINCT gene_symbol) AS annotated_genes,
    COUNT(DISTINCT pathway_id) AS represented_pathways
  FROM gene_pathways;
")

DBI::dbDisconnect(con)



# Connect to the database
con <- DBI::dbConnect(
  RSQLite::SQLite(),
  dbname = file.path(project_folder, "data/processed/psoriasis.sqlite")
)

# Read and execute the saved SQL query
pathway_query <- paste(
  readLines(
    file.path(project_folder, "sql/04_pathway_analysis.sql"),
    warn = FALSE
  ),
  collapse = "\n"
)

pathway_membership <- DBI::dbGetQuery(con, pathway_query)

print(pathway_membership)

DBI::dbDisconnect(con)



# Use tested genes with a pathway annotation as the background
background_genes <- intersect(
  de_results$gene_symbol,
  gene_pathways$gene_symbol
)

# Identify threshold-meeting genes within that background
selected_genes <- intersect(
  de_results$gene_symbol[
    !is.na(de_results$adj.P.Val) &
      de_results$adj.P.Val < 0.05 &
      abs(de_results$logFC) >= 1
  ],
  background_genes
)

# Count background genes and selected genes
N <- length(background_genes)
n <- length(selected_genes)

# Group annotated genes by pathway
pathway_sets <- split(
  gene_pathways$gene_symbol,
  gene_pathways$pathway_id
)

# Exploratory size filter: 10–500 tested, annotated genes per pathway
pathway_sets <- lapply(pathway_sets, unique)
pathway_sizes <- lengths(pathway_sets)
pathway_sets <- pathway_sets[
  pathway_sizes >= 10 & pathway_sizes <= 500
]

# One-sided hypergeometric tests for overrepresentation
enrichment_results <- do.call(
  rbind,
  lapply(names(pathway_sets), function(id) {
    pathway_genes <- pathway_sets[[id]]
    K <- length(pathway_genes)
    k <- sum(pathway_genes %in% selected_genes)
    
    data.frame(
      pathway_id = id,
      background_pathway_genes = K,
      selected_pathway_genes = k,
      fold_enrichment = (k / n) / (K / N),
      p_value = phyper(
        k - 1, K, N - K, n,
        lower.tail = FALSE
      )
    )
  })
)

# Adjust across all tested pathways
enrichment_results$adjusted_p_value <- p.adjust(
  enrichment_results$p_value,
  method = "BH"
)

# Add pathway names
enrichment_results$pathway_name <- pathways$pathway_name[
  match(enrichment_results$pathway_id, pathways$pathway_id)
]

# Rank results
enrichment_results <- enrichment_results[
  order(enrichment_results$adjusted_p_value),
]

# Inspect the background and top results
c(background_genes = N, selected_genes = n)
head(enrichment_results, 10)



# Save all tested pathways
write.csv(
  enrichment_results,
  file.path(
    project_folder,
    "results/tables/pathway_enrichment_results.csv"
  ),
  row.names = FALSE
)

# Select pathways meeting the adjusted p-value threshold
significant_pathways <- subset(
  enrichment_results,
  adjusted_p_value < 0.05
)

# Count and save them
nrow(significant_pathways)

write.csv(
  significant_pathways,
  file.path(
    project_folder,
    "results/tables/significant_pathways.csv"
  ),
  row.names = FALSE
)



