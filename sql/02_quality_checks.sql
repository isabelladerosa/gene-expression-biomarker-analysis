
-- Count samples in each disease group
SELECT
    disease_group,
    COUNT(*) AS sample_count
FROM samples
GROUP BY disease_group;

-- Check how many expression records each sample has
SELECT
    MIN(gene_count) AS minimum_genes_per_sample,
    MAX(gene_count) AS maximum_genes_per_sample
FROM (
    SELECT sample_id, COUNT(*) AS gene_count
    FROM expression
    GROUP BY sample_id
);

-- Check for broken foreign-key relationships
PRAGMA foreign_key_check;





