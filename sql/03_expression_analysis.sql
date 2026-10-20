
-- Summarize mean expression for each gene by group
-- Preview average expression by gene and disease group
SELECT
    e.gene_symbol,
    s.disease_group,
    COUNT(*) AS sample_count,
    AVG(e.rpkm) AS mean_rpkm
FROM expression AS e
JOIN samples AS s
    ON e.sample_id = s.sample_id
GROUP BY e.gene_symbol, s.disease_group
ORDER BY e.gene_symbol, s.disease_group
LIMIT 10;

-- Summarize both group averages on one row per gene
-- These are descriptive averages, not significance tests
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