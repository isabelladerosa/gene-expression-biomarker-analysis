-- Count threshold-meeting genes in each pathway
SELECT
    p.pathway_id,
    p.pathway_name,
    COUNT(DISTINCT gp.gene_symbol) AS selected_gene_count
FROM pathways AS p
JOIN gene_pathways AS gp
    ON p.pathway_id = gp.pathway_id
JOIN differential_expression AS de
    ON gp.gene_symbol = de.gene_symbol
WHERE de."adj.P.Val" < 0.05
    AND ABS(de.logFC) >= 1
GROUP BY p.pathway_id, p.pathway_name
ORDER BY selected_gene_count DESC, p.pathway_name
LIMIT 10;