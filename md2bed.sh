#!/bin/bash
DATE=$(date +%Y%m%d)
FILENAME="/home/mobidetails/mddev_export_vars_${DATE}.bed"

echo "Export vers $FILENAME..."

psql -p 5433 -d mobidetails -c "\copy (
    SELECT
        CONCAT('chr', chr) AS chrom,
        CASE
            WHEN g_name LIKE '%del%' THEN pos
            ELSE pos - 1
        END AS chromStart,
        CASE
            WHEN g_name LIKE '%ins%' THEN pos
            WHEN g_name LIKE '%dup%' THEN pos + LENGTH(COALESCE(pos_alt, '')) - 1
            WHEN g_name LIKE '%del%' THEN pos + LENGTH(COALESCE(pos_ref, '')) - 1
            ELSE pos + GREATEST(LENGTH(COALESCE(pos_ref, '')), LENGTH(COALESCE(pos_alt, ''))) - 1
        END AS chromEnd,
        CONCAT('chr', chr, ':g.', g_name) AS name,
        0, '+',
        CONCAT('https://mobidetails.chu-montpellier.fr/api/variant/', feature_id, '/browser')
    FROM variant
    WHERE pos IS NOT NULL
        AND genome_version = 'hg38'
    ORDER BY chrom, chromStart
) TO '$FILENAME' WITH (FORMAT csv, DELIMITER E'\t', HEADER false)"

if [ $? -eq 0 ]; then
    echo "Succès !"
else
    echo "Erreur lors de l'export."
fi
