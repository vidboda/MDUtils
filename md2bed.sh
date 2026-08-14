#!/bin/bash
DATE=$(date +%Y%m%d)
FILENAME="mddev_export_vars_${DATE}.bed"

echo "Export vers ${FILENAME}..."

# echo "track type=bigBed name='MobiDetails variants' description='A bigBed File with all MobiDetails variants, date ${DATE}' bigDataUrl=https://mobidetails.chu-montpellier.fr/mdbigbed.bb url=https://mobidetails.chu-montpellier.fr/api/variant/\${id}/browser" > "${FILENAME}"

psql -p 5433 -d mobidetails -c "\copy (
    SELECT
        CONCAT('chr', a.chr) AS chrom,
        CASE
            WHEN a.g_name LIKE '%del%' THEN a.pos
            ELSE a.pos - 1
        END AS chromStart,
        CASE
            WHEN a.g_name LIKE '%dup%' THEN a.pos + LENGTH(COALESCE(a.pos_alt, '')) - 1
            WHEN a.g_name LIKE '%del%' THEN a.pos + LENGTH(COALESCE(a.pos_ref, '')) - 1
            WHEN a.g_name LIKE '%ins%' THEN a.pos
            ELSE a.pos + GREATEST(LENGTH(COALESCE(a.pos_ref, '')), LENGTH(COALESCE(a.pos_alt, ''))) - 1
        END AS chromEnd,
        CONCAT('chr', a.chr, ':g.', a.g_name) AS name,
        0,
        b.strand,
        a.feature_id,
        c.gene_symbol,
        CONCAT(c.refseq, ':c.', c.c_name)
    FROM variant a, gene b, variant_feature c
    WHERE a.feature_id = c.id
        AND b.refseq = c.refseq
        AND a.pos IS NOT NULL
        AND a.genome_version = 'hg38'
    ORDER BY chrom, chromStart
) TO STDOUT WITH (FORMAT csv, DELIMITER E'\t', HEADER false)" > "${FILENAME}"

if [ $? -eq 0 ]; then
    echo "Succès !"
else
    echo "Erreur lors de l'export."
fi
