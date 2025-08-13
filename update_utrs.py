import psycopg2
import psycopg2.extras
import json
from precompute_spipv2 import get_db, log
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import md_utilities


def main():
    # script meant to update UTR variants as the paradigm has changed
    # start- and end_segment_type are only exon/intron now
    # so 5 and 3UTR variants reported before shoudl be changed as exons and
    # exon numbering should be updated accordingly
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)

    # get utr variants
    curs.execute(
        """
        SELECT *
        FROM variant a, variant_feature b
        WHERE a.feature_id = b.id
            AND (b.start_segment_type IN ('5UTR', '3UTR')
            OR b.end_segment_type IN ('5UTR', '3UTR'))
            AND a.genome_version = 'hg38'
        """
    )
    utr_vars = curs.fetchall()
    for var in utr_vars:
        log("debug", f"UTR variants: {var['c_name']} - {var['refseq']}")
        try:
            json_file = open('{0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                var['gene_symbol']
            ))
        except IOError:
            # print('file_not_found_error')
            return 'file_not_found_error'
        try:
            vv_json = json.load(json_file)
        finally:
            json_file.close()
        ncbi_chr = md_utilities.get_ncbi_chr_name(db, f"chr{var['chr']}", 'hg38')
        # log("debug", f"NCBI chr: {ncbi_chr[0]}")
        for vv_transcript in vv_json['transcripts']:
            if vv_transcript['reference'] == var['refseq']:
                if ncbi_chr[0] in vv_transcript['genomic_spans']:
                    if md_utilities.get_var_genic_csq(var['c_name']) == 'intron':
                        log("warning", f"Variant to be updated by hand: {var['id']}")
                        break
                    # we loop on exons
                    for exon in vv_transcript['genomic_spans'][ncbi_chr[0]]['exon_structure']:
                        if exon['genomic_start'] <= var['pos'] and \
                                exon['genomic_end'] >= var['pos']:
                            # log("debug", f"csq: {md_utilities.get_var_genic_csq(var['c_name'])}")
                            if (var['variant_size'] == 1) or  \
                                    (var['start_segment_type'] == var['end_segment_type'] and \
                                        var['start_segment_number'] == var['end_segment_number']):
                                # update possible
                                log("debug", f"UPDATE variant_feature SET start_segment_type = 'exon', end_segment_type = 'exon', start_segment_number = '{exon['exon_number']}', end_segment_number = '{exon['exon_number']}' WHERE id = '{var['id']}'")
                                curs.execute(
                                    """
                                    UPDATE variant_feature 
                                    SET start_segment_type = 'exon', 
                                        end_segment_type = 'exon', 
                                        start_segment_number = %s,
                                        end_segment_number = %s
                                    WHERE id = %s
                                    """,
                                    (exon['exon_number'], exon['exon_number'], var['id'])
                                )
                                break
                            else:
                                log("warning", f"Variant to be updated by hand: {var['id']}")
                                break

    db.commit()
    db_pool.putconn(db)


if __name__ == '__main__':
    main()