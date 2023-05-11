import os
import json
import psycopg2
import psycopg2.extras
from precompute_spipv2 import get_db, log
from MobiDetailsApp import md_utilities

def main():
    # script to get a list of gene symbols for which the MD canonical <> MANE Select
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)

    curs.execute(  # get genes - one transcript per gene (canonical)
        """
        SELECT gene_symbol, refseq
        FROM gene
        WHERE canonical = 't'
        ORDER by gene_symbol
        """
    )
    res_gene = curs.fetchall()
    for gene in res_gene:
        vv_file = open('{0}{1}.json'.format(
            md_utilities.local_files['variant_validator']['abs_path'],
            gene['gene_symbol']
        ), 'rb')
        vv_json = json.load(vv_file)
        if 'error' in vv_json or \
                ('message' in vv_json and
                vv_json['message'] == 'Internal Server Error'):
            continue
        update_required = any(
            vv_transcript['reference'] == gene['refseq']
            and 'mane_select' in vv_transcript['annotations']
            and vv_transcript['annotations']['mane_select'] is False
            for vv_transcript in vv_json['transcripts']
        )
        if update_required:
            # loop again on all transcripts
            for vv_transcript in vv_json['transcripts']:
                if 'mane_select' in vv_transcript['annotations'] and \
                        vv_transcript['annotations']['mane_select']:
                    # check if variant_creation is ok for the new MANE
                    curs.execute(
                        """
                        SELECT variant_creation
                        FROM gene
                        WHERE refseq = %s
                        """,
                        (vv_transcript['reference'],)
                    )
                    res_create = curs.fetchone()
                    if res_create and res_create['variant_creation'] in [
                        'ok',
                        'hg19_mapping_default',
                    ]:
                        noupdate = any(
                            'mane_plus_clinical'
                            in vv_transcript2['annotations']
                            and vv_transcript2['annotations'][
                                'mane_plus_clinical'
                            ]
                            and vv_transcript2['reference'] == gene['refseq']
                            for vv_transcript2 in vv_json['transcripts']
                        )
                        if not noupdate:
                            log('INFO', 'NonMANECanonical:{0} - MANE: {1}'.format(gene['gene_symbol'], vv_transcript['reference']))
                            curs.execute(
                                """
                                UPDATE gene
                                SET canonical = 'f'
                                WHERE refseq = %s
                                """,
                                (gene['refseq'],)
                            )
                            curs.execute(
                                """
                                UPDATE gene
                                SET canonical = 't'
                                WHERE refseq = %s
                                """,
                                (vv_transcript['reference'],)
                            )
                            db.commit()
                        else:
                            log('INFO', 'NonMANECanonicalNoUpdatePlusClinical:{0} - MANE: {1}'.format(gene['gene_symbol'], vv_transcript['reference']))
                    else:
                        log('WARNING', 'NonMANENoCreationCanonical:{0} - MANE: {1}'.format(gene['gene_symbol'], vv_transcript['reference']))


    db_pool.putconn(db)

if __name__ == '__main__':
    main()
