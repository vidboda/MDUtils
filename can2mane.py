import json
import psycopg2
import psycopg2.extras
from precompute_spipv2 import get_db, log
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import md_utilities


# script to set MD canonical to MANE or refseqSelect if no mane?
def update_canonical(gene_symbol, transcript, curs, db):
    curs.execute(
        """
        UPDATE gene
        SET canonical = 'f'
        WHERE name[1] = %s
        """,
        (gene_symbol,)
    )
    curs.execute(
        """
        UPDATE gene
        SET canonical = 't'
        WHERE name[2] = %s
        """,
        (transcript,)
    )
    db.commit()
    log('INFO', 'Changed canonical of {0} to {1}'.format(gene_symbol, transcript))


def main():
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    curs.execute(  # get each gene
        """
        SELECT DISTINCT(name[1]) AS gene_symbol
        FROM gene
        ORDER BY name[1]
        """
    )
    genes = curs.fetchall()
    for gene in genes:
        curs.execute(  # get each transcript
            "SELECT name, canonical FROM gene WHERE name[1] = %s",
            (gene['gene_symbol'],)
        )
        transcripts = curs.fetchall()
        # get refseq select, mane, etc from vv json file
        no_vv_file = 0
        try:
            json_file = open('{0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                gene['gene_symbol']
            ))
        except IOError:
            no_vv_file = 1
        if no_vv_file == 0:
            vv_json = json.load(json_file)
            transcript_checked = 0
            if 'transcripts' in vv_json:
                for vv_transcript in vv_json['transcripts']:
                    if transcript_checked == 1:
                        break
                    for md_transcript in transcripts:
                        # need to check vv isoforms against MD isoforms to keep only relevant ones
                        if vv_transcript['reference'] == md_transcript['name'][1]:
                            # log('DEBUG', 'VV:{0}-MD:{1}'.format(vv_transcript['reference'], md_transcript['name'][1]))
                            # MANE and canonical
                            if 'mane_select' in vv_transcript['annotations'] and \
                                    vv_transcript['annotations']['mane_select'] is True and \
                                    md_transcript['canonical'] is True:
                                transcript_checked = 1  # ends 1st loop
                                break  # ends 2nd loop
                            if 'mane_select' in vv_transcript['annotations'] and \
                                    vv_transcript['annotations']['mane_select'] is True and \
                                    md_transcript['canonical'] is False:
                                # update
                                update_canonical(gene['gene_symbol'], md_transcript['name'][1], curs, db)
                                transcript_checked = 1
                                break
                            if 'mane_select' not in vv_transcript['annotations']:
                                log('WARNING', 'No MANE in VV file for {0}'.format(gene['gene_symbol']))
                            # if vv_transcript['annotations']['mane_select'] is False and \
                            #         vv_transcript['annotations']['refseq_select'] is True and \
                            #         md_transcript['canonical'] is False:
                            #     # update
                            #     update_canonical(gene['gene_symbol'], md_transcript['name'][1], curs, db)
                            #     break
                                # we end only the 1st loop, in case there is a MANE <> refseqSelect (should not be possible, just in case)
            else:
                log('WARNING', 'No transcript in VV file for {0}'.format(gene['gene_symbol']))
        else:
            log('WARNING', 'No VV file for {0}'.format(gene['gene_symbol']))
    db.close()


if __name__ == '__main__':
    main()
