import os
import sys
import re
import json
import psycopg2
import psycopg2.extras
import time
from precompute_spipv2 import get_db, log
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config, md_utilities


# script to set MD canonical to MANE or refseqSelect if no mane?
def update_canonical(gene_symbol, transcript, curs, db):
    curs.execute(
    #     "UPDATE gene set canonical = 'f' WHERE gene_name[1] = %s",
    #     (gene_symbol,)
    # )
    # curs.execute(
    #     "UPDATE gene set canonical = 't' WHERE gene_name[2] = %s",
    #     (transcript,)
    # )
    # db.commit()
    log('INFO', 'Changed canonical of {0} to {1}'.format(gene_symbol, transcript))

def main():
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    curs.execute(  # get each gene
        "SELECT name FROM gene WHERE canonical = 't' ORDER BY name[1]"
    )
    genes = curs.fetchall()
    for gene in genes:
        curs.execute(  # get each transcript
            "SELECT name, canonical FROM gene WHERE name[1] = %s",
            (gene['name'][1])
        )
        transcripts = curs.fetchall()
        # get refseq select, mane, etc from vv json file
        no_vv_file = 0
        try:
            json_file = open('{0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                gene['name'][0]
            ))
        except IOError:
            no_vv_file = 1
        if no_vv_file == 0:
            vv_json = json.load(json_file)
            transcript_checked = 0
            for vv_transcript in vv_json['transcripts']:
                if transcript_checked == 1:
                    break
                for md_transcript in transcripts:
                    # need to check vv isoforms against MD isoforms to keep only relevant ones
                    if vv_transcript['reference'] == md_transcript['name'][1]:
                        # MANE and canonical
                        if vv_transcript['annotations']['mane_select'] is True and \
                                md_transcript['canonical'] is True:
                            transcript_checked = 1  # ends 1st loop
                            break  # ends 2nd loop
                        if vv_transcript['annotations']['mane_select'] is True and \
                                md_transcript['canonical'] is False:
                            # update
                            update_canonical(gene['name'][0], md_transcript['name'][1], curs, db)
                            transcript_checked = 1
                            break
                        if vv_transcript['annotations']['mane_select'] is False and \
                                vv_transcript['annotations']['refseq_select'] is True and \
                                md_transcript['canonical'] is False:
                            # update
                            update_canonical(gene['name'][0], md_transcript['name'][1], curs, db)
                            break
                            # we end only the 1st loop, in case there is a MANE <> refseqSelect (should not be possible, just in case)
        else:
            log('WARNING', 'No VV file for {0}'.format(gene['name'][1]))


if __name__ == '__main__':
    main()
