import re
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import urllib
import certifi
import json
from precompute_spipv2 import get_db, log
from MobiDetailsApp import md_utilities

# removes all c.2del then checks that these variants can be created in genes using API
# to be used on private dev server


def main():
    parser = argparse.ArgumentParser(description='Checks that genes accept variant creation', usage='python check_variant_creation.py [-r remote_server_url]')
    parser.add_argument('-r', '--remote-server', default='', required=True, help='base URL of the remote server')
    parser.add_argument('-k', '--api-key', default='', required=True, help='Your API key visible on your profile page on the website.')

    args = parser.parse_args()
    remote_addr = args.remote_server
    if re.search(r'mobidetails\.iurc', remote_addr):
        log('ERROR', 'This script is not intended to work with the production server')
    if len(args.api_key) != 43:
        log('ERROR', 'Invalid API key, please check it')
    else:
        api_key = args.api_key
    print()
    log('INFO', f'Working with server {remote_addr}')

    # get db connector and cursor
    # db = get_db()
    # curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    # curs.execute(
    #    """
    #    DELETE FROM variant_feature
    #    WHERE c_name = '2del'
    #    """
    # )
    # db.commit()
    # # reinitialise gene state
    # curs.execute(
    #    """
    #    UPDATE gene
    #    SET variant_creation = 'ok'
    #    """
    # )
    # db.commit()

    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())

    curs.execute(
        """
        SELECT gene_symbol, refseq, variant_creation
        FROM gene
        WHERE variant_creation = 'ok'
        ORDER BY name
        """
    )
    #  AND variant_creation IN ('hg19_mapping_default', 'hg38_mapping_default')
    transcripts = curs.fetchall()
    num_transcripts = curs.rowcount
    i = 0
    j = 0
    k = 0
    # variant = 'c.1A>T'
    failed_genes = []
    for transcript in transcripts:
        print('.', end="", flush=True)
        i += 1
        if i % 500 == 0:
            log('INFO', '{0}/{1} transcripts checked'.format(i, num_transcripts))
        # variant = '{0}.{1}:c.1A>T'.format(gene['name'][1], gene['nm_version'])
        # md_url = '{0}/api/variant/create/{1}/{2}'.format(remote_addr, variant, api_key)
        md_url = '{0}/api/variant/create'.format(remote_addr)
        variant_chgvs = '{0}:c.2del'.format(transcript['refseq'])
        data = {
            'variant_chgvs': urllib.parse.quote(variant_chgvs),
            'caller': 'cli',
            'api_key': api_key
        }
        # reinitialise gene state before query
        curs.execute(
            """
            UPDATE gene
            SET variant_creation = 'ok'
            WHERE refseq = %s
            """,
            (transcript['refseq'],)
        )
        db.commit()
        try:
            md_response = json.loads(http.request('POST', md_url, headers=md_utilities.api_agent, fields=data).data.decode('utf-8'))
        except Exception:
            log('WARNING', f'failed MD API call for {variant_chgvs}')
            # put back transcript state as it was
            curs.execute(
                """
                UPDATE gene
                SET variant_creation = %s
                WHERE refseq = %s
                """,
                (transcript['variant_creation'], transcript['refseq'])
            )
            db.commit()
            k += 1
            failed_genes.append(f"{transcript['refseq']}")
            continue
        if 'mobidetails_error' in md_response:
            j += 1
            log('WARNING', 'variant creation failed for transcript {0} with error {1}'.format(transcript['refseq'], md_response['mobidetails_error']))
            # put back transcript state as it was
            variant_creation_status = transcript['variant_creation']
            if variant_creation_status == 'ok':
                variant_creation_status = 'mapping_default'
            if re.search(r'Internal\sServer\sError', md_response['mobidetails_error']):
                variant_creation_status = 'vv_server_error'
            curs.execute(
                """
                UPDATE gene
                SET variant_creation = %s
                WHERE refseq = %s
                """,
                (variant_creation_status, transcript['refseq'])
            )
            db.commit()
        elif 'mobidetails_id' in md_response and transcript['variant_creation'] != 'ok':
            curs.execute(
                """
                UPDATE gene
                SET variant_creation = 'ok'
                WHERE refseq = '{}'
                """.format(transcript['refseq'])
            )
            db.commit()
    db_pool.putconn(db)
    log('INFO', '{0}/{1} transcripts reported a VV error'.format(j, num_transcripts))
    log('INFO', '{0}/{1} transcripts triggered an MD error:'.format(k, num_transcripts))
    log('INFO', failed_genes)


if __name__ == '__main__':
    main()
