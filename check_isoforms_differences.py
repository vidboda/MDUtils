import re
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
# import datetime
import json
import pprint
import subprocess
from precompute_spipv2 import get_db, log


def main():
    parser = argparse.ArgumentParser(description='Check isoforms differences between 2 DBs',
                                     usage='python check_isoforms_differences.py -k md_api_key -nk ncbi_api_key')
    parser.add_argument('-k', '--api-key', default=None, required=True,
                        help='Your API key visible on your profile page on the website.')
    parser.add_argument('-nk', '--ncbi-api-key', default=None, required=True,
                        help='NCBI Entrez API key.')
    parser.add_argument('-md', '--diff-md', default='', required=False,
                        help='Check differences between MD Dev and Prod', action='store_true')
    # parser.add_argument('-ncbi', '--diff-ncbi', default='', required=False,
    #                     help='Check differences between NCBI RefSeq and MD', action='store_true')
    # ncbi_api_key = None
    args = parser.parse_args()
    # if args.ncbi_api_key is not None:
    #     if not re.search(r'\w+', args.ncbi_api_key):
    #         log('ERROR', 'Invalid NCBI API key, please check')
    #     else:
    #         ncbi_api_key = args.ncbi_api_key
    # get db connector and cursor
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    username = None
    if len(args.api_key) != 43:
        log('ERROR', 'Invalid API key, please check it')
    else:
        api_key = args.api_key
        # user
        curs.execute(
            """
            SELECT username
            FROM mobiuser
            WHERE api_key = %s
            """,
            (api_key,)
        )
        res_user = curs.fetchone()
        if res_user is None:
            log('ERROR', 'Unknown API key')
        username = res_user['username']
        log('INFO', 'User: {}'.format(username))

    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    if args.diff_md:
        # meant to be ran on prod server
        diff = {}
        # get genes w/ more than one isoform
        curs.execute(
            """
            SELECT gene_symbol, refseq
            FROM gene
            WHERE canonical = 't'
                AND gene_symbol IN (
                    SELECT gene_symbol
                    FROM gene
                    GROUP BY gene_symbol
                    HAVING COUNT(gene_symbol) > 1
                )
            ORDER by gene_symbol
            """
        )
        res = curs.fetchall()
        i = 0
        for gene in res:
            i += 1
            log('INFO', 'Treating gene {0} - #{1}'.format(gene['gene_symbol'], i))
            full_nm = gene['refseq']
            base_url = "http://10.34.20.79"
            md_url = '{0}/MD/api/gene/{1}'.format(base_url, gene['gene_symbol'])
            http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
            try:
                md_data = json.loads(http.request('GET', md_url).data.decode('utf-8'))
            except Exception:
                log('WARNING', 'MD not responding for {}'.format(gene['hgnc']))
            if full_nm in md_data and \
                    md_data[full_nm]['canonical'] is True:
                # log('INFO', 'No change for {}'.format(gene['hgnc']))
                continue
            elif full_nm not in md_data:
                log('DEBUG', '{0}-{1}'.format(md_data, full_nm))
            for key in md_data:
                matchobj = re.search(r'^(NM_\d+)\.\d+$', key)
                if matchobj:
                    new_nm = matchobj.group(1)
                    if md_data[key]['canonical'] is True:
                        diff[gene['gene_symbol']] = {
                            'old_can': full_nm,
                            'new_can': key
                        }
                        log('INFO', 'updating canonical for {0}: {1} instead of {2}'.format(gene['gene_symbol'], key, full_nm))
                        curs.execute(
                            """
                            UPDATE gene
                            SET canonical = 'f'
                            WHERE gene_symbol = %s
                            """,
                            (gene['gene_symbol'],)
                        )
                        curs.execute(
                            """
                            UPDATE gene
                            SET canonical = 't'
                            WHERE refseq = %s
                            """,
                            (new_nm,)
                        )
                        db.commit()
                        cmd = "python3 update_vars_when_iso_change.py -k {0} -g {1}".format(api_key, gene['gene_symbol'])
                        returned_value = subprocess.call(cmd, shell=True)
                        log('INFO', 'Variants update returned value for {0}: {1}'.format(gene['gene_symbol'], returned_value))
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(diff)
    db.close()


if __name__ == '__main__':
    main()
