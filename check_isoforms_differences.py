import re
import sys
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import time
import datetime
import json
from insert_genes import get_db
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config


def log(level, text):
    localtime = time.asctime( time.localtime(time.time()) )
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))

def main():
    # script meant to update variants when the canonical form of a gene is changed
    parser = argparse.ArgumentParser(description='Define a canonical transcript per gene when several are defined',
                                     usage='python check_isoforms_differences.py -k md_api_key -nk ncbi_api_key')
    parser.add_argument('-k', '--api-key', default=None, required=True,
                        help='Your API key visible on your profile page on the website.')
    parser.add_argument('-nk', '--ncbi-api-key', default=None, required=True,
                        help='NCBI Entrez API key.')
    parser.add_argument('-md', '--diff-md', default='', required=False,
                        help='Check differences between MD Dev and Prod', action='store_true')
    parser.add_argument('-ncbi', '--diff-ncbi', default='', required=False,
                        help='Check differences between NCBI RefSeq and MD', action='store_true')
    ncbi_api_key = None
    args = parser.parse_args()
    if args.ncbi_api_key is not None:
        if not re.search(r'\w+', args.ncbi_api_key):
            log('ERROR: Invalid NCBI API key, please check')
        else:
            ncbi_api_key = args.ncbi_api_key
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
            "SELECT username FROM mobiuser WHERE api_key = %s",
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
            "SELECT name[1] AS hgnc, name[2] AS nm, nm_version FROM gene WHERE canonical = 't' AND name[1] IN (SELECT name[1] FROM gene GROUP BY name[1] HAVING COUNT(name[1]) > 1) ORDER by name[1]"
        )
        res = curs.fetchall()
        i = 0
        for gene in res:
            i += 1
            log('INFO', 'Treating gene {0} - #{1}'.format(gene['hgnc'], i))
            full_nm = '{0}.{1}'.format(gene['nm'], gene['nm_version'])
            base_url = "http://10.34.20.79"
            md_url = '{0}/MD/api/gene/{1}'.format(base_url, gene['hgnc'])
            http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
            try:
                md_data = json.loads(http.request('GET', md_url).data.decode('utf-8'))
            except Exception:
                log('WARNING', 'MD not responding for {}'.format(gene['hgnc']))
            if full_nm in md_data and \
                    md_data[full_nm]['canonical'] is True:
                # log('INFO', 'No change for {}'.format(gene['hgnc']))
                continue
            else:
                log('DEBUG', '{0}-{1}'.format(md_data, full_nm))
            for key in md_data:
                if re.search('NM_.+', key):
                    if md_data[key]['canonical'] is True:
                        log('INFO', 'New canonical for {0}: {1} instead of {2}'.format(gene['hgnc'], key, full_nm))
                        diff[gene['hgnc']] = {
                            'old_can': full_nm,
                            'new_can': key
                        }
        print(diff)

if __name__ == '__main__':
    main()