# from os import listdir
# from os.path import isfile, join
# import os
import re
import sys
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import time
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config


def log(level, text):
    localtime = time.asctime( time.localtime(time.time()) )
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


def get_db():
    try:
        # read connection parameters
        params = config.mdconfig()
        db = psycopg2.connect(**params)
    except (Exception, psycopg2.DatabaseError) as error:
        log('ERROR', error)
    return db


def main():
    parser = argparse.ArgumentParser(description='Define a canonical transcript per gene when several are defined',
                                     usage='python update_canonical_when_several.py -k ncbi_api_key')
    parser.add_argument('-k', '--ncbi-api-key', default=None, required=True,
                        help='NCBI Entrez API key.')
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

    i = 0
    curs.execute(
        "select name[1] as name from gene where canonical = 't' group by name[1] having count(name[1]) > 1;"
    )
    res = curs.fetchall()
    if res is not None:
        # get canonical from NCBI
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        for hgnc_name in res:
            # print(hgnc_name)
            log('INFO', "Examining {}".format(hgnc_name['name']))
            curs.execute(
               "select * from gene where name[1] = %s",
               (hgnc_name['name'],)
            )
            full_gene = curs.fetchall()
            for acc in full_gene:
                ncbi_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={0}&api_key={1}'.format(acc['name'][1], ncbi_api_key)
                # log('DEBUG', ncbi_url)
                eutils_response = http.request('GET', ncbi_url).data.decode('utf-8')
                # if acc['name'][1] == 'NM_018257':
                # log('DEBUG', eutils_response)
                    # log('DEBUG', acc['canonical'])
                    # log('DEBUG', re.search(r'"RefSeq\sSelect\scriteria"', eutils_response))
                if re.search(r'"RefSeq\sSelect\scriteria"', eutils_response):
                    curs.execute(
                        "UPDATE gene SET canonical = 'f' WHERE name[1] = %s",
                        (acc['name'][0],)
                    )
                    log('INFO', "UPDATE gene SET canonical = 'f' WHERE name[1] = '{}'".format(acc['name'][0]))
                    curs.execute(
                        "UPDATE gene SET canonical = 't' WHERE name[2] = %s",
                        (acc['name'][1],)
                    )
                    log('INFO', "UPDATE gene SET canonical = 't' WHERE name[2] = '{}'".format(acc['name'][1]))
                    i += 1
                    log('INFO', 'Updated gene {}'.format(acc['name'][0]))
        log('INFO', 'Updated {} genes'.format(i))
    else:
        log('INFO', 'No genes to update. Everything is ok.')
    db.commit()


if __name__ == '__main__':
    main()
