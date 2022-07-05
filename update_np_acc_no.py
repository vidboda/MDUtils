import sys
import re
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
from insert_genes import get_db


# update genes for which np acc no is NP_000000.0
# Still useful?????????


def log(level, text):
    if level == 'ERROR':
        sys.exit('[{0}]: {1}'.format(level, text))
    print('[{0}]: {1}'.format(level, text))


def main():
    parser = argparse.ArgumentParser(description='Defines NP RefSeq acc_no when lacking', usage='python update_np_acc_no.py -k ncbi_api_key')
    parser.add_argument('-k', '--ncbi-api-key', default=None, required=True, help='NCBI Entrez API key.')
    args = parser.parse_args()
    # get file

    ncbi_api_key = None
    if args.ncbi_api_key is not None:
        if not re.search(r'\w+', args.ncbi_api_key):
            log('ERROR', 'Invalid NCBI API key, please check')
        else:
            ncbi_api_key = args.ncbi_api_key

    # get db connector and cursor
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)

    i = 0

    if ncbi_api_key is not None:
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        # get list of remaining genes with no canonical defined
        curs.execute(
            """
            SELECT gene_symbol, refseq, np
            FROM gene
            WHERE np = 'NP_000000.0'
            ORDER by name
            """
        )
        res = curs.fetchall()
        for acc in res:
            # ncbi
            ncbi_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={0}&api_key={1}'.format(acc['refseq'], ncbi_api_key)
            eutils_response = http.request('GET', ncbi_url).data.decode('utf-8')
            match_object = re.search(r'accession\s"(NP_\d+)",\s+version\s(\d+)$', eutils_response, re.MULTILINE)
            if match_object:
                curs.execute(
                    """
                    UPDATE gene
                    SET np = '{0}.{1}'
                    WHERE refseq = '{2}'
                    """.format(match_object.group(1), match_object.group(2), acc['refseq'])
                )
                log('INFO', 'Updated gene NP acc no of {0} to {1}.{2}'.format(acc['gene_symbol'], match_object.group(1), match_object.group(2)))
                i += 1
    log('INFO', '{} genes updated'.format(i))

    db.commit()
    db.close()


if __name__ == '__main__':
    main()
