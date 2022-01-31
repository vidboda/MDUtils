import sys
import re
import psycopg2
import psycopg2.extras
import argparse
import urllib3
import certifi
from precompute_spipv2 import get_db, log

# check UNIPROT IDs


def main():
    parser = argparse.ArgumentParser(description='Update UNIPROT ids and protein size',
                                     usage='python check_uniprot_ids.py [-k NCBI_API_KEY]')
    parser.add_argument('-k', '--ncbi-api-key', default=None, required=True,
                        help='NCBI Entrez API key. If not provided, 3rd method is not executed')
    args = parser.parse_args()
    ncbi_api_key = None
    if args.ncbi_api_key is not None:
        if not re.search(r'\w+', args.ncbi_api_key):
            sys.exit('ERROR: Invalid NCBI API key, please check')
        else:
            ncbi_api_key = args.ncbi_api_key
    # get db connector and cursor
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)

    i = 0

    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    curs.execute(
        """
        SELECT DISTINCT(np), prot_size, uniprot_id
        FROM gene
        ORDER BY np
        """
    )
    res = curs.fetchall()
    count = curs.rowcount
    i = 0
    for gene in res:
        # ncbi
        print('.', end="", flush=True)
        i += 1
        if i % 500 == 0:
            log('INFO', '{0}/{1} transcripts checked'.format(i, count))
        if gene['np']:
            # get prot size from eutils
            ncbi_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={0}&rettype=gp&complexity=3&api_key={1}'.format(gene['np'], ncbi_api_key)
            prot_size = -1
            try:
                eutils_response = http.request('GET', ncbi_url).data.decode('utf-8')
                # log('DEBUG', eutils_response)
                prot_match = re.search(r'Protein\s+1\.\.(\d+)', eutils_response)  # Protein\s+1\.\.(\d+)$
                if prot_match:
                    # log('DEBUG', 'ouhou')
                    prot_size = int(prot_match.group(1))
                    # log('DEBUG', prot_size)
            except Exception:
                log('WARNING', 'no protein size w/ eutils NP acc no {0}, eutils URL:{1}'.format(gene['np'], ncbi_url))
            # log('DEBUG', prot_size)
            if int(prot_size) != -1 and \
                    ((gene['prot_size'] is not None and
                        int(prot_size) != int(gene['prot_size'])) or
                        gene['prot_size'] is None):
                curs.execute(
                    """
                    UPDATE gene
                    SET prot_size = '{0}'
                    WHERE np = '{1}'
                    """.format(prot_size, gene['np'])
                )
                log('WARNING', 'Updated protein size for protein {0} - {1} to {2}'.format(
                    gene['np'],
                    gene['uniprot_id'],
                    prot_size
                ))
    log('INFO', '{} isoforms updated'.format(i))

    db.commit()
    db.close()


if __name__ == '__main__':
    main()
