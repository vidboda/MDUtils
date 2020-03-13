import sys
import re
import psycopg2
import psycopg2.extras
import argparse
import urllib3
import certifi
import json
from insert_genes import get_db
import psycopg2.extras
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config


# check UNIPROT IDs

def log(level, text):
    print()
    if level == 'ERROR':
        sys.exit('[{0}]: {1}'.format(level, text))
    print('[{0}]: {1}'.format(level, text))


def main():
    parser = argparse.ArgumentParser(description='Define a canonical transcript per gene',
                                     usage='python define_canonical.py [-r path/to/refGeneCanonical_2019_09_23.txt]')
    parser.add_argument('-k', '--ncbi-api-key', default=None, required=False,
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
        "SELECT name, np, uniprot_id, prot_size FROM gene ORDER BY name"
    )
    res = curs.fetchall()
    count = curs.rowcount
    i = 0
    for gene in res:
        # ncbi
        print('.', end="", flush=True)
        i += 1
        if i % 500 == 0:
            log('INFO', '{0}/{1} genes checked'.format(i, count))
        match_obj = re.search(r'(NP_\d+)\.\d', gene['np'])
        if match_obj:
            # log('DEBUG', gene['name'][0])
            np = match_obj.group(1)
            uniprot_url = 'https://www.ebi.ac.uk/proteins/api/proteins/refseq:{}?offset=0&size=100&reviewed=true'.format(np)
            uniprot_response = json.loads(http.request('GET', uniprot_url, headers={'Accept': 'application/json'}).data.decode('utf-8'))
            # print(uniprot_response[0]['accession'])
            try:
                if uniprot_response[0]['accession']:
                    # get uniport id prot size
                    # print('{0}-{1}'.format(gene['uniprot_id'], uniprot_response[0]['sequence']['length']))
                    if gene['uniprot_id'] == uniprot_response[0]['accession']:
                        # print('INFO: RefSeq: {0} - {1} - {2} OK'.format(gene['np'], gene['name'][1], gene['name'][0]))
                        pass
                    else:
                        curs.execute(
                            "UPDATE gene SET uniprot_id = '{0}' WHERE name[2] = '{1}'".format(uniprot_response[0]['accession'], gene['name'][1])
                        )
                        # print("UPDATE gene SET uniprot_id = '{0}' WHERE name[2] = '{1}'".format(uniprot_response[0]['accession'], gene['name'][1]))
                        log('WARNING', 'Updated gene UNIPROT ID of {0} - {1} from {2} to {3}'.format(
                            gene['name'][0],
                            gene['name'][1],
                            gene['uniprot_id'],
                            uniprot_response[0]['accession']
                        ))
                        i += 1
                else:
                    log('WARNING', 'md_uniprot_id: {0} - RefSeq: {1} - {2} - {3} :not checked'.format(
                            gene['uniprot_id'],
                            gene['np'],
                            gene['name'][1],
                            gene['name'][0]
                        ))
            except Exception:
                log('WARNING', 'no UNIPROT ID {0} for {1} - {2}'.format(uniprot_response, gene['name'][1], gene['name'][0]))
            # get prot size from eutils
            ncbi_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={0}&rettype=gp&complexity=3&api_key={1}'.format(np, ncbi_api_key)           
            # ncbi_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={0}&api_key={1}'.format(gene['name'][1], ncbi_api_key)
            prot_size = -1
            try:
                eutils_response = http.request('GET', ncbi_url).data.decode('utf-8')
                # log('DEBUG', eutils_response)
                prot_match = re.search(r'Protein\s+1\.\.(\d+)', eutils_response)#Protein\s+1\.\.(\d+)$
                if prot_match:
                    # log('DEBUG', 'ouhou')
                    prot_size = prot_match.group(1)
                    # log('DEBUG', prot_size)
            except Exception:
                log('WARNING', 'no protein size w/ got from eutils NP acc no {0}, eutils URL:{1}'.format(gene['np'], ncbi_url))
            # log('DEBUG', prot_size)
            if int(prot_size) != -1 and \
                    gene['prot_size'] is not None and \
                    int(prot_size) != int(gene['prot_size']):
                curs.execute(
                    "UPDATE gene SET prot_size = '{0}' WHERE name[2] = '{1}'".format(prot_size, gene['name'][1])
                )
                log('WARNING', 'Updated protein size for gene {0} - {1} - {2} to {3}'.format(
                    gene['name'][0],
                    gene['name'][1],
                    gene['uniprot_id'],
                    prot_size
                ))
        else:
            log('WARNING', 'pb w/ NP acc no {}'.format(gene['np']))
    log('INFO', '{} isoforms updated'.format(i))

    db.commit()


if __name__ == '__main__':
    main()
