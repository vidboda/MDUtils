import sys
import os
import re
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import time
from insert_genes import get_db
import psycopg2.extras
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config

# first genes in database did not have canonical transcripts.
# fixed with the refGenecanonical file


def log(level, text):
    localtime = time.asctime( time.localtime(time.time()) )
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


def main():
    parser = argparse.ArgumentParser(description='Define a canonical transcript per gene',
                                     usage='python define_canonical.py [-r path/to/refGeneCanonical_2019_09_23.txt]')
    parser.add_argument('-r', '--refgene', default='', required=True,
                        help='Path to the file containing the canonical refSeq IDs per gene (from UCSC)')
    parser.add_argument('-k', '--ncbi-api-key', default=None, required=False,
                        help='NCBI Entrez API key. If not provided, 3rd method is not executed')
    parser.add_argument('-u', '--update-refgene', default=None, required=False,
                        help='Update RefGene (canonical) for genes w/ on variants based on NCBI (requires NCBI API key)', action='store_true')
    args = parser.parse_args()
    # get file
    if os.path.isfile(args.refgene):
        refgeneFile = args.refgene
    else:
        sys.exit('ERROR: Invalid input path, please check your command')
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
    # lacking_nm = []

    # first when only one isoform => canonical

    curs.execute(
        "SELECT name FROM gene WHERE canonical = 'f' AND name[1] IN (SELECT name[1] FROM gene GROUP BY name[1] HAVING COUNT (name[1]) = 1)"
    )
    res = curs.fetchall()
    for acc in res:
        curs.execute(
            "UPDATE gene SET canonical = 't' WHERE name[2] = '{}'".format(acc['name'][1])
        )
        # lacking_nm.append(acc['name'][0])
        log('INFO', 'Updated gene {} (1st method)'.format(acc['name'][0]))
        i += 1
    db.commit()
    # second check the refgene file

    for geneLine in open(refgeneFile).readlines():
        # ENST - NM - gene
        geneLineList = geneLine.rstrip().split("\t")
        # print(geneLineList[2])
        if geneLineList[2] != 'n/a' and geneLineList[2] != 'hg38.refGene.name2':
            curs.execute(  # gene exists in MD (no main already set)
                "SELECT DISTINCT(name[1]) FROM gene WHERE name[1] = '{}' AND name[1] NOT IN \
                                                                    (SELECT name[1] FROM gene WHERE canonical = 't') ORDER BY name".format(geneLineList[2])
            )
            mdgene = curs.fetchone()

            if mdgene is not None:
                # nm exists in md?
                curs.execute(
                    "SELECT name FROM gene WHERE name[2] = '{0}'".format(geneLineList[1])
                )  # exists in table gene_annotation? get a nm
                mdnm = curs.fetchone()
                if mdnm is not None:
                    # ok => canonical
                    i += 1
                    postGene = '{"' + mdnm['name'][0] + '","' + mdnm['name'][1] + '"}'
                    # print("UPDATE gene SET canonical = 't' WHERE name = '{}'".format(postGene))
                    curs.execute(
                         "UPDATE gene SET canonical = 't' WHERE name = '{}'".format(postGene)
                    )
                    log('INFO', 'Updated gene {} (2nd method)'.format(mdnm['name'][0]))
                # else:
                    # lacking_nm.append(geneLineList[2])
    # print(lacking_nm)
    db.commit()
    # 3rd get info at NCBI
    # API key mandatory
    if ncbi_api_key is not None:
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        # get list of remaining genes with no canonical defined
        curs.execute(
            "SELECT name, np, canonical FROM gene WHERE name[1] NOT IN (SELECT name[1] FROM gene WHERE canonical='t') ORDER BY name"
        )
        res = curs.fetchall()
        for acc in res:
            # ncbi
            ncbi_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={0}&api_key={1}'.format(acc['name'][1], ncbi_api_key)
            eutils_response = http.request('GET', ncbi_url).data.decode('utf-8')
            if re.search(r'"RefSeq\sSelect"', eutils_response):
                curs.execute(
                    "UPDATE gene SET canonical = 't' WHERE name[2] = '{}'".format(acc['name'][1])
                )
                i += 1
                log('INFO', 'Updated gene {} (3rd method)'.format(acc['name'][0]))
            if acc['np'] == 'NP_000000.0':
                if re.search(r'accession\s"NP_\d+",\s+version\s\d$', eutils_response, re.MULTILINE):
                    match_object = re.search(r'accession\s"(NP_\d+)",\s+version\s(\d+)$', eutils_response, re.MULTILINE)
                    curs.execute(
                        "UPDATE gene SET np = '{0}.{1}' WHERE name[2] = '{2}'".format(match_object.group(1), match_object.group(2), acc['name'][1])
                    )
                    log('INFO', 'Updated gene NP acc no of {0} to {1}.{2}'.format(acc['name'][0], match_object.group(1), match_object.group(2)))
        if args.update_refgene:
            # get genes w/ no variants, and at least 2 isoforms to check which one should be canonical
            curs.execute(
                "SELECT name, canonical FROM gene WHERE (name[1] NOT IN \
                (SELECT gene_name[1] FROM variant_feature)) AND \
                (name[1] IN (SELECT name[1] FROM gene GROUP BY name[1] \
                HAVING COUNT (name[1]) > 1)) ORDER BY name"
            )
            res = curs.fetchall()
            for acc in res:
                ncbi_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={0}&api_key={1}'.format(
                    acc['name'][1],
                    ncbi_api_key
                )
                # log('DEBUG', ncbi_url)
                eutils_response = http.request('GET', ncbi_url).data.decode('utf-8')
                # if acc['name'][1] == 'NM_018257':
                    # log('DEBUG', eutils_response)
                    # log('DEBUG', acc['canonical'])
                    # log('DEBUG', re.search(r'"RefSeq\sSelect\scriteria"', eutils_response))
                if re.search(r'"RefSeq\sSelect\scriteria"', eutils_response) and acc['canonical'] is False:
                    curs.execute(
                        "UPDATE gene SET canonical = 'f' WHERE name[1] = '{}'".format(acc['name'][0])
                    )
                    # log('INFO', "UPDATE gene SET canonical = 'f' WHERE name[1] = '{}'".format(acc['name'][0]))
                    curs.execute(
                        "UPDATE gene SET canonical = 't' WHERE name[2] = '{}'".format(acc['name'][1])
                    )
                    # log('INFO', "UPDATE gene SET canonical = 't' WHERE name[2] = '{}'".format(acc['name'][1]))
                    i += 1
                    log('INFO', 'Updated gene {} (4th method)'.format(acc['name'][0]))
    log('INFO', '{} genes modified'.format(i))

    db.commit()


if __name__ == '__main__':
    main()
