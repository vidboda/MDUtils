import re
import sys
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import json
# import time
# from precompute_spipv2 import log
from precompute_spipv2 import get_db, log
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import md_utilities

# update genes from dev MDAPI on prod after a vv transcript update on the dev
# see update_md_transcripts.py script for details


# def log(level, text):
#     localtime = time.asctime(time.localtime(time.time()))
#     if level == 'ERROR':
#         sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
#     print('[{0}]: {1} - {2}'.format(level, localtime, text))

# def get_db():
#     try:
#         # read connection parameters
#         db_pool = psycopg2.pool.ThreadedConnectionPool(
#             1,
#             5,
#             user='mobidetails',
#             password='bio;info1',
#             host='localhost',
#             port='5433',
#             database='mobidetailsstage'
#         )
#         # db = psycopg2.connect(**params)
#     except (Exception, psycopg2.DatabaseError) as error:
#         log('ERROR', error)
#     return db_pool, db_pool.getconn()


def main():
    parser = argparse.ArgumentParser(description='Define new transcripts and optionally updates various fields',
                                     usage='python update_genes_from_remote.py [-r remote_server_url]')
    parser.add_argument('-r', '--remote-server', default='', required=True,
                        help='base URL of the remote server')
    # parser.add_argument('-k', '--ncbi-api-key', default=None, required=True,
    #                     help='NCBI Entrez API key')
    args = parser.parse_args()
    remote_addr = args.remote_server
    # if not re.search(r'\w+', args.ncbi_api_key):
    #     sys.exit('ERROR: Invalid NCBI API key, please check')
    # else:
    #     ncbi_api_key = args.ncbi_api_key
    # args = parser.parse_args(['-np'])
    print()
    log('INFO', 'Working with server {}'.format(remote_addr))

    # headers
    header = {
        'Accept': 'application/json',
        'User-Agent': 'python-requests Python/{}.{}.{}'.format(sys.version_info[0], sys.version_info[1], sys.version_info[2]),
    }
    # get db connector and cursor
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    # get local list of genes
    curs.execute(
        """
        SELECT DISTINCT(name[1]) AS gene_symbol, second_name, chr, hgnc_id, ng, strand, hgnc_name
        FROM gene
        ORDER BY name[1]
        """
    )
    genes = curs.fetchall()
    count = curs.rowcount
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    i = 0
    for gene in genes:
        i += 1
        if i % 500 == 0:
            log('INFO', '{0}/{1} genes checked'.format(i, count))
        # call MDAPI
        req_url = '{0}/api/gene/{1}'.format(remote_addr, gene['hgnc_id'])
        api_response = json.loads(http.request('GET', req_url, headers=header).data.decode('utf-8'))
        # log('DEBUG', 'req_url:{0}'.format(req_url))
        # log('DEBUG', 'api_response:{0}'.format(api_response))
        update_sql_gene = []
        for key in api_response:
            # check general features such as hgnc id, second_name, strand and update all transcript of this gene at once for these features
            api_refgene = 'NULL'
            if api_response['RefGene']:
                api_refgene = 'NG_000000.0' if api_response['RefGene'] == 'No RefGene in MobiDetails' else api_response['RefGene']
            if key == 'Chr':
                if api_response[key] != gene['chr']:
                    update_sql_gene.append(' chr = \'{0}\' '.format(api_response[key]))
            elif key == 'HGNCID':
                if api_response[key] != gene['hgnc_id']:
                    update_sql_gene.append(' hgnc_id = {0} '.format(api_response[key]))
            elif key == 'HGNCName':
                if api_response[key] != gene['hgnc_name']:
                    update_sql_gene.append(' hgnc_name = E\'{0}\' '.format(api_response[key].replace("'", "\\\'")))
            elif key == 'HGNCSymbol':
                if api_response[key] != gene['gene_symbol']:
                    update_sql_gene.append(' name[1] = \'{0}\' '.format(api_response[key]))
            elif key == 'RefGene':
                if api_refgene != gene['ng']:
                    update_sql_gene.append(' ng = \'{0}\' '.format(api_refgene))
            elif key == 'Strand':
                if api_response[key] != gene['strand']:
                    update_sql_gene.append(' strand = \'{0}\' '.format(api_response[key]))
            # transcripts
            ncbi_transcript_regexp = md_utilities.regexp['ncbi_transcript']
            nm_obj = re.search(rf'^({ncbi_transcript_regexp})$', key)
            if nm_obj:
                # check if already exists
                # if yes => update if different => check each attribute? if one is differet, update all?
                # if no create transcript
                currrent_nm = nm_obj.group(1)
                # log('DEBUG', nm_obj.group(1))
                curs.execute(
                    "SELECT * FROM gene WHERE name[2] = %s",
                    (currrent_nm,)
                )
                prod_nm = curs.fetchone()
                api_canonical = 'f' if api_response[currrent_nm]['canonical'] is False else 't'
                update_sql_transcript = []
                if prod_nm:
                    for key2 in api_response[key]:
                        # check if need to update
                        if key2 == 'refProtein':
                            if api_response[key][key2] != prod_nm['np']:
                                update_sql_transcript.append(' np = \'{0}\' '.format(api_response[key][key2]))
                        if key2 == 'UNIPROT':
                            if api_response[key][key2] != prod_nm['uniprot_id']:
                                update_sql_transcript.append(' uniprot_id = \'{0}\' '.format(api_response[key][key2]))
                        if key2 == 'proteinSize':
                            if api_response[key][key2] != prod_nm['prot_size']:
                                update_sql_transcript.append(' prot_size = {0} '.format(api_response[key][key2]))
                        if key2 == 'canonical':
                            if api_response[currrent_nm]['canonical'] != prod_nm['canonical']:
                                update_sql_transcript.append(' canonical = \'{0}\' '.format(api_canonical))
                        if key2 == 'ensemblProtein':
                            if api_response[key][key2] != prod_nm['ensp']:
                                update_sql_transcript.append(' ensp = \'{0}\' '.format(api_response[key][key2]))
                        if key2 == 'ensemblTranscript':
                            if api_response[key][key2] != prod_nm['enst']:
                                update_sql_transcript.append(' enst = \'{0}\' '.format(api_response[key][key2]))
                        if key2 == 'numberOfExons':
                            if api_response[key][key2] != prod_nm['number_of_exons']:
                                update_sql_transcript.append(' number_of_exons = {0} '.format(api_response[key][key2]))
                        if key2 == 'variantCreationTag':
                            if api_response[key][key2] != prod_nm['variant_creation']:
                                update_sql_transcript.append(' variant_creation = \'{0}\' '.format(api_response[key][key2]))
                    # update ?
                    if update_sql_transcript:
                        update_sql_transcript = "UPDATE gene SET {0} WHERE name[2] = '{1}'".format(','.join(update_sql_transcript), currrent_nm)
                        # log('DEBUG', update_sql_transcript)
                        curs.execute(update_sql_transcript)
                        db.commit()
                else:
                    # insert
                    insert_gene = """
                    INSERT INTO gene (name, second_name, chr, strand, number_of_exons, hgnc_name, prot_size, uniprot_id, ng, np, enst, ensp, canonical, variant_creation, hgnc_id )
                    VALUES ('{{\"{0}\",\"{1}\"}}', '{2}', '{3}', '{4}', {5}, '{6}', {7}, '{8}', '{9}', '{10}', '{11}', '{12}', '{13}', '{14}', {15})
                    """.format(
                        api_response['HGNCSymbol'],
                        currrent_nm,
                        gene['second_name'],
                        api_response['Chr'],
                        api_response['Strand'],
                        api_response[currrent_nm]['numberOfExons'],
                        api_response['HGNCName'].replace("'", "\\\'"),
                        api_response[currrent_nm]['proteinSize'],
                        api_response[currrent_nm]['UNIPROT'],
                        api_refgene,
                        api_response[currrent_nm]['RefProtein'],
                        api_response[currrent_nm]['ensemblTranscript'],
                        api_response[currrent_nm]['ensemblProtein'],
                        api_canonical,
                        api_response[currrent_nm]['variantCreationTag'],
                        api_response['HGNCID']
                    ).replace("None", "NULL")
                    insert_gene = insert_gene.replace("'NULL'", "NULL")
                    log('DEBUG', insert_gene)
                    db.commit()
        # update all transcripts at once
        if update_sql_gene:
            update_sql_gene = "UPDATE gene SET {0} WHERE name[1] = '{1}'".format(','.join(update_sql_gene), gene['gene_symbol'])
            # log('DEBUG', update_sql_gene)
            curs.execute(update_sql_gene)
            db.commit()
        print('.', end="", flush=True)
    db_pool.putconn(db)


if __name__ == '__main__':
    main()
