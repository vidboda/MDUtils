import re
import sys
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import json
# import time
from precompute_spipv2.py import get_db, log
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import md_utilities

# update genes from dev MDAPI on prod after a vv transcript update on the dev
# see update_md_transcripts.py script for details


# def log(level, text):
#     localtime = time.asctime(time.localtime(time.time()))
#     if level == 'ERROR':
#         sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
#     print('[{0}]: {1} - {2}'.format(level, localtime, text))


def main():
    parser = argparse.ArgumentParser(description='Define new transcripts and optionally updates various fields',
                                     usage='python update_genes_from_remote.py [-r remote_server_url]')
    parser.add_argument('-r', '--remote-server', default='', required=True,
                        help='base URL of the remote server')
    parser.add_argument('-k', '--ncbi-api-key', default=None, required=True,
                        help='NCBI Entrez API key')
    args = parser.parse_args()
    remote_addr = args.remote_server
    if not re.search(r'\w+', args.ncbi_api_key):
        sys.exit('ERROR: Invalid NCBI API key, please check')
    else:
        ncbi_api_key = args.ncbi_api_key
    # args = parser.parse_args(['-np'])
    print()
    log('INFO', 'Working with server {}'.format(remote_addr))

    # headers
    header = {
        'Accept': 'application/json',
        'User-Agent': 'python-requests Python/{}.{}.{}'.format(sys.version_info[0], sys.version_info[1], sys.version_info[2]),
    }
    # get db connector and cursor
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    # get local list of genes
    curs.execute(
        """
        SELECT DISTINCT(name[1]), second_name, chr, hgnc_id, ng, strand, hgnc_name AS gene_symbol
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
        log('DEBUG', 'req_url:{0}'.format(req_url))
        # log('DEBUG', 'api_response:{0}'.format(api_response))
        for key in api_response:
            # check general features such as hgnc id, second_name, strand and update all transcript of this gene at once for these features
            update_sql_gene = None
            if key == 'Chr':
                if api_response[key] != gene['chr']:
                    update_sql_gene += ' chr = {0} '.format(api_response[key])
            if key == 'HGNCID':
                if api_response[key] != gene['hgnc_id']:
                    update_sql_gene += ' hgnc_id = {0} '.format(api_response[key])
            if key == 'HGNCName':
                if api_response[key] != gene['hgnc_name']:
                    update_sql_gene += ' hgnc_name = {0} '.format(api_response[key])
            if key == 'HGNCSymbol':
                if api_response[key] != gene['name'][0]:
                    update_sql_gene += ' name[1] = {0} '.format(api_response[key])
            if key == 'RefGene':
                if api_response[key] != gene['ng'][0]:
                    update_sql_gene += ' ng = {0} '.format(api_response[key])
            if key == 'Strand':
                if api_response[key] != gene['strand'][0]:
                    update_sql_gene += ' strand = {0} '.format(api_response[key])
            # update all transcripts at once
            if update_sql_gene:
                update_sql_gene = "UPDATE gene SET {0} WHERE name[1] = '{1}'".format(update_sql_gene.join(','), gene['name'][0])
                log('DEBUG', update_sql_gene)
                # curs.execute(update_sql_gene)
                # db.commit()
            # transcripts
            ncbi_transcript_regexp = md_utilities.regexp.ncbi_transcript
            nm_obj = re.search(rf'^{ncbi_transcript_regexp}$', key)
            if nm_obj:
                # check if already exists
                # if yes => update if different => check each attribute? if one is differet, update all?
                # if no create transcript
                curs.execute(
                    "SELECT * FROM gene WHERE name[2] = %s",
                    (nm_obj.group(1))
                )
                prod_nm = curs.fetchone()
                if prod_nm:
                    # check if need to update
                    update_sql_transcript = None
                    if key == 'refProtein':
                        if api_response[key] != prod_nm['np']:
                            update_sql_transcript += ' np = {0} '.format(api_response[key])
                    if key == 'UNIPROT':
                        if api_response[key] != prod_nm['uniprot_id']:
                            update_sql_transcript += ' uniprot_id = {0} '.format(api_response[key])
                    if key == 'proteinSize':
                        if api_response[key] != prod_nm['prot_size']:
                            update_sql_transcript += ' prot_size = {0} '.format(api_response[key])
                    if key == 'canonical':
                        if api_response[key] != prod_nm['canonical']:
                            update_sql_transcript += ' canonical = {0} '.format(api_response[key])
                    if key == 'ensemblProtein':
                        if api_response[key] != prod_nm['ensp']:
                            update_sql_transcript += ' ensp = {0} '.format(api_response[key])
                    if key == 'ensemblTranscript':
                        if api_response[key] != prod_nm['enst']:
                            update_sql_transcript += ' enst = {0} '.format(api_response[key])
                    if key == 'numberOfExons':
                        if api_response[key] != prod_nm['number_of_exons']:
                            update_sql_transcript += ' number_of_exons = {0} '.format(api_response[key])
                    if key == 'variantCreationTag':
                        if api_response[key] != prod_nm['variantCreationTag']:
                            update_sql_transcript += ' variantCreationTag = {0} '.format(api_response[key])
                    # update ?
                    if update_sql_transcript:
                        update_sql_transcript = "UPDATE gene SET {0} WHERE name[2] = '{1}'".format(update_sql_transcript.join(','), prod_nm)
                        log('DEBUG', update_sql_transcript)
                        # curs.execute(update_sql_transcript)
                        # db.commit()
                else:
                    # insert
                    insert_gene = """
                    INSERT INTO gene (name, second_name, chr, strand, number_of_exons, hgnc_name, prot_size, uniprot_id, ng, np, enst, ensp, canonical, variant_creation, hgnc_id )
                    VALUES ('{{\"{0}\",\"{1}\"}}', '{2}', '{3}', '{4}', {5}, '{6}', {7}, '{8}', '{9}', '{10}', '{11}', '{12}', '{13}', '{14}', {15})
                    """.format(
                        api_response['HGNCSymbol'],
                        nm_obj,
                        gene['second_name'],
                        api_response['Chr'],
                        api_response['Strand'],
                        api_response['numberOfExons'],
                        api_response['HGNCName'],
                        api_response['proteinSize'],
                        api_response['UNIPROT'],
                        api_response['RefGene'],
                        api_response['RefProtein'],
                        api_response['ensemblTranscript'],
                        api_response['ensemblProtein'],
                        api_response['canonical'],
                        api_response['variantCreationTag'],
                        api_response['HGNCID']
                    ).replace("'NULL'", "NULL")
                    log('DEBUG', insert_gene)
                    # curs.execute(insert_gene)
                    # db.commit()
        print('.', end="", flush=True)
    db.close()


if __name__ == '__main__':
    main()
