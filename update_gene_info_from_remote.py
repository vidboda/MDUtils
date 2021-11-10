import re
import sys
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import json
import time
from insert_genes import get_db
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config

# fix genes in database did not have canonical transcripts.
# fix from remote MD server which has the information (typically dev server), using the API
# fix UNIPORT IDs, using the MDAPI and Uniprot API


def log(level, text):
    localtime = time.asctime( time.localtime(time.time()) )
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


def main():
    parser = argparse.ArgumentParser(description='Define a canonical transcript per gene and optionally updates various fields',
                                     usage='python update_canonical_from_remote.py [-r remote_server_url]')
    parser.add_argument('-r', '--remote-server', default='', required=True,
                        help='base URL of the remote server')
    parser.add_argument('-uca', '--update-can-all', default='', required=False,
                        help='Optionally update canonical for all genes w/ no variants', action='store_true')
    parser.add_argument('-np', '--update-np', default='', required=False,
                        help='Optionally update NP for genes lacking NP', action='store_true')
    parser.add_argument('-uu', '--update-uniprot', default='', required=False,
                        help='Optionally update UNIPROT IDs', action='store_true')
    parser.add_argument('-uc', '--update-creation', default='', required=False,
                        help='Optionally update variant_creation tag', action='store_true')
    parser.add_argument('-un', '--update-nm', default='', required=False,
                        help='Optionally update RefSeq NM accession number tag', action='store_true')
    parser.add_argument('-npf', '--update-np-full', default='', required=False,
                        help='Optionally update NP for all genes', action='store_true')
    parser.add_argument('-ue', '--update-exons', default='', required=False,
                        help='Optionally update number of exons for all genes', action='store_true')
    parser.add_argument('-ugn', '--update-gene-names', default='', required=False,
                        help='Optionally update HGNC gene names for all genes', action='store_true')

    args = parser.parse_args()
    remote_addr = args.remote_server
    # args = parser.parse_args(['-np'])
    print()
    log('INFO', 'Working with server {}'.format(remote_addr))

    #headers
    header = {
        'Accept': 'application/json',
        'User-Agent': 'python-requests Python/{}.{}.{}'.format(sys.version_info[0], sys.version_info[1], sys.version_info[2]),
    }
    # get db connector and cursor
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    # get local list of genes with no canonical isoform defined
    curs.execute(
        "SELECT DISTINCT(name[1]) as hgnc FROM gene WHERE name[1] NOT IN (SELECT name[1] FROM gene WHERE canonical='t') ORDER BY name[1]"
    )
    no_can = curs.fetchall()

    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())

    i = 0
    # lacking_nm = []

    for gene in no_can:
        req_url = '{0}/api/gene/{1}'.format(remote_addr, gene['hgnc'])
        api_response = json.loads(http.request('GET', req_url, headers=header).data.decode('utf-8'))
        log('DEBUG', 'req_url:{0}'.format(req_url))
        # log('DEBUG', 'api_response:{0}'.format(api_response))
        for keys in api_response:
            # log('DEBUG', 'key:{0} - value:{1}'.format(keys, api_response[keys]))
            if isinstance(keys, dict) and \
                    'canonical' in api_response[keys]:
                if api_response[keys]['canonical'] is True:
                    if re.search(r'NM_\d+\.\d+', keys):
                        match_obj = re.search(r'(NM_\d+)\.\d+', keys)
                        nm_acc = match_obj.group(1)
                        curs.execute(
                            "UPDATE gene set canonical = 't' WHERE name[2] = %s",
                            (nm_acc,)
                        )
                        log('INFO', 'Updating {}'.format(nm_acc))
                        i += 1
    db.commit()
    log('INFO', '{} genes modified (canonical)'.format(i))
    i = 0
    if args.update_can_all:
        # get genes with no variants and at least 2 isoforms to see if we need to update canonical
        curs.execute(
            "SELECT name, canonical FROM gene WHERE (name[1] NOT IN (SELECT gene_name[1] FROM variant_feature)) \
            AND (name[1] IN (SELECT name[1] FROM gene GROUP BY name[1] HAVING COUNT (name[1]) > 1)) ORDER BY name"
        )
        res = curs.fetchall()
        for acc in res:
            req_url = '{0}/api/gene/{1}'.format(remote_addr, acc['name'][0])
            api_response = json.loads(http.request('GET', req_url, headers=header).data.decode('utf-8'))
            for keys in api_response:
                 if isinstance(keys, dict) and \
                        'canonical' in api_response[keys]:
                    if api_response[keys]['canonical'] is True and acc['canonical'] == 0:
                        if re.search(r'NM_\d+\.\d+', keys):
                            match_obj = re.search(r'(NM_\d+)\.\d+', keys)
                            nm_acc = match_obj.group(1)
                            # double check
                            if nm_acc == acc['name'][1]:
                                curs.execute(
                                    "UPDATE gene SET canonical = 'f' WHERE name[1] = %s",
                                    (acc['name'][0],)
                                )
                                # log('INFO', "UPDATE gene SET canonical = 'f' WHERE name[1] = '{}'".format(acc['name'][0]))
                                curs.execute(
                                    "UPDATE gene SET canonical = 't' WHERE name[2] = %s",
                                    (acc['name'][1],)
                                )
                                # log('INFO', "UPDATE gene SET canonical = 't' WHERE name[2] = '{}'".format(acc['name'][1]))
                                i += 1
                                log('INFO', 'Updated gene {}'.format(acc['name'][0]))
        db.commit()

        log('INFO', '{} genes modified (canonical all)'.format(i))

    if args.update_np:
        curs.execute(
            "SELECT DISTINCT(name[1]) as hgnc FROM gene WHERE np = 'NP_000000.0'"
        )
        no_np = curs.fetchall()
        j = 0
        for gene in no_np:
            req_url = '{0}/api/gene/{1}'.format(remote_addr, gene['hgnc'])
            api_response = json.loads(http.request('GET', req_url, headers=header).data.decode('utf-8'))
            for keys in api_response:
                if isinstance(keys, dict):
                    if 'RefProtein' in api_response[keys] and \
                            api_response[keys]['RefProtein'] != 'NP_000000.0':
                        if re.search(r'NP_\d+\.\d+', api_response[keys]['RefProtein']):
                            match_obj = re.search(r'(NM_\d+)\.\d+', keys)
                            nm_acc = match_obj.group(1)
                            np_acc = api_response[keys]['RefProtein']
                            curs.execute(
                                "UPDATE gene set np = %s WHERE name[2] = %s",
                                (np_acc, nm_acc)
                            )
                            log('INFO', 'Updating gene NP acc no of {0} to {1}'.format(nm_acc, np_acc))
                            j += 1
        db.commit()
        log('INFO', '{} NP acc no modified'.format(j))
    if args.update_uniprot or args.update_creation or args.update_nm or args.update_np_full or args.update_exons or args.update_gene_names:
        curs.execute(
            "SELECT  name[1] as HGNC, name[2] as nm, nm_version, np, uniprot_id, variant_creation FROM gene ORDER BY name"
        )
        res = curs.fetchall()
        k = l = m = n = p = q = r = 0
        o = curs.rowcount
        for gene in res:
            req_url = '{0}/api/gene/{1}'.format(remote_addr, gene['hgnc'])
            api_response = json.loads(http.request('GET', req_url, headers=header).data.decode('utf-8'))
            l += 1
            if l % 1000 == 0:
                log('INFO', '{0}/{1} isoforms checked'.format(l, o))
            for keys in api_response:
                if 'HGNC Name' in api_response[keys] and args.update_gene_names:
                    gene_name = api_response[keys]['HGNC Name']
                    if gene_name != gene['prot_name']:
                        curs.execute(
                            "UPDATE gene set prot_name = %s WHERE name[2] = %s",
                            (gene_name, nm_acc)
                        )
                        log('INFO', 'Updating gene name of {0} to {1}'.format(nm_acc, ungene_nameiprot))
                        r += 1
                match_obj = re.search(r'^(NM_\d+)\.(\d+)$', keys)
                if match_obj:
                    nm_acc = match_obj.group(1)
                    # check again
                    if nm_acc == gene['nm']:
                        if args.update_nm:
                            nm_version = match_obj.group(2)
                            # log('DEBUG', '{0}dev:{1}-prod:{2}'.format(gene['hgnc'], int(nm_version), int(gene['nm_version'])))
                            if int(nm_version) != int(gene['nm_version']):
                                # no downgrade? y => downgrade
                                curs.execute(
                                    "UPDATE gene set nm_version = %s WHERE name[2] = %s",
                                    (nm_version, nm_acc)
                                )
                                log('INFO', 'Updating gene RefSeq NM accession version of {0} from {1} to {2}'.format(nm_acc, gene['nm_version'], nm_version))
                                n += 1
                        if 'UNIPROT' in api_response[keys] and args.update_uniprot:
                            uniprot = api_response[keys]['UNIPROT']
                            if uniprot != gene['uniprot_id']:
                                curs.execute(
                                    "UPDATE gene set uniprot_id = %s WHERE name[2] = %s",
                                    (uniprot, nm_acc)
                                )
                                log('INFO', 'Updating gene UNIPROT id of {0} to {1}'.format(nm_acc, uniprot))
                                k += 1
                        if 'variantCreationTag' in api_response[keys] and args.update_creation:
                            tag = api_response[keys]['variantCreationTag']
                            if tag != gene['variant_creation']:
                                curs.execute(
                                    "UPDATE gene set variant_creation = %s WHERE name[2] = %s",
                                    (tag, nm_acc)
                                )
                                log('INFO', 'Updating gene variantCreationTag of {0} to {1}'.format(nm_acc, tag))
                                m += 1
                        if 'RefProtein' in api_response[keys] and args.update_np_full:
                            np = api_response[keys]['RefProtein']
                            if np != gene['np']:
                                curs.execute(
                                    "UPDATE gene set np = %s WHERE name[2] = %s",
                                    (np, nm_acc)
                                )
                                log('INFO', 'Updating gene np of {0} to {1}'.format(nm_acc, np))
                                p += 1
                        if 'total exons' in api_response[keys] and args.update_np_full:
                            exons = api_response[keys]['total exons']
                            if exons != gene['number_of_exons']:
                                curs.execute(
                                    "UPDATE gene set number_of_exons = %s WHERE name[2] = %s",
                                    (exons, nm_acc)
                                )
                                log('INFO', 'Updating gene total exons of {0} to {1}'.format(nm_acc, exons))
                                q += 1
        db.commit()
        
        log('INFO', '{} UNIPROT IDs modified'.format(k))
        log('INFO', '{} variantCreationTag modified'.format(m))
        log('INFO', '{} RefSeq NM accession version modified'.format(n))
        log('INFO', '{} NP version modified'.format(p))
        log('INFO', '{} Exons number modified'.format(q))
        log('INFO', '{} Gene names modified'.format(r))


if __name__ == '__main__':
    main()
