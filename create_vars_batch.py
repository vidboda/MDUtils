import os
import sys
import re
import argparse
import urllib3
import certifi
import json
import time
import urllib.parse

# get a txt file list of variants in the format NM_XXXXX.Y:c.HGVScdna or NC_XXXXXX.Y:g.HGVChg38genomic;HGNCgeneName
# and calls MD API to create the variants


def log(level, text):
    localtime = time.asctime( time.localtime(time.time()) )
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


def main():
    parser = argparse.ArgumentParser(
        description='Creates variants in MobiDetails (batch mode)',
        usage='python create_vars_batch.py -l path/to/variant_list.txt -k api_key'
    )
    parser.add_argument('-l', '--varlist', default='', required=True,
                        help='Path to the text (.txt) file containing the list of variants\
                             in the format NM_XXXXX.Y:c.hgvscdna.')
    parser.add_argument('-k', '--api-key', default='', required=True,
                        help='Your API key visible on your profile page on the website.')
    parser.add_argument('-u', '--url', default='',
                        help='actions the provided server url instead of the base production server.')
    parser.add_argument('-f', '--format', default='transcript',
                        help="Whether the input format is NM_XXXXX.Y:c.HGVScdna or NC_XXXXXX.Y:g.HGVChg38genomic;HGNC. Values 'transcript' or 'genomic'.")
    args = parser.parse_args()
    # get file
    if os.path.isfile(args.varlist):
        batchFile = args.varlist
    else:
        log('ERROR', 'Invalid input path, please check your command')

    if len(args.api_key) != 43:
        log('ERROR', 'Invalid API key, please check it')
    else:
        api_key = args.api_key

    md_base_url = 'https://mobidetails.iurc.montp.inserm.fr/MD'
    if args.url:
        md_base_url = args.url
    input_format = 'transcript'
    if args.format and \
           args.format == 'genomic':
        input_format = 'genomic'
    #headers
    header = {
        'Accept': 'application/json',
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_3) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/35.0.1916.47 Safari/537.36'
    }
        
    print()
    log('INFO', 'Working on server {0} with format {1}'.format(md_base_url, input_format))

    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    for variant in open(batchFile).readlines():
        variant = variant.rstrip()
        variant = variant.replace(' ', '')
        if variant != '':
            match_obj = None
            if input_format == 'transcript':
                match_obj = re.search(r'^(NM_\d+)\.(\d{1,2}):(c\..+)$', variant)
            if input_format == 'genomic':
                match_obj = re.search(r'^([Nn][Cc]_\d+\.\d{1,2}:g\.[\dATGCagctdelinsup_>]+);(.+)$', variant)
            if match_obj is not None:
                # check if RefSeq accession number and version is suitable for variant creation:
                if input_format == 'transcript':
                    acc_number = match_obj.group(1)
                    acc_version = match_obj.group(2)
                    var = match_obj.group(3)
                    md_url_check = '{0}/api/gene/{1}'.format(md_base_url, match_obj.group(1))
                    try:
                        md_check_response = json.loads(http.request('GET', md_url_check, headers=header).data.decode('utf-8'))
                        for flag in md_check_response:
                            match_response = re.search(r'(NM_\d+)\.(\d)', flag)
                            if match_response and \
                                    match_response.group(2) != acc_version and \
                                    match_response.group(1) == acc_number:
                                # get correct accession version
                                acc_version = match_response.group(2)
                    except Exception:
                        log('WARNING', 'RefSeq accession number {} is not available in MobiDetails'.format(match_obj.group(1)))
                        next

                    # md_url = '{0}/api/variant/create/{1}/{2}'.format(md_base_url, variant, api_key)
                    md_url = '{0}/api/variant/create/{1}.{2}:{3}/{4}'.format(md_base_url, urllib.parse.quote(acc_number), urllib.parse.quote(acc_version), urllib.parse.quote(var), api_key)
                    log('INFO', 'Submitting variant {0}.{1}:{2} to MobiDetails: {3}'.format(acc_number, acc_version, var, md_url))
    
                    try:
                        md_response = json.loads(http.request('GET', md_url, headers=header).data.decode('utf-8'))
                    except Exception:
                        log('WARNING', 'VariantValidator call failed for variant {}'.format(variant))
                        continue
                elif input_format == 'genomic':
                    g_var = match_obj.group(1)
                    gene = match_obj.group(2)
                    md_url_check = '{0}/api/gene/{1}'.format(md_base_url, gene)
                    semaph = 'no'
                    try:
                        md_check_response = json.loads(http.request('GET', md_url_check, headers=header).data.decode('utf-8'))
                        for flag in md_check_response:
                            if flag == 'HGNC':
                                semaph = 'yes'
                                continue
                    except Exception:
                        log('WARNING', 'The gene {} could not be checked for some reason in MobiDetails'.format(gene))
                    if semaph == 'yes':
                        md_url = '{0}/api/variant/create_g/{1}/{2}/cli/{3}'.format(md_base_url, g_var, gene, api_key)
                        log('INFO', 'Submitting variant {0} for gene in {1} to MobiDetails: {2}'.format(g_var, gene, md_url))
                        try:
                            md_response = json.loads(http.request('GET', md_url, headers=header).data.decode('utf-8'))
                        except Exception:
                            log('WARNING', 'VariantValidator call failed for genomic variant {}'.format(variant))
                            continue
                    else:
                        log('WARNING', 'the gene {} is unknown in MobiDetails'.format(gene))
                        continue
                log('INFO', md_response)
                print()
            else:
                 log('ERROR', "the format {} does not seem to fit with your file. If you choose the transcript format, you file whould be of type: NM_XXXXX.Y:c.HGVScdna, one variant per line, and if the format is 'genomic', your file should look like NC_XXXXXX.Y:g.HGVChg38genomic;HGNCgeneName".format(input_format))


if __name__ == '__main__':
    main()
