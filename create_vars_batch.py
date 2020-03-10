import os
import sys
import re
import argparse
import urllib3
import certifi
import json

# get a txt file list of variants in the format NM_XXXXX.Y:c.hgvscdna
# and calls MD API to create the variants


def log(level, text):
    if level == 'ERROR':
        sys.exit('[{0}]: {1}'.format(level, text))
    print('[{0}]: {1}'.format(level, text))


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
    print()
    log('INFO', 'Working on server {}'.format(md_base_url))

    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    for variant in open(batchFile).readlines():
        variant = variant.rstrip()
        variant = variant.replace(' ', '')
        if variant != '':
            match_obj = re.search(r'(NM_\d+)\.(\d):(c\..+)', variant)
            if match_obj:
                # check if RefSeq accession number and version is suitable for variant creation:
                acc_number = match_obj.group(1)
                acc_version = match_obj.group(2)
                var = match_obj.group(3)
                md_url_check = '{0}/api/gene/{1}'.format(md_base_url, match_obj.group(1))
                try:
                    md_check_response = json.loads(http.request('GET', md_url_check).data.decode('utf-8'))
                    for flag in md_check_response:
                        match_response = re.search(r'(NM_\d+)\.(\d)', flag)
                        if match_response and \
                                match_response.group(2) != acc_version and \
                                match_response.group(1) == acc_number:
                            # get good accession version
                            acc_version = match_response.group(2)
                except Exception:
                    log('WARNING', 'RefSeq accession number {} is not available in MobiDetails'.format(match_obj.group(1)))
                    next

                # md_url = '{0}/api/variant/create/{1}/{2}'.format(md_base_url, variant, api_key)
                md_url = '{0}/api/variant/create/{1}.{2}:{3}/{4}'.format(md_base_url, acc_number, acc_version, var, api_key)
                log('INFO', 'Submitting variant {0}.{1}:{2} to MobiDetails: {3}'.format(acc_number, acc_version, var, md_url))

                try:
                    md_response = json.loads(http.request('GET', md_url).data.decode('utf-8'))
                except Exception:
                    log('WARNING', 'VariantValidator call failed for variant {}'.format(variant))
                    continue

                log('INFO', md_response)
                print()


if __name__ == '__main__':
    main()
