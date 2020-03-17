import os
import sys
import re
import argparse
import urllib3
import certifi
import requests
import json
# import hashlib

# connect to a distant resource and check whether or
# not we should update - and updates


def log(level, text):
    if level == 'ERROR':
        sys.exit('[{0}]: {1}'.format(level, text))
    print('[{0}]: {1}'.format(level, text))


def get_last_clinvar_md5_file(resources_dir):
    files = os.listdir(resources_dir)
    dates = []
    for current_file in files:
        # print(current_file)
        match_obj = re.search(r'clinvar_(\d+).vcf.gz.md5$', current_file)
        if match_obj:
            dates.append(match_obj.group(1))
    current_clinvar = '{0}clinvar_{1}.vcf.gz.md5'.format(resources_dir, max(dates))
    with open(current_clinvar, 'r') as clinvar_file:
        # print(clinvar_file.read())
        match_obj = re.search(r'^(\w+)\s', clinvar_file.read())
        if match_obj:
            return match_obj.group(1)
    return 'f'

# from https://www.techcoil.com/blog/how-to-download-a-file-via-http-post-and-http-get-with-python-3-requests-library/


def download_file_from_server_endpoint(server_endpoint, local_file_path): 
    # Send HTTP GET request to server and attempt to receive a response
    response = requests.get(server_endpoint)     
    # If the HTTP GET request can be served
    if response.status_code == 200:         
        # Write the file contents in the response to a file specified by local_file_path
        log('INFO', 'Downloading Clinvar file as {}'.format(local_file_path))
        with open(local_file_path, 'wb') as local_file:
            for chunk in response.iter_content(chunk_size=128):
                local_file.write(chunk)


def main():
    parser = argparse.ArgumentParser(
        description='Checks for MD resources distant updates',
        usage='python update_resources.py -c path/to/variant_list.txt -k api_key'
    )
    parser.add_argument('-c', '--clinvar', default='', required=False,
                        help='Optionally updates clinvar vcf', action='store_true')
    parser.add_argument('-d', '--dbsnp', default='', required=False,
                        help='Optionally updates dbsnp vcf')
    parser.add_argument('-ip', '--prod-ip', default='194.167.35.207', required=False,
                        help='Production server IP for rsync copy')
    args = parser.parse_args()
    clinvar_url = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/' 
    
    dbsnp_url = ''
    resources_path = '/home/adminbioinfo/Devs/MobiDetails/MobiDetailsApp/static/resources/'
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())

    if args.clinvar:
        distant_md5 = None
        download_semaph = None
        clinvar_dir_content = None
        try:
            # Get all file names from clinvar website (html)
            clinvar_dir_html = http.request('GET', clinvar_url).data.decode('utf-8')
            clinvar_dir_content = re.split('\n', clinvar_dir_html)
            for html in clinvar_dir_content:
                match_obj = re.search(r'\"clinvar_(\d+).vcf.gz\"', html)
                if match_obj:
                    break
                # log('DEBUG', html)
            # log('DEBUG', match_obj)
        except Exception:
            log('WARNING', 'Unable to contact clinvar {}'.format(clinvar_url))
        if match_obj:
            # Get current clinvar date            
            # log('DEBUG', match_obj)
            if match_obj:
                clinvar_date = match_obj.group(1)
                # Read current clinvar md5
                try:
                    clinvar_md5 = http.request('GET', '{0}clinvar_{1}.vcf.gz.md5'.format(clinvar_url, clinvar_date)).data.decode('utf-8')
                    match_obj = re.search(r'^(\w+)\s', clinvar_md5)
                    if match_obj:
                        distant_md5 = match_obj.group(1)
                        log('INFO', 'Distant md5: {}'.format(distant_md5))
                except Exception:
                    log('WARNING', 'Unable to contact clinvar {}'.format(clinvar_url))
                if distant_md5:
                    # Get md5 from local file
                    current_md5_value = get_last_clinvar_md5_file('{}clinvar/hg38/'.format(resources_path))
                    log('INFO', 'Local md5: {}'.format(current_md5_value))
                    if current_md5_value != distant_md5:
                        # Download remote file
                        try:
                            download_file_from_server_endpoint('{0}clinvar_{1}.vcf.gz'.format(clinvar_url, clinvar_date), '{0}clinvar/hg38/clinvar_{1}.vcf.gz'.format(resources_path, clinvar_date))
                            download_file_from_server_endpoint('{0}clinvar_{1}.vcf.gz.md5'.format(clinvar_url, clinvar_date), '{0}clinvar/hg38/clinvar_{1}.vcf.gz.md5'.format(resources_path, clinvar_date))
                            download_file_from_server_endpoint('{0}clinvar_{1}.vcf.gz.tbi'.format(clinvar_url, clinvar_date), '{0}clinvar/hg38/clinvar_{1}.vcf.gz.tbi'.format(resources_path, clinvar_date))
                            download_semaph = 1
                        except Exception:
                            log('WARNING', 'Unable to download new clinvar file {0}clinvar_{1}.vcf.gz'.format(clinvar_url, clinvar_date))
                        # Then check clinvar md5 and test w/ a variant


if __name__ == '__main__':
    main()


# if we need to check a md5
# with open('MobiDetailsApp/static/resources/clinvar/hg38/clinvar_20200310.vcf.gz', 'rb') as clinvar_file:
#      buf = clinvar_file.read(BLOCKSIZE)
#      while len(buf) > 0:
#          hasher.update(buf)
#          buf = clinvar_file.read(BLOCKSIZE)
# 
# print(hasher.hexdigest())
