import os
import sys
import re
import argparse
import urllib3
import certifi
import requests
import hashlib
import time
import datetime
import json
import shutil

# connect to a distant resource and check whether or
# not we should update - and updates


def log(level, text):
    localtime = time.asctime(time.localtime(time.time()))
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


def get_last_md5_file(resource_dir, resource_type, resource_regexp, target_suffix, suffix):
    files = os.listdir(resource_dir)
    dates = []
    for current_file in files:
        # print(current_file)
        match_obj = re.search(rf'{resource_regexp}{target_suffix}{suffix}.md5$', current_file)
        if match_obj:
            dates.append(match_obj.group(1))
    current_resource = '{0}{1}_{2}{3}{4}.md5'.format(resource_dir, resource_type, max(dates), target_suffix, suffix)
    # log('DEBUG', current_resource)
    with open(current_resource, 'r') as current_file:
        # print(clinvar_file.read())
        match_obj = re.search(r'^(\w+)\s', current_file.read())
        if match_obj:
            return match_obj.group(1)
    return 'f'

# from https://www.techcoil.com/blog/how-to-download-a-file-via-http-post-and-http-get-with-python-3-requests-library/


def download_file_from_server_endpoint(server_endpoint, local_file_path):
    # Send HTTP GET request to server and attempt to receive a response
    response = requests.get(server_endpoint)
    # If the HTTP GET request can be served
    # log('DEBUG', response.status_code)
    if response.status_code == 200:
        # Write the file contents in the response to a file specified by local_file_path
        try:
            log('INFO', 'Downloading file as {}'.format(local_file_path))
            with open(local_file_path, 'wb') as local_file:
                for chunk in response.iter_content(chunk_size=128):
                    local_file.write(chunk)
            # log('DEBUG', 'Downloaded file as {}'.format(local_file_path))
        except Exception:
            log('WARNING', 'Unable to download {}'.format(server_endpoint))
    else:
        log('WARNING', 'Unable to contact {}'.format(server_endpoint))


def get_new_ncbi_resource_file(http, resource_type, resource_dir, regexp, label, url, target_suffix):
    distant_md5 = None
    download_semaph = None
    resource_dir_content = None
    # log('DEBUG', '{0}{1}.gz'.format(regexp, target_suffix))
    try:
        # Get all file names from clinvar website (html)
        resource_dir_html = http.request('GET', url).data.decode('utf-8')
        resource_dir_content = re.split('\n', resource_dir_html)
        for html in resource_dir_content:
            match_obj = re.search(rf'\"{regexp}{target_suffix}.gz\"', html)
            if match_obj:
                break
            # log('DEBUG', html)
        # log('DEBUG', match_obj)
    except Exception:
        log('WARNING', 'Unable to contact {0} {1}'.format(label, url))
    if match_obj:
        # Get current resource date
        # log('DEBUG', match_obj)
        if match_obj:
            resource_date = match_obj.group(1)
            # Read current clinvar md5
            try:
                resource_md5 = http.request('GET', '{0}{1}_{2}{3}.gz.md5'.format(url, resource_type, resource_date, target_suffix)).data.decode('utf-8')
                match_obj = re.search(r'^(\w+)\s', resource_md5)
                if match_obj:
                    distant_md5 = match_obj.group(1)
                    log('INFO', '{0} distant md5: {1}'.format(label, distant_md5))
            except Exception:
                log('WARNING', 'Unable to contact {0} md5 {1}'.format(label, url))
            if distant_md5:
                # Get md5 from local file
                # current_md5_value = get_last_clinvar_md5_file('{}clinvar/hg38/'.format(resources_path))
                current_md5_value = get_last_md5_file(resource_dir, resource_type, regexp, target_suffix, '.gz')
                log('INFO', '{0} local md5: {1}'.format(label, current_md5_value))
                exit
                if current_md5_value != distant_md5:
                    # Download remote file
                    # log('DEBUG', '{0}{1}_{2}{3}.gz'.format(url, resource_type, resource_date, target_suffix))
                    # log('DEBUG', resource_dir)
                    # log('DEBUG', '{0}{1}_{2}{3}.gz'.format(resource_dir, resource_type, resource_date, target_suffix))
                    # try:
                    download_file_from_server_endpoint(
                        '{0}{1}_{2}{3}.gz'.format(url, resource_type, resource_date, target_suffix),
                        '{0}{1}_{2}{3}.gz'.format(resource_dir, resource_type, resource_date, target_suffix)
                    )
                    download_file_from_server_endpoint(
                        '{0}{1}_{2}{3}.gz.md5'.format(url, resource_type, resource_date, target_suffix),
                        '{0}{1}_{2}{3}.gz.md5'.format(resource_dir, resource_type, resource_date, target_suffix)
                    )
                    download_file_from_server_endpoint(
                        '{0}{1}_{2}{3}.gz.tbi'.format(url, resource_type, resource_date, target_suffix),
                        '{0}{1}_{2}{3}.gz.tbi'.format(resource_dir, resource_type, resource_date, target_suffix)
                    )
                    download_semaph = 1
                    # except Exception:
                    #     log('WARNING', 'Unable to download new {0} file {1}{2}_{3}{4}.gz'.format(label, url, resource_type, resource_date, target_suffix))
                    # Then check new md5 and test w/ a variant
                    if download_semaph == 1:
                        with open('{0}{1}_{2}{3}.gz'.format(resource_dir, resource_type, resource_date, target_suffix), 'rb') as new_resource_file:
                            BLOCKSIZE = 65536
                            buf = new_resource_file.read(BLOCKSIZE)
                            hasher = hashlib.md5()
                            while len(buf) > 0:
                                hasher.update(buf)
                                buf = new_resource_file.read(BLOCKSIZE)
                            # log('DEBUG', hasher.hexdigest() )
                            if hasher.hexdigest() == distant_md5:
                                # Download successful
                                log('INFO', 'Successfully downloaded and checked {0} file {1}_{2}{3}.gz'.format(label, resource_type, resource_date, target_suffix))
                            else:
                                # Remove file
                                os.remove('{0}{1}_{2}{3}.gz'.format(resource_dir, resource_type, resource_date, target_suffix))
                                os.remove('{0}{1}_{2}{3}.gz.md5'.format(resource_dir, resource_type, resource_date, target_suffix))
                                os.remove('{0}{1}_{2}{3}.gz.tbi'.format(resource_dir, resource_type, resource_date, target_suffix))
                                log('WARNING', 'Error in md5 sum for {0} file {0}_{1}{2}.gz'.format(resource_type, resource_date, target_suffix))


def main():
    parser = argparse.ArgumentParser(
        description='Checks for MD resources distant updates',
        usage='python update_resources.py <-c>'
    )
    parser.add_argument('-c', '--clinvar', default='', required=False,
                        help='Optionally updates clinvar vcf', action='store_true')
    parser.add_argument('-cg', '--clingen-criteria', default='', required=False,
                        help='Optionally updates clingen criteria summary json and generates a txt file containing the list of concerned genes', action='store_true')
    parser.add_argument('-d', '--dbsnp', default='', required=False,
                        help='Optionally updates dbsnp vcf', action='store_true')
    args = parser.parse_args()

    clinvar_url = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/'
    dbsnp_url = 'https://ftp.ncbi.nih.gov/snp/latest_release/'
    clingen_criteria_url = 'https://cspec.genome.network/cspec/CriteriaCode/summary'
    resources_path = '/home/adminbioinfo/Devs/MobiDetails/MobiDetailsApp/static/resources/'

    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    match_obj = None
    if args.clinvar:
        # get_new_resource_file(resource_type, resources_dir, regexp, label, url, target_path, suffix, http object)
        get_new_ncbi_resource_file(http, 'clinvar', '{}clinvar/hg38/'.format(resources_path), r'clinvar_(\d+)', 'ClinVar', clinvar_url, '.vcf')
    if args.clingen_criteria:
        # current date
        current_date = datetime.datetime.today()
        date_string = '{0}{1}{2}'.format(current_date.year, current_date.strftime('%m'), current_date.strftime('%d'))
        clingen_new_file = 'clingenCriteriaSpec_{}'.format(date_string)
        download_file_from_server_endpoint(
            clingen_criteria_url,
            '{0}clingen/{1}.json'.format(resources_path, clingen_new_file)
        )
        # we need to md5 the new file and the current one
        current_md5_value = get_last_md5_file('{}clingen/'.format(resources_path), 'clingenCriteriaSpec', r'clingenCriteriaSpec_(\d+)', '', '')
        BLOCKSIZE = 65536
        hasher = hashlib.md5()
        with open('{0}clingen/{1}.json'.format(resources_path, clingen_new_file), 'rb') as clingen_file:
            buf = clingen_file.read(BLOCKSIZE)
            while len(buf) > 0:
                hasher.update(buf)
                buf = clingen_file.read(BLOCKSIZE)
        clingen_file.close()
        if hasher.hexdigest() != current_md5_value:
            log('INFO', 'New Clingen file found')
            # create the md5 file
            new_md5 = open('{0}clingen/{1}.md5'.format(resources_path, clingen_new_file), "w")
            new_md5.write('{0} {1}'.format(hasher.hexdigest(), clingen_new_file))
            new_md5.close()
            log('INFO', 'New Clingen file ready and hashed')
            # generate txt file containing the gene symbols list ("index")
            with open('{0}clingen/{1}.json'.format(resources_path, clingen_new_file), 'r') as clingen_file:
                clingen_json = json.load(clingen_file)
            genes_symbols_file = open('{0}clingen/{1}.txt'.format(resources_path, clingen_new_file), "w")
            if 'data' in clingen_json:
                for rule in clingen_json['data']:
                    if 'genes' in rule:
                        for gene in rule['genes']:
                            if 'label' in gene:
                                genes_symbols_file.write('{0}\n'.format(gene['label']))
                            else:
                                log('WARNING', 'No gene symbol in gene key {0} from rule {1}'.format(gene, rule))
                    else:
                        log('WARNING', 'No gene in key {0}'.format(rule))
            else:
                log('ERROR', 'No data in the new Cligen File - erasing')
                os.remove('{0}clingen/{1}.md5'.format(resources_path, clingen_new_file))
                os.remove('{0}clingen/{1}'.format(resources_path, clingen_new_file))
            genes_symbols_file.close()
            shutil.copyfile('{0}clingen/{1}.json'.format(resources_path, clingen_new_file), '{0}clingen/clingenCriteriaSpec_last.json'.format(resources_path))
            shutil.copyfile('{0}clingen/{1}.md5'.format(resources_path, clingen_new_file), '{0}clingen/clingenCriteriaSpec_last.md5'.format(resources_path))
            shutil.copyfile('{0}clingen/{1}.txt'.format(resources_path, clingen_new_file), '{0}clingen/clingenCriteriaSpec_last.txt'.format(resources_path))
            log('INFO', 'New files processed and ready for use.')

    if args.dbsnp:
        # get dbsnp version from https://ftp.ncbi.nih.gov/snp/latest_release/release_notes.txt
        download_file_from_server_endpoint(
            '{}release_notes.txt'.format(dbsnp_url),
            '{}dbsnp/release_notes.txt'.format(resources_path)
        )
        with open('{}dbsnp/release_notes.txt'.format(resources_path), 'r') as f:
            match_obj = re.search(r'dbSNP build (\d+) release notes', f.readline())
            semaph = 0
            if match_obj:
                dbsnp_version = match_obj.group(1)
                log('INFO', 'dbSNP version file found: v{0}'.format(dbsnp_version))
                if not os.path.exists('{0}/dbsnp/hg38/v{1}'.format(resources_path, dbsnp_version)):
                    os.makedirs('{0}/dbsnp/hg38/v{1}'.format(resources_path, dbsnp_version))
                else:
                    log('INFO', 'dbSNP version file found: v{0} same as current'.format(dbsnp_version))
                    semaph = 1
                # os.mkdir('{0}/dbsnp/v{1}'.format(resources_path, dbsnp_version))
            else:
                log('ERROR', 'Unable to donwload/read dbSNP release file from {}release_notes.txt'.format(dbsnp_url))
            if semaph == 0:
                # we can proceed
                get_new_ncbi_resource_file(http, 'GCF', '{0}dbsnp/hg38/v{1}/'.format(resources_path, dbsnp_version), r'GCF_(\d+)', 'dbSNP', '{0}VCF/'.format(dbsnp_url), '.38')


if __name__ == '__main__':
    main()

# From https://www.pythoncentral.io/hashing-files-with-python/
# if we need to check a md5
# BLOCKSIZE = 65536
# hasher = hashlib.md5()
# with open('MobiDetailsApp/static/resources/clinvar/hg38/clinvar_20200310.vcf.gz', 'rb') as clinvar_file:
#      buf = clinvar_file.read(BLOCKSIZE)
#      while len(buf) > 0:
#          hasher.update(buf)
#          buf = clinvar_file.read(BLOCKSIZE)

# print(hasher.hexdigest())
