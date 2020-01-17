import os
import sys
import re
import argparse
import urllib3
import certifi
import json

#get a txt file list of variants in the format NM_XXXXX.Y:c.hgvscdna
#and calls MD API to create the variants

def main():
	parser = argparse.ArgumentParser(description='Creates variants in MobiDetails (batch mode)', usage='python create_vars_batch.py -l path/to/variant_list.txt')
	parser.add_argument('-l', '--varlist', default='', required=True, help='Path to the text (.txt) file containing the list of variants in the format NM_XXXXX.Y:c.hgvscdna.')
	parser.add_argument('-k', '--api-key', default='', required=True, help='Your API key visible on your profile page on the website.')
	parser.add_argument('-u', '--url', default='', help='actions the provided server url instead of the base production server.')
	args = parser.parse_args()
	#get file 
	if os.path.isfile(args.varlist):
		batchFile = args.varlist
	else:
		sys.exit('ERROR: Invalid input path, please check your command')
	
	if len(args.api_key) != 43:
		sys.exit('ERROR: Invalid API key, please check it')
	else:
		api_key = args.api_key
	
	md_base_url = 'https://mobidetails.iurc.montp.inserm.fr/MD'
	if args.url:
		md_base_url = args.url
	print()
	print('INFO: Working on server {}'.format(md_base_url))
	
	http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
	for variant in open(batchFile).readlines():
		variant = variant.rstrip()
		variant = variant.replace(' ', '')
		if variant != '':
			if re.search('NM_\d+\.\d:c\..+', variant):
				md_url = '{0}/api/variant/create/{1}/{2}'.format(md_base_url, variant, api_key)
				print('INFO: Submitting variant {0} to MobiDetails: {1}'.format(variant, md_url))
				
				try:
					md_response = json.loads(http.request('GET', md_url).data.decode('utf-8'))
				except:
					print('ERROR: Call failed for variant {}'.format(variant))
					continue
					
				print(md_response)
				print()

if __name__ == '__main__':
	main()