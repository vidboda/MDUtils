import os
import argparse
import urllib3
import certifi
import json

#get a txt file list of variants in the format NM_XXXXX.Y:c.hgvscdna
#and calls MD API to create the variants

def main():
	parser = argparse.ArgumentParser(description='Creates variants in MobiDetails (batch mode)', usage='python create_vars_batch.py -l path/to/variant_list.txt')
	parser.add_argument('-l', '--varlist', default='', required=True, help='Path to the text (.txt) file containing the list of variants in the format NM_XXXXX.Y:c.hgvscdna')
	parser.add_argument('-u', '--url', default='', help='actions the provided server url instead of the base production server.')
	args = parser.parse_args()
	#get file 
	if os.path.isfile(args.varlist):
		batchFile = args.varlist
	else:
		sys.exit('Invalid input path, please check your command')
	
	md_base_url = 'https://mobidetails.iurc.montp.inserm.fr/MD'
	if args.url:
		md_base_url = args.url
	print()
	print('Working on server {}'.format(md_base_url))
	
	http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
	for variant in open(batchFile).readlines():
		variant = variant.rstrip()
		md_url = '{0}/api/variant/create/{1}'.format(md_base_url, variant)
		print('Submitting variant {0} to MobiDetails: {1}'.format(variant, md_url))
		
		#try:
		md_response = json.loads(http.request('GET', md_url).data.decode('utf-8'))
		#except:
		#	print('Call failed')
		#	continue
			
		print(md_response)
		print()

if __name__ == '__main__':
	main()