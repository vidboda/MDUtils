import re
import sys
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import json
from insert_genes import get_db
#requires MobiDetails config module + database.ini file
from MobiDetailsApp import config

#removes all c.1A>T then checks that these variants can be created in genes using API
#to be used on private dev server

def main():
	parser = argparse.ArgumentParser(description='Checks that genes accept variant creation', usage='python check_variant_creation.py [-r remote_server_url]')
	parser.add_argument('-r', '--remote-server', default='', required=True, help='base URL of the remote server')
	parser.add_argument('-k', '--api-key', default='', required=True, help='Your API key visible on your profile page on the website.')
	
	args = parser.parse_args()
	remote_addr = args.remote_server
	if re.search('mobidetails\.iurc', remote_addr):
		sys.exit('ERROR: This script is not untended to work with the production server')
	if len(args.api_key) != 43:
		sys.exit('ERROR: Invalid API key, please check it')
	else:
		api_key = args.api_key
	print()
	print('INFO: Working with server {}'.format(remote_addr))
	
	#get db connector and cursor
	db = get_db()
	curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
	#get local list of genes with no canonical isoform defined
	curs.execute(
		"DELETE FROM variant_feature WHERE c_name = 'c.1A>T'"
	)
	db.commit()
	
	http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
	
	curs.execute(
		"SELECT DISTINCT(name), nm_version FROM gene WHERE canonical = 't' ORDER BY name"
	)
	can = curs.fetchall()	
	i = 0
	j = 0
	variant = 'c.1A>T'
	failed_genes = []
	for gene in can:
		i += 1
		if i % 500 == 0:
			print('INFO: {0}/{1} genes checked)'.format(l, curs.rowcount))
		variant = '{0}.{1}:c.1A>T'.format(gene['name'][1], gene['nm_version'])
		md_url = '{0}/api/variant/create/{1}/{2}'.format(remote_addr, variant, api_key)
		try:
			md_response = json.loads(http.request('GET', md_url).data.decode('utf-8'))
			for key in md_response:
				if key == 'mobidetails_error':
					j += 1
					print('ERROR: variant creation failed for gene {0} with error {1}'.format(gene['name'], md_response[key]))
					failed_genes.append('{}-1'.format(gene['name'][0]))
		except:
			j += 1
			failed_genes.append('{}-2'.format(gene['name'][0]))
			print('ERROR: Call failed for gene {}'.format(gene['name']))
			continue
	print('{0}/{1} genes failed to create variants:'.format(j, curs.rowcount))
	print(failed_genes)
		
if __name__ == '__main__':
	main()