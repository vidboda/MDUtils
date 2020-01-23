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

def log(level, text):
	print()
	if level == 'ERROR':
		sys.exit('[{0}]: {1}'.format(level, text))
	print('[{0}]: {1}'.format(level, text))

def main():
	parser = argparse.ArgumentParser(description='Checks that genes accept variant creation', usage='python check_variant_creation.py [-r remote_server_url]')
	parser.add_argument('-r', '--remote-server', default='', required=True, help='base URL of the remote server')
	parser.add_argument('-k', '--api-key', default='', required=True, help='Your API key visible on your profile page on the website.')
	
	args = parser.parse_args()
	remote_addr = args.remote_server
	if re.search('mobidetails\.iurc', remote_addr):
		log('ERROR', 'This script is not untended to work with the production server')
	if len(args.api_key) != 43:
		log('ERROR', 'Invalid API key, please check it')
	else:
		api_key = args.api_key
	print()
	log('INFO','Working with server {}'.format(remote_addr))
	
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
		"SELECT DISTINCT(name), nm_version, variant_creation FROM gene WHERE canonical = 't' ORDER BY name"
	)
	can = curs.fetchall()
	num_can = curs.rowcount
	i = 0
	j = 0
	k = 0
	variant = 'c.1A>T'
	failed_genes = []
	for gene in can:
		print('.', end="", flush=True)
		i += 1
		if i % 500 == 0:
			log('INFO', '{0}/{1} genes checked'.format(i, num_can))
		variant = '{0}.{1}:c.1A>T'.format(gene['name'][1], gene['nm_version'])
		md_url = '{0}/api/variant/create/{1}/{2}'.format(remote_addr, variant, api_key)
		try:
			md_response = json.loads(http.request('GET', md_url).data.decode('utf-8'))			
			if 'mobidetails_error' in md_response:
				j += 1
				log('WARNING', 'variant creation failed for gene {0} with error {1}'.format(gene['name'], md_response['mobidetails_error']))
				new_nm_match_obj = re.search(r'A more recent version of the selected reference sequence NM_\d+\.\d is available \((NM_\d+)\.(\d)\)', md_response['mobidetails_error'])
				if new_nm_match_obj:
					nm_to_check = new_nm_match_obj.group(1)
					new_ver = new_nm_match_obj.group(2)
					if nm_to_check == gene['name'][1]:
						curs.execute(
							"UPDATE gene SET nm_version = '{0}' WHERE name[2] = '{}'".format(new_ver, gene['name'][1])
						)
					#recheck
					variant_2 = '{0}.{1}:c.1A>T'.format(gene['name'][1], new_ver)
					md_url_2 = '{0}/api/variant/create/{1}/{2}'.format(remote_addr, variant_2, api_key)
					try:
						md_response_2 = json.loads(http.request('GET', md_url_2).data.decode('utf-8'))
						if 'mobidetails_id' in md_response_2 and gene['variant_creation'] != 'ok':
							curs.execute(
								"UPDATE gene SET variant_creation = 'ok' WHERE name[2] = '{}'".format(gene['name'][1])
							)
							continue
					except:
						k += 1
						failed_genes.append('{}'.format(gene['name'][0]))
						continue
				if re.search(r'cannot be mapped directly to genome build GRCh38', md_response['mobidetails_error']):
					curs.execute(
						"UPDATE gene SET variant_creation = 'hg38_mapping_default' WHERE name[2] = '{}'".format(gene['name'][1])
					)
					log('INFO','MD gene table updated with variant_creation = hg38_mapping_default')
				elif re.search(r'does not seem to map correctly to hg19', md_response['mobidetails_error']):
					curs.execute(
						"UPDATE gene SET variant_creation = 'hg19_mapping_default' WHERE name[2] = '{}'".format(gene['name'][1])
					)
					log('INFO','MD gene table updated with variant_creation = hg19_mapping_default')
			elif 'mobidetails_id' in md_response and gene['variant_creation'] != 'ok':
				curs.execute(
					"UPDATE gene SET variant_creation = 'ok' WHERE name[2] = '{}'".format(gene['name'][1])
				)
			db.commit()
		except:
			k += 1
			failed_genes.append('{}'.format(gene['name'][0]))
			#### TODO  put gene variant creation tag as f
			#print('ERROR: Call failed for gene {}'.format(gene['name']))
			continue
	log('INFO','{0}/{1} genes reported a VV error'.format(j, num_can))
	log('INFO','{0}/{1} genes triggered an MD error:'.format(k, num_can))
	log('INFO',failed_genes)
		
if __name__ == '__main__':
	main()