import re
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import json
from insert_genes import get_db
#requires MobiDetails config module + database.ini file
from MobiDetailsApp import config

#fix genes in database did not have canonical transcripts.
#fix from remote MD server which has the information (typically dev server), using the API
#fix UNIPORT IDs, using the MDAPI and Uniprot API


def main():
	parser = argparse.ArgumentParser(description='Define a canonical transcript per gene and optionally updates various fields', usage='python update_canonical_from_remote.py [-r remote_server_url]')
	parser.add_argument('-r', '--remote-server', default='', required=True, help='base URL of the remote server')
	parser.add_argument('-np', '--update-np', default='', required=False, help='Optionally update NP for genes', action='store_true')
	parser.add_argument('-uu', '--update-uniprot', default='', required=False, help='Optionally update UNIPROT IDs', action='store_true')
	parser.add_argument('-uc', '--update-creation', default='', required=False, help='Optionally update variant_creation tag', action='store_true')
	parser.add_argument('-un', '--update-nm', default='', required=False, help='Optionally update RefSeq Nm accession number tag', action='store_true')
	
	args = parser.parse_args()
	remote_addr = args.remote_server
#	args = parser.parse_args(['-np'])
	print()
	print('INFO: Working with server {}'.format(remote_addr))
	
	#get db connector and cursor
	db = get_db()
	curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
	#get local list of genes with no canonical isoform defined
	curs.execute(
		"SELECT DISTINCT(name[1]) as hgnc FROM gene WHERE name[1] NOT IN (SELECT name[1] FROM gene WHERE canonical='t') ORDER BY name[1]"
	)
	no_can = curs.fetchall()
	
	http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
	
	i = 0
	lacking_nm = []
	
	for gene in no_can:
		req_url = '{0}/api/gene/{1}'.format(remote_addr, gene['hgnc'])
		api_response = json.loads(http.request('GET', req_url).data.decode('utf-8'))
		for keys in api_response:
			if 'canonical' in api_response[keys]:
				if api_response[keys]['canonical'] is True:
					if re.search(r'NM_\d+\.\d', keys):
						match_obj = re.search(r'(NM_\d+)\.\d', keys)
						nm_acc = match_obj.group(1)
						curs.execute(
							"UPDATE gene set canonical = 't' WHERE name[2] = '{}'".format(nm_acc)
						)
						print("INFO: Updating {}".format(nm_acc))
						i += 1
	db.commit()
	print("INFO: {} genes modified".format(i))
	
	if args.update_np:
		curs.execute(
			"SELECT DISTINCT(name[1]) as hgnc FROM gene WHERE np = 'NP_000000.0'"
		)
		no_np = curs.fetchall()
		j = 0
		for gene in no_np:
			req_url = '{0}/api/gene/{1}'.format(remote_addr, gene['hgnc'])
			api_response = json.loads(http.request('GET', req_url).data.decode('utf-8'))
			for keys in api_response:
				if 'RefProtein' in api_response[keys] and api_response[keys]['RefProtein'] != 'NP_000000.0':
					if re.search(r'NP_\d+\.\d', api_response[keys]['RefProtein']):
						match_obj = re.search(r'(NM_\d+)\.\d', keys)
						nm_acc = match_obj.group(1)
						np_acc = api_response[keys]['RefProtein']
						curs.execute(
							"UPDATE gene set np = '{0}' WHERE name[2] = '{1}'".format(np_acc, nm_acc)
						)
						print("INFO: Updating gene NP acc no of {0} to {1}".format(nm_acc, np_acc))
						j += 1
		db.commit()
		print("INFO: {} NP acc no modified".format(j))
	if args.update_uniprot or args.update_creation or args.update_nm:
		curs.execute(
			"SELECT  name[1] as HGNC, name[2] as nm, nm_version, np, uniprot_id, variant_creation FROM gene ORDER BY name"
		)
		res = curs.fetchall()
		k = l = m = n = 0
		#l = 0
		#m = 0
		#n = 0
		o = curs.rowcount
		for gene in res:
			req_url = '{0}/api/gene/{1}'.format(remote_addr, gene['hgnc'])
			api_response = json.loads(http.request('GET', req_url).data.decode('utf-8'))
			l += 1
			if l % 1000 == 0:
				print('INFO: {0}/{1} isoforms checked)'.format(l, o))
			for keys in api_response:
				match_obj = re.search(r'(NM_\d+)\.(\d)', keys)
				if match_obj:
					nm_acc = match_obj.group(1)
					#check again
					if nm_acc == gene['nm']:
						if args.update_nm:
							nm_version = match_obj.group(2)
							if int(nm_version) > gene['nm_version']:
								#no downgrade?
								curs.execute(
									"UPDATE gene set nm_version = '{0}' WHERE name[2] = '{1}'".format(nm_version, nm_acc)
								)
								print("INFO: Updating gene RefSeq NM accession version of {0} to {1}".format(nm_acc, nm_version))
								n += 1
						if 'UNIPROT' in api_response[keys] and args.update_uniprot:					
							uniprot = api_response[keys]['UNIPROT']
							if uniprot != gene['uniprot_id']:
								curs.execute(
									"UPDATE gene set uniprot_id = '{0}' WHERE name[2] = '{1}'".format(uniprot, nm_acc)
								)
								print("INFO: Updating gene UNIPROT id of {0} to {1}".format(nm_acc, uniprot))
								k += 1
						if 'variantCreationTag' in api_response[keys] and args.update_creation:
							tag = api_response[keys]['variantCreationTag']
							if tag != gene['variant_creation']:
								curs.execute(
									"UPDATE gene set variant_creation = '{0}' WHERE name[2] = '{1}'".format(tag, nm_acc)
								)
								print("INFO: Updating gene variantCreationTag of {0} to {1}".format(nm_acc, tag))
								m += 1
		db.commit()
		print("INFO: {} UNIPROT IDs modified".format(k))
		print("INFO: {} variantCreationTag modified".format(m))
		print("INFO: {} RefSeq NM accession version modified".format(n))
	
if __name__ == '__main__':
	main()