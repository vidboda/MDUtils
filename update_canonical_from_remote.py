import sys
import os
import re
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import json
from insert_genes import get_db
import psycopg2.extras
#requires MobiDetails config module + database.ini file
from MobiDetailsApp import config

#fix genes in database did not have canonical transcripts.
#fix from remote MD server which has the information (typically dev server), using the API


def main():
	parser = argparse.ArgumentParser(description='Define a canonical transcript per gene', usage='python define_canonical.py [-r remote_server_url]')
	parser.add_argument('-r', '--remote-server', default='', required=True, help='base URL of the remote server')
	args = parser.parse_args()
	
	remote_addr = args.remote_server
	
	print()
	print('Working with server {}'.format(remote_addr))
	
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
				print(api_response[keys]['canonical'])
				if api_response[keys]['canonical'] == True:
					curs.execute(
						"UPDATE gene set canonical = 't' WHERE name[2] = '{}'".format(keys)
					)
	db.commit()
	print("{} genes modified".format(i))		
	
if __name__ == '__main__':
	main()