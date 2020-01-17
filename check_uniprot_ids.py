import sys
import os
import re
#import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import json
from insert_genes import get_db
import psycopg2.extras
#requires MobiDetails config module + database.ini file
from MobiDetailsApp import config


#check UNIPROT IDs

def main():
	#get db connector and cursor
	db = get_db()
	curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
	
	i = 0
	
	http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
	curs.execute(
		"SELECT name, np, uniprot_id FROM gene ORDER BY name"
	)
	res = curs.fetchall();
	for gene in res:
		#ncbi
		#print('INFO: processing RefSeq: {0} - {1} - {2}'.format(gene['np'], gene['name'][1], gene['name'][0]))
		uniprot_url = 'https://www.ebi.ac.uk/proteins/api/proteins/refseq:{}?offset=0&size=100&reviewed=true'.format(gene['np'])
		uniprot_response = json.loads(http.request('GET', uniprot_url, headers={'Accept':'application/json'}).data.decode('utf-8'))
		#print(uniprot_response[0]['accession'])
		try:
			if uniprot_response[0]['accession']:
				if gene['uniprot_id'] == uniprot_response[0]['accession']:
					#print('INFO: RefSeq: {0} - {1} - {2} OK'.format(gene['np'], gene['name'][1], gene['name'][0]))
					continue
				else:
					#print('ERROR: md_uniprot_id: {0} - ncbi uniprot_id: {1} - RefSeq: {2} - {3} - {4} : discrepancy'.format(gene['uniprot_id'],uniprot_response[0]['accession'], gene['np'], gene['name'][1], gene['name'][0]))
					curs.execute(
						"UPDATE gene SET uniprot_id = '{0}' WHERE name[2] = '{1}'".format(uniprot_response[0]['accession'], gene['name'][1])
					)
					#print("UPDATE gene SET uniprot_id = '{0}' WHERE name[2] = '{1}'".format(uniprot_response[0]['accession'], gene['name'][1]))
					print('WARNING Updated gene UNIPROT ID of {0} - {1} from {2} to {3}'.format(gene['name'][0], gene['name'][1], gene['uniprot_id'], uniprot_response[0]['accession']))
					i += 1
			else:
				print('WARNING: md_uniprot_id: {0} - RefSeq: {1} - {2} - {3} :not checked'.format(gene['uniprot_id'], gene['np'], gene['name'][1], gene['name'][0]))
		except:
			print('ERROR: no UNIPROT ID {0} for {1} - {2}'.format(uniprot_response, gene['name'][1], gene['name'][0]))
	print("{} isoforms updated".format(i))
	
	db.commit()
	
		
if __name__ == '__main__':
	main()