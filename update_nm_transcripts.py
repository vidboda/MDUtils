import re
import urllib3
import certifi
import json
import psycopg2
import psycopg2.extras

from insert_genes import get_db
#requires MobiDetails config module + database.ini file
from MobiDetailsApp import config

def main():
	#script meant to be croned to update NM acc versions in MD according to VariantValidator
	#to be ran after uta docker update for example
	#uses VV API genes2transcript
	#https://rest.variantvalidator.org:443/tools/gene2transcripts/{gene_name}
	#vv_url_base = "https://rest.variantvalidator.org:443"
	vv_url_base = "http://0.0.0.0:8000/"
	
	db = get_db()
	curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
	
	curs.execute(#get genes
		"SELECT name, nm_version FROM gene"
	)
	genes = curs.fetchall()
	i = 0
	for gene in genes:
		print("{}-{}".format(gene['name'][0], i))
		i += 1
		#print("MD------{}".format(gene['name'][1]))
		#get VV info for the gene
		http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
		vv_url = "{0}/tools/gene2transcripts/{1}".format(vv_url_base, gene['name'][1])
		#return intervar_url
		try:
			vv_data = json.loads(http.request('GET', vv_url).data.decode('utf-8'))
			if 'transcripts' in vv_data:
				for transcript in vv_data['transcripts']:
					#print("VV------{}".format(transcript['reference']))
					match_object = re.search('^(N[MR]_\d+)\.(\d{1,2})', transcript['reference'])
					if match_object is not None:
						nm_acc = match_object.group(1)
						#if nm_acc == gene['name'][1]:
						nm_version = match_object.group(2)
						if nm_acc == gene['name'][1] and int(nm_version) > int(gene['nm_version']):
							#print("IN--{}-{}-{}-{}-{}-".format(vv_data['current_symbol'], transcript['reference'], gene['name'][1], nm_version, gene['nm_version']))
							curs.execute(
								"UPDATE gene SET nm_version = '{0}' WHERE name[2] = '{1}'".format(nm_version, gene['name'][1])
							)
							print("UPDATE gene SET nm_version = '{0}' WHERE name[2] = '{1}'".format(nm_version, gene['name'][1]))
		except:
			print('No value for {0}'.format(gene['name'][0]))
	
	db.commit()
	
	
		
if __name__ == '__main__':
	main()