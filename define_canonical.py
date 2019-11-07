import sys
import os
import re
import argparse
import psycopg2
import psycopg2.extras
from insert_genes import get_db
import psycopg2.extras
#requires MobiDetails config module + database.ini file
from MobiDetailsApp import config

#first genes in databse did not have cnaonical transcripts.
#fixed with the refGenecanonical file

def main():
	parser = argparse.ArgumentParser(description='Define a canonical transcript per gene', usage='python define_canonical.py [-r path/to/refGeneCanonical_2019_09_23.txt]')
	parser.add_argument('-r', '--refgene', default='', required=True, help='Path to the file containing the canonical refSeq IDs per gene')
	args = parser.parse_args()
	#get file 
	if os.path.isfile(args.refgene):
		refgeneFile = args.refgene
	else:
		sys.exit('Invalid input path, please check your command')
	
	#get db connector and cursor
	db = get_db()
	curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
	i = 0
	lacking_nm = []

	for geneLine in open(refgeneFile).readlines():
		#ENST - NM - gene
		geneLineList = geneLine.rstrip().split("\t")
		#print(geneLineList[2])
		if geneLineList[2] != 'n/a' and geneLineList[2] != 'hg38.refGene.name2':
			curs.execute(#gene exists in MD (no main already set)
				"SELECT DISTINCT(name[1]) FROM gene WHERE name[1] = '{}' AND name[1] NOT IN (SELECT name[1] FROM gene WHERE canonical = 't');".format(geneLineList[2])
			)
			mdgene = curs.fetchone()
			
			if mdgene is not None:
				#nm exists in md?
				curs.execute(
					"SELECT name FROM gene WHERE name[2] = '{0}'".format(geneLineList[1])
				)#exists in table gene_annotation? get a nm
				mdnm = curs.fetchone()
				if mdnm is not None:
					#ok => canonical
					i += 1
					postGene = '{"' + mdnm['name'][0] + '","' + mdnm['name'][1] + '"}'
					#print("UPDATE gene SET canonical = 't' WHERE name = '{}'".format(postGene))
					curs.execute(
					 	"UPDATE gene SET canonical = 't' WHERE name = '{}'".format(postGene)
					)
				else:
					lacking_nm.append(geneLineList[2])
	print(lacking_nm)
	print("Genes modified {}".format(i))			
	
	db.commit()
	
	
		
if __name__ == '__main__':
	main()