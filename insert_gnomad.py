import sys
import os
import re
import argparse
import psycopg2
import psycopg2.extras
from insert_genes import get_db
#requires MobiDetails config module + database.ini file
from MobiDetailsApp import config


def main():
	parser = argparse.ArgumentParser(description='Insert gnomAD data into MD', usage='python insert_gnomad.py [-d path/to/dir/containing/gnomad.v2.1.1.lof_metrics.by_gene.txt]')
	parser.add_argument('-d', '--directory', default='', required=True, help='Path to the directory containing the gnomAD metrics by gene file')
	args = parser.parse_args()
	#get file 
	if os.path.isfile(args.directory):
		gnomadFile = args.directory
	else:
		sys.exit('Invalid input path, please check your command')
	
	#get db connector and cursor
	db = get_db()
	curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
	i = 0

	for geneLine in open(gnomadFile).readlines():
		geneLineList = geneLine.split("\t")
		#print(geneLineList[0])
		curs.execute(#exists in MD?
			"SELECT name FROM gene WHERE name[1] = '{0}' AND number_of_exons IN (SELECT MAX(number_of_exons) FROM gene WHERE name[1] = '{0}')".format(geneLineList[0])
		)
		mdNMFirst = curs.fetchone()
		if mdNMFirst is not None:
			#print(mdNMFirst['nm'])
			curs.execute(
				"SELECT DISTINCT(gene_name[2]) FROM gene_annotation WHERE gene_name[1] = '{}'".format(geneLineList[0])
			)#exists in table gene_annotation? get a nm
			mdNMSecond = curs.fetchone()
			if mdNMSecond is None:
				#does not exists => creation
				i += 1
				postGene = '{"' + mdNMFirst['name'][0] + '","' + mdNMFirst['name'][1] + '"}'
				#print(postGene,geneLineList[13],geneLineList[24],geneLineList[25],geneLineList[4],geneLineList[26],geneLineList[27],geneLineList[23],geneLineList[28],geneLineList[29])
				oeValues = {
					'synoe': geneLineList[13], 'synlower': geneLineList[24], 'synupper': geneLineList[25], 'misoe': geneLineList[4], 'mislower': geneLineList[26], 'misupper': geneLineList[27], 'lofoe': geneLineList[23], 'loflower': geneLineList[28], 'lofupper': geneLineList[29]
				}
				for oeval in oeValues:
					try:
						oeValues[oeval] = float(oeValues[oeval])
						oeValues[oeval] = "{:.2f}".format(oeValues[oeval])
					except:
						next

				#print("INSERT INTO gene_annotation VALUES('{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}')".format(postGene,oeValues['synoe'],oeValues['synlower'],oeValues['synupper'],oeValues['misoe'],oeValues['mislower'],oeValues['misupper'],oeValues['lofoe'],oeValues['loflower'],oeValues['lofupper']))
				curs.execute(
				 	"INSERT INTO gene_annotation VALUES('{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}')".format(postGene,oeValues['synoe'],oeValues['synlower'],oeValues['synupper'],oeValues['misoe'],oeValues['mislower'],oeValues['misupper'],oeValues['lofoe'],oeValues['loflower'],oeValues['lofupper'])
				)

	print(i)			
	
	db.commit()
	
	
		
if __name__ == '__main__':
	main()