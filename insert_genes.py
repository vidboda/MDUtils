import sys
#from os import listdir
#from os.path import isfile, join
import os
import re
import argparse
import psycopg2
import psycopg2.extras
#requires MobiDetails config module + database.ini file
from MobiDetailsApp import config



def get_db():
	try:
		# read connection parameters
		params = config.mdconfig()
		db = psycopg2.connect(**params)
	except (Exception, psycopg2.DatabaseError) as error:
		sys.exit(error)
	return db

def get_file_list(sqlPath):
	sqlFiles = []
	for f in os.listdir(sqlPath):
		if os.path.isfile(os.path.join(sqlPath, f)) and f.lower().endswith('.sql'):
			tmpFile = [f]
			sqlFiles.extend(tmpFile)
	if sqlFiles:
		return(sqlFiles)
	else:
		sys.exit('No SQL files in path {}'.format(sqlPath))

def main():
	parser = argparse.ArgumentParser(description='Insert MD SQL gene files', usage='python insert_genes.py [-d path/to/dir/containing/md/sql/genes/files/]')
	parser.add_argument('-d', '--directory', default='', required=True, help='Path to the directory containing the SQL genes files for MobiDetails')
	args = parser.parse_args()
	#get sql file list
	if os.path.isdir(args.directory):
		sqlPath = args.directory
	else:
		sys.exit('Invalid input path, please check your command')
	sqlFiles = get_file_list(sqlPath)
	#print(sqlFiles)
	
	#get db connector and cursor
	db = get_db()
	curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
	i = 0
	for sqlFile in sorted(sqlFiles):
		i += 1
		for line in open(os.path.join(sqlPath, sqlFile)).readlines():
			if re.match('INSERT INTO (gene|segment|protein_domain)', line):
				try:
					line = re.sub(",'',", ",NULL,", line)
					line = re.sub(",'NULL',", ",NULL,", line)
					curs.execute(line)
				except psycopg2.Error as e:
					print(sqlFile + "\n" + line + "\n" + e)
			else:
				print('non matching line ' + line + 'in file ' + sqlFile)
		
		print(sqlFile + ' inserted - #' + str(i))
		
	db.commit()
	
	
		
if __name__ == '__main__':
	main()