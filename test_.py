#import pytest
from . import insert_genes

#test db connection

def test_db_conn():
	test_value = False
	try:
		db = insert_genes.get_db()
		db.close()
		test_value = True
	except:
		test_value = False
	assert test_value == True
	

def test_file_list():
	files = insert_genes.get_file_list('test_gene')
	assert files[0] == 'XPC_SQL.sql'