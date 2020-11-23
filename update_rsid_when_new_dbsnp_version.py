import os
import re
import argparse
import psycopg2
import psycopg2.extras
import time
import tabix
from insert_genes import get_db
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config

# script meant to be ran manually after a change in dbSNP version
# will update rsids for old variants w/ no rsids - should be ran on prod server


def log(level, text):
    localtime = time.asctime( time.localtime(time.time()) )
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


def main():
    # script meant to update variants when the canonical form of a gene is changed
    parser = argparse.ArgumentParser(description='Define a canonical transcript per gene when several are defined',
                                     usage='python update_canonical_when_several.py -k md_api_key -g gene_hgnc')
    parser.add_argument('-d', '--dbsnp-file', default='', required=True,
                        help='Path to the dbSNP file')
    args = parser.parse_args()
    if os.path.isfile(args.dbsnp_file):
        db = get_db()
        curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
        curs.execute(
            "SELECT a.id, b.pos, b.pos_ref, b.pos_alt, c.ncbi_name FROM variant_feature a, variant b, chromosomes c WHERE a.id = b.feature_id AND b.chr = c.name AND c.genome_version = 'hg38' AND b.genome_version = 'hg38' ORDER BY a.id"
        ) #  AND a.dbsnp_id is NULL
        res = curs.fetchall()
        count = curs.rowcount
        log('INFO', 'Found {0} variants to check'.format(count))
        i = 0
        j = 0
        for var in res:
            j += 1
            if j % 500 == 0:
                log('INFO', '{0}/{1} variant checked'.format(j, count))
            tb = tabix.open(args.dbsnp_file)
            query = "{0}:{1}-{2}".format(var['ncbi_name'], var['pos'], var['pos'])
            records = tb.querys(query)
            for record in records:
                match_object = re.search(r'RS=(\d+);', record[7])
                if match_object:
                    pos_ref_list = re.split(',', record[3])
                    pos_alt_list = re.split(',', record[4])
                    if var['pos_ref'] in pos_ref_list and \
                            var['pos_alt'] in pos_alt_list:
                        new_rs_id = match_object.group(1)
                        # need to update var entry
                        curs.execute(
                            "UPDATE variant_feature SET dbsnp_id = %s WHERE id = %s",
                            (new_rs_id, var['id'])
                        )
                        log('INFO', 'Adding rsid for variant {0} - ref:{1} - alt:{2} to rs{3}'.format(var['id'], var['pos_ref'], var['pos_alt'], new_rs_id))
                        i += 1
                        db.commit()
        log('INFO', '{0} rsids added'.format(i))
    else:
        log('ERROR', 'Your input file was not found {0}'.format(args.dbsnp_file))


if __name__ == '__main__':
    main()