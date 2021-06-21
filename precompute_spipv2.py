
import os
import sys
import time
import tempfile
import subprocess
import psycopg2
import psycopg2.extras

# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config, md_utilities

# script to precompute existing variants


def log(level, text):
    localtime = time.asctime(time.localtime(time.time()))
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


def get_db():
    try:
        # read connection parameters
        params = config.mdconfig()
        db = psycopg2.connect(**params)
    except (Exception, psycopg2.DatabaseError) as error:
        log('ERROR', error)
    return db


def main():
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    curs.execute(
        "SELECT id, gene_name[1] as symbol, gene_name[2] as nm, c_name FROM variant_feature ORDER BY RANDOM() LIMIT 3000"
    )
    res = curs.fetchall()
    total = curs.rowcount
    precomputed = 0
    j = 0
    for var in res:
        j += 1
        if j % 500 == 0:
            log('INFO', '{0}/{1} variant checked'.format(j, total))
        if not os.path.exists(
                '{0}{1}.txt'.format(md_utilities.local_files['spip']['abs_path'], var['id'])
                ):
            tf = tempfile.NamedTemporaryFile(suffix='.txt')
            # tfout = tempfile.NamedTemporaryFile()
            tf.write(
                bytes('gene\tvarID\n{0}\t{1}:c.{2}'.format(
                    var['symbol'], var['nm'], var['c_name']
                ), encoding='utf-8')
            )
            # tf.write(b"gene    varID\nUSH2A  NM_206933:c.2276G>T")
            tf.seek(0)
            result = subprocess.run(
                [
                    md_utilities.ext_exe['Rscript'],
                    '{}'.format(md_utilities.ext_exe['spip']),
                    '-I',
                    '{}'.format(tf.name),
                    '-O',
                    '{0}{1}.txt'.format(md_utilities.local_files['spip']['abs_path'], var['id']),
                    '-g',
                    'hg38'
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT
            )
            if result.returncode == 0:
                with open(r'{0}{1}.txt'.format(md_utilities.local_files['spip']['abs_path'], var['id'])) as spip_file:
                    num_lines = len(spip_file.readlines())
                if num_lines != 2:
                    # we don't keep if the number of lines is wrong
                    os.remove(r'{0}{1}.txt'.format(md_utilities.local_files['spip']['abs_path'], var['id']))
                    continue
                precomputed += 1
    log('INFO', 'Pre-computed {0} variants / {1}'.format(precomputed, total))


if __name__ == '__main__':
    main()
