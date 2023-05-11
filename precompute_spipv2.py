import os
import sys
import time
import tempfile
import subprocess
import psycopg2
import psycopg2.extras

# requires MobiDetails config module + database.ini file
# from MobiDetailsApp import config, md_utilities, get_db
from MobiDetailsApp import md_utilities, configuration

# script to precompute existing variants


def log(level, text):
    localtime = time.asctime(time.localtime(time.time()))
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('\n[{0}]: {1} - {2}'.format(level, localtime, text))


def get_db():
    try:
        # read connection parameters
        params = configuration.mdconfig()
        db_pool = psycopg2.pool.ThreadedConnectionPool(
            1,
            5,
            user=params['user'],
            password=params['password'],
            host=params['host'],
            port=params['port'],
            database=params['database']
        )
        # db = psycopg2.connect(**params)
    except (Exception, psycopg2.DatabaseError) as error:
        log('ERROR', error)
    return db_pool, db_pool.getconn()


def main():
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    curs.execute(
        """
        SELECT id, gene_symbol, refseq, c_name
        FROM variant_feature
        ORDER BY RANDOM()
        LIMIT 2
        """
    )
    res = curs.fetchall()
    total = curs.rowcount
    precomputed = 0
    j = 0
    for var in res:
        j += 1
        # log('DEBUG', j)
        if j % 500 == 0:
            log('INFO', '{0}/{1} variant checked'.format(j, total))
        if not os.path.exists(
                '{0}{1}.txt'.format(md_utilities.local_files['spip']['abs_path'], var['id'])
                ):
            tf = tempfile.NamedTemporaryFile(suffix='.txt')
            # tfout = tempfile.NamedTemporaryFile()
            tf.write(
                bytes('gene\tvarID\n{0}\t{1}:c.{2}'.format(
                    var['gene_symbol'], var['refseq'], var['c_name']
                ), encoding='utf-8')
            )
            # tf.write(b"gene    varID\nUSH2A  NM_206933:c.2276G>T")
            tf.seek(0)
            result = subprocess.run(
                [
                    md_utilities.ext_exe['Rscript'],
                    f"{md_utilities.ext_exe['spip']}",
                    '-I',
                    f'{tf.name}',
                    '-O',
                    '{0}{1}.txt'.format(
                        md_utilities.local_files['spip']['abs_path'], var['id']
                    ),
                    '-g',
                    'hg38',
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
            # log('DEBUG', result.returncode)
            if result.returncode == 0:
                with open(r'{0}{1}.txt'.format(md_utilities.local_files['spip']['abs_path'], var['id'])) as spip_file:
                    num_lines = len(spip_file.readlines())
                if num_lines != 2:
                    # we don't keep if the number of lines is wrong
                    os.remove(r'{0}{1}.txt'.format(md_utilities.local_files['spip']['abs_path'], var['id']))
                    continue
                precomputed += 1
    log('INFO', 'Pre-computed {0} variants / {1}'.format(precomputed, total))
    # db.close()
    db_pool.putconn(db)


if __name__ == '__main__':
    main()
