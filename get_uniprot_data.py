import os
import sys
import re
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import argparse
from precompute_spipv2 import get_db, log
from MobiDetailsApp import md_utilities

# get UNIPROT GFF file and populate uniprot_domain Table


def main():
    parser = argparse.ArgumentParser(description='Manage UNIPROT data', usage='python get_uniprot_data.py [-f path/to/dir/containing/uniprot/ids.txt]')
    parser.add_argument('-f', '--file', default='', required=False, help='Path to the UNIPROT ID file to be added to MD')
    args = parser.parse_args()
    # get a file of uniprot id list as input, otherwise queries the whole database
    id_file = None
    id_list = []
    if args.file and \
            os.path.isfile(args.file):
        # log('DEBUG', 'File: {0}'.format(args.file))
        id_file = args.file
        for line in open(id_file).readlines():
            if re.search(r'^\w+$', line):
                id_list.append(line.rstrip())
    elif args.file:
        log('ERROR', 'Invalid input path for gene file, please check your command')
    # get db connector and cursor
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    if not id_file:
        curs.execute(
            "SELECT id FROM uniprot"
        )
        res_ids = curs.fetchall()
        for id in res_ids:
            id_list.append(id['id'])
    header = {
        'User-Agent': 'python-requests Python/{}.{}.{} - MDUtils'.format(sys.version_info[0], sys.version_info[1], sys.version_info[2]),
    }
    for id in id_list:
        # exists in MD?
        curs.execute(
            "SELECT id FROM uniprot WHERE id = %s",
            (id,)
        )
        res_id = curs.fetchone()
        if not res_id:
            # insetr value
            curs.execute(
                "INSERT INTO uniprot (id) VALUES (%s)",
                (id,)
            )
            db.commit()
        uniprot_url = 'https://www.uniprot.org/uniprot/{0}.gff'.format(id)
        # log('DEBUG', 'URL: {0}'.format(uniprot_url))
        try:
            uniprot_response = http.request('GET', uniprot_url, headers=header).data.decode('utf-8')
        except Exception:
            log('WARNING', 'No value for {0}'.format(id))
            continue
        if not os.path.isfile(
            '{0}{1}.json'.format(
                md_utilities.local_files['uniprot']['abs_path'],
                id
            )
                ):
            # copy in file system
            gff_file = open(
                '{0}{1}.gff'.format(
                    md_utilities.local_files['uniprot']['abs_path'],
                    id
                ),
                "w",
                encoding='utf-8'
            )
            gff_file.write(uniprot_response)
        for line in open('{0}{1}.gff'.format(
            md_utilities.local_files['uniprot']['abs_path'],
            id
        )).readlines():
            # get domain info
            if re.search(rf'^{id}\t', line):
                info = re.split('\t', line)
                # log('DEBUG', 'info: {0}'.format(info))
                if info[2] == 'Domain':
                    name = re.split('=', re.split(';', info[8])[0])[1]
                    # exists?
                    curs.execute(
                        "SELECT name FORM uniprot_domain WHERE uniprot_id = %s AND name = %s AND aa_start = %s",
                        (
                            id,
                            name,
                            info[3]
                        )
                    )
                    res_exists = curs.fetchone()
                    if not res_exists:
                        # log('DEBUG', "INSERT INTO uniprot_domain (uniprot_id, aa_start, aa_end, name) VALUES ('{0}', '{1}', '{2}', '{3}')".format(
                        #     id,
                        #     info[3],
                        #     info[4],
                        #     name
                        # ))
                        curs.execute(
                            "INSERT INTO uniprot_domain (uniprot_id, aa_start, aa_end, name) VALUES (%s, %s, %s, %s)",
                            (
                                id,
                                info[3],
                                info[4],
                                name
                            )
                        )
                        db.commit()


if __name__ == '__main__':
    main()
