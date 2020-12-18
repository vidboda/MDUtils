import os
import re
import sys
import argparse
import psycopg2
import psycopg2.extras
import time
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config


def log(level, text):
    localtime = time.asctime( time.localtime(time.time()) )
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


def get_aliases(previous_symbols, current_aliases):
    aliases = '{0};{1}'.format(previous_symbols.replace(' ', ''), current_aliases.replace(' ', ''))
    aliases_obj = re.search(r'^;(.+)$', aliases)
    if aliases_obj:
        aliases = aliases_obj.group(1)
    aliases_obj = re.search(r'^(.+);$', aliases)
    if aliases_obj:
        aliases = aliases_obj.group(1)
    return aliases


def main():
    parser = argparse.ArgumentParser(description='Inset HGNC IDs in MD', usage='python add_hgnc_ids.py [-f path/to/hgnc_custom_file.txt/]')
    parser.add_argument('-f', '--file', default='', required=True, help='Path to the HGNC custom file')
    # HGNC file structure:
    # HGNC ID	Approved symbol	Approved name	Status	Previous symbols	Alias symbols	Chromosome	RefSeq IDs	Alias names
    args = parser.parse_args()
    # get sql file list
    if os.path.isfile(args.file):
        hgnc_file = args.file
    else:
        log('ERROR', 'Invalid input path, please check your command')
    i = 0
    j = 0
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    header = {
        'HGNC ID': 0,
        'Approved symbol': 1,
        'Approved name': 2,
        'Status': 3,
        'Previous symbols': 4,
        'Alias symbols': 5,
        'Chromosome': 6,
        'RefSeq IDs': 7,
        'Alias names': 8,
    }
    for line in open(hgnc_file).readlines():        
        i += 1
        if i == 1:
            continue
        hgnc = re.split("\t", line)
        hgnc_chr = None
        hgnc_id = None
        if hgnc[header['Status']] == 'Symbol Withdrawn' or \
                hgnc[header['Status']] == 'Entry Withdrawn':
            continue
        hgnc_chr_obj = re.search(r'^([\dXY]{1,2})[pq]*.*', hgnc[header['Chromosome']])
        if hgnc_chr_obj:
            hgnc_chr = hgnc_chr_obj.group(1)
        else:
            log('WARNING', 'Cannot identify Chr: line {0}: {1}'.format(i, line))
        hgnc_id_obj = re.search(r'^HGNC:(\d+)$', hgnc[header['HGNC ID']])
        if hgnc_id_obj:
            hgnc_id = hgnc_id_obj.group(1)
        if hgnc_id:
            curs.execute(
                "SELECT name, second_name, chr, hgnc_id FROM gene WHERE name[1] = %s AND canonical = 't'",
                (hgnc[header['Approved symbol']],)
            )
            md_gene = curs.fetchone()
            if md_gene:
                # check cannonical and print warning if differ
                # if md_gene['name'][1] != hgnc[header['RefSeq IDs']]:
                #     log('WARNING', 'RefSeq IDs discrepancy: HGNC: {0} - MD: {1}'.format(hgnc[header['RefSeq IDs']], md_gene['name'][1]))
                # check again good gene
                # log('DEBUG', 'MD gene: {0} - HGNC gene {1}'.format(md_gene['name'][0], hgnc[header['Approved symbol']]))
                if md_gene['name'][0] == hgnc[header['Approved symbol']] and \
                        int(md_gene['hgnc_id']) != int(hgnc_id):
                    if hgnc_chr and \
                            hgnc_chr == md_gene['chr']:
                        aliases = get_aliases(hgnc[header['Previous symbols']], hgnc[header['Alias symbols']])
                        # aliases = '{0}/{1}'.format(hgnc[header['Previous symbols']].replace(' ', ''), hgnc[header['Alias symbols']].replace(' ', ''))
                        # aliases_obj = re.search(r'^\/(.+)$', aliases)
                        # if aliases_obj:
                        #     aliases = aliases_obj.group(1)
                        # aliases_obj = re.search(r'^(.+)\/$', aliases)
                        # if aliases_obj:
                        #     aliases = aliases_obj.group(1)
                        # log('INFO', 'Aliases: MD: {0} - HGNC {1}-{2}'.format(md_gene['second_name'], hgnc[header['Previous symbols']], aliases))
                        log('INFO', 'Aliases: MD: {0} - HGNC {1}'.format(md_gene['second_name'], aliases))
                        curs.execute(
                            "UPDATE gene SET hgnc_id = %s WHERE name[1] = %s",
                            (hgnc_id, md_gene['name'][0])
                        )
                        j += 1
                        log('INFO', "UPDATE gene SET hgnc_id = {0} WHERE name[1] = '{1}'".format(hgnc_id, md_gene['name'][0]))
                        if aliases != ';' and \
                                aliases != md_gene['second_name']:
                            curs.execute(
                                "UPDATE gene SET second_name = %s WHERE name[1] = %s",
                                (aliases, md_gene['name'][0])
                            )
                            log('INFO', "UPDATE gene SET second_name = '{0}' WHERE name[1] = '{1}'".format(aliases, md_gene['name'][0]))
                        db.commit()
                    else:
                        log('WARNING', 'Chr discrepancy: MD chr: {0} - HGNC chr: {1}'.format(md_gene['chr'], hgnc_chr))
            else:
                previous_names_list = re.split(',', hgnc[header['Previous symbols']].replace(' ', ''))
                # update_semaph = 0
                for name in previous_names_list:
                    curs.execute(
                        "SELECT name, second_name, chr, hgnc_id FROM gene WHERE name[1] = %s AND canonical = 't'",
                        (name,)
                    )
                    md_gene_to_update = curs.fetchone()                    
                    if md_gene_to_update:
                        # update_semaph = 1
                        log('INFO', 'MD Gene name {0} needs to be updated to {1}'.format(md_gene_to_update['name'][0], hgnc[header['Approved symbol']]))
                        curs.execute(
                            "UPDATE gene SET name[1] = %s, hgnc_id = %s WHERE name[1] = %s",
                            (hgnc[header['Approved symbol']], hgnc_id, md_gene_to_update['name'][0])
                        )
                        db.commit()
                        j += 1
                        log('INFO', "UPDATE gene SET name[1] = '{0}', hgnc_id = {1} WHERE name[1] = '{2}'".format(hgnc[header['Approved symbol']], hgnc_id, md_gene_to_update['name'][0]))
                        aliases = get_aliases(hgnc[header['Previous symbols']], hgnc[header['Alias symbols']])
                        if aliases != ';' and \
                                aliases != md_gene_to_update['second_name']:
                            curs.execute(
                                "UPDATE gene SET second_name = %s WHERE name[1] = %s",
                                (aliases, hgnc[header['Approved symbol']])
                            )
                            log('INFO', "UPDATE gene SET second_name = '{0}' WHERE name[1] = '{1}'".format(aliases, hgnc[header['Approved symbol']]))
                            db.commit()
                        break
                    # if update_semaph == 0:
                    log('WARNING', 'HGNC not in MD: {0}-{1}'.format(hgnc[header['HGNC ID']], hgnc[header['Approved symbol']]))
        else:
            log('WARNING', 'Cannot identify HGNC ID: {}'.format(hgnc[header['HGNC ID']]))
        # if i > 50:
        #     break
    log('INFO', '{0} lines checked and {1} genes updated'.format(i, j))

if __name__ == '__main__':
    main()