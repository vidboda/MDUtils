import re
import json
import psycopg2
import psycopg2.extras
from precompute_spipv2 import get_db, log
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import md_utilities
from check_transcripts_in_vvjson import download_vv_file


def main():
    # script to add genes in MD from an HGNC full set download
    # vv_url_base = "https://www608.lamp.le.ac.uk"
    ncbi_chr = {}
    ncbi_chr_hg19 = {}
    # hgnc file in resources/hgnc/hgnc_complete_set.txt
    # beginning of line example
    # HGNC:24086      A1CF    APOBEC1 complementation factor  protein-coding gene     gene with protein product       Approved
    # so we need id [0].split(':')[1], symbol [1], and validation [3] = protein-coding gene and [5] = Approved
    hgnc_file = open(md_utilities.local_files['hgnc_full_set']['abs_path'], 'r')
    i = 0
    for line in hgnc_file:
        i += 1
        if i % 500 == 0:
            log('INFO', '{0}/{1} lines checked'.format(i, 42946))
        gene_info = re.split('\t', line)
        match_id = re.search(r'HGNC:(\d+)', gene_info[0])
        if match_id:
            hgnc_id = match_id.group(1)
            # log('DEBUG', 'Status:{0}'.format(gene_info[4]))
            if gene_info[3] == 'protein-coding gene' and \
                    gene_info[5] == 'Approved':
                # log('DEBUG', 'HGNC:{0}'.format(hgnc_id))
                # suitable for MD
                # check ID VS MD
                # if ok, check current symbol and change it if needed
                # else download file from VV and insert new gene and transcript
                hgnc_current_symbol = gene_info[1]
                db = get_db()
                curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
                curs.execute(  # get genes - one transcript per gene (canonical) - allows update of all trasncripts
                    """
                    SELECT name, hgnc_name, second_name
                    FROM gene
                    WHERE hgnc_id = %s
                    """,
                    (hgnc_id,)
                )
                md_gene = curs.fetchone()
                if md_gene:
                    # check symbol and name
                    if md_gene['name'][0] != hgnc_current_symbol:
                        # log('DEBUG', 'Symbol differs for HGNC {0}, MD {1}'.format(hgnc_current_symbol, md_gene['name'][0]))
                        # get the one VV is using
                        download_vv_file(hgnc_current_symbol, hgnc_current_symbol)
                        vv_file = open('{0}{1}.json'.format(
                            md_utilities.local_files['variant_validator']['abs_path'],
                            hgnc_current_symbol
                        ), 'rb')
                        vv_json = json.load(vv_file)
                        if 'error' in vv_json or \
                                ('message' in vv_json and
                                    vv_json['message'] == 'Internal Server Error'):
                            continue
                        else:
                            # log('DEBUG', 'Symbol change for HGNC {0}, MD {1}'.format(hgnc_current_symbol, md_gene['name'][0]))
                            curs.execute(
                                """
                                UPDATE gene
                                SET name[1] = %s, second_name = %s || ',' || %s
                                WHERE hgnc_id = %s
                                """,
                                (
                                    hgnc_current_symbol,
                                    md_gene['name'][0],
                                    md_gene['second_name'],
                                    hgnc_id
                                )
                            )
                            db.commit()
                    if md_gene['hgnc_name'] != gene_info[2]:
                        # log('DEBUG', 'Name differs for HGNC:{0}, MD {1}'.format(gene_info[2], md_gene['name'][0]))
                        curs.execute(
                            """
                            UPDATE gene
                            SET hgnc_name = %s
                            WHERE hgnc_id = %s
                            """,
                            (
                                gene_info[2],
                                hgnc_id
                            )
                        )
                        db.commit()
                else:
                    # gene not yet in MD
                    download_vv_file(hgnc_current_symbol, hgnc_current_symbol)
                    vv_file = open('{0}{1}.json'.format(
                        md_utilities.local_files['variant_validator']['abs_path'],
                        hgnc_current_symbol
                    ), 'rb')
                    vv_json = json.load(vv_file)
                    insert_dict = {}
                    insert_dict['hgnc_id'] = hgnc_id
                    insert_dict['second_name'] = '{0};{1}'.format(gene_info[10], gene_info[8]).replace('"', '')
                    insert_dict['hgnc_name'] = gene_info[2]
                    insert_dict['uniprot_id'] = gene_info[25]
                    if insert_dict['uniprot_id'] == '':
                        insert_dict['uniprot_id'] = 'NULL'
                    if 'error' in vv_json or \
                            ('message' in vv_json and
                                vv_json['message'] == 'Internal Server Error'):
                        continue
                    if vv_json['current_symbol'] == hgnc_current_symbol:
                        if 'transcripts' in vv_json:
                            for vv_transcript in vv_json['transcripts']:
                                if 'reference' in vv_transcript and \
                                        re.search(r'^NM_\d+\.\d{1,2}', vv_transcript['reference']):
                                    insert_dict['chr'] = vv_transcript['annotations']['chromosome']
                                    if isinstance(insert_dict['chr'], list):
                                        insert_dict['chr'] = insert_dict['chr'][0]
                                    if insert_dict['chr'] not in ncbi_chr:
                                        curs.execute(
                                            """
                                            SELECT ncbi_name, genome_version
                                            FROM chromosomes
                                            WHERE name = %s
                                            """,
                                            (insert_dict['chr'],)
                                        )
                                        ncbi_name = curs.fetchall()
                                        for chrom in ncbi_name:
                                            if chrom['genome_version'] == 'hg19':
                                                ncbi_chr_hg19[insert_dict['chr']] = chrom['ncbi_name']
                                            elif chrom['genome_version'] == 'hg38':
                                                ncbi_chr[insert_dict['chr']] = chrom['ncbi_name']
                                        db.commit()
                                    # chek that we have hg19 and hg38
                                    if ncbi_chr[insert_dict['chr']] in vv_transcript['genomic_spans'] and \
                                            ncbi_chr_hg19[insert_dict['chr']] in vv_transcript['genomic_spans']:
                                        insert_dict['strand'] = '+'
                                        if int(vv_transcript['genomic_spans'][ncbi_chr[insert_dict['chr']]]['orientation']) == -1:
                                            insert_dict['strand'] = '-'
                                        insert_dict['number_of_exons'] = vv_transcript['genomic_spans'][ncbi_chr[insert_dict['chr']]]['total_exons']
                                        insert_dict['ng'] = 'NG_000000.0'
                                        for acc in vv_transcript:
                                            ng_match = re.search(r'^(NG_\d+\.\d{1,2})$', acc)
                                            if ng_match:
                                                insert_dict['ng'] = ng_match.group(1)
                                        insert_dict['np'] = vv_transcript['translation']
                                        with open('gene2ensembl_hs', 'r') as f:
                                            for line in f:
                                                line = line.rstrip('\n')
                                                if re.search(re.split(r'\.', vv_transcript['reference'])[0], line):
                                                    to_ensembl = re.split('\t', line)
                                                    insert_dict['enst'] = re.split(r'\.', to_ensembl[4])[0]
                                                    insert_dict['ensp'] = re.split(r'\.', to_ensembl[6])[0]
                                                    if not re.search(r'^ENSP\d+$', insert_dict['ensp']):
                                                        insert_dict['ensp'] = 'NULL'
                                                    if not re.search(r'^ENST\d+$', insert_dict['enst']):
                                                        insert_dict['enst'] = 'NULL'
                                                    break
                                        insert_dict['canonical'] = 'f'
                                        if vv_transcript['annotations']['mane_select'] is True:
                                            insert_dict['canonical'] = 't'
                                        # afterwards need to look for genes with no canonical
                                        s = ", "
                                        t = "', '"
                                        log('INFO', "INSERT INTO gene (name, {0}) VALUES ('{{\"{1}\",\"{2}\"}}', '{3}')".format(
                                                s.join(insert_dict.keys()),
                                                hgnc_current_symbol,
                                                vv_transcript['reference'],
                                                t.join(map(str, insert_dict.values()))
                                            ).replace("'NULL'", "NULL")
                                            )
                                        curs.execute(
                                            """
                                            INSERT INTO gene (name, {0})
                                            VALUES ('{{\"{1}\",\"{2}\"}}', '{3}')
                                            """.format(
                                                s.join(insert_dict.keys()),
                                                hgnc_current_symbol,
                                                vv_transcript['reference'],
                                                t.join(map(str, insert_dict.values()))
                                            ).replace("'NULL'", "NULL")
                                        )
                                        db.commit()
                                    else:
                                        if re.search(r'NM_\d+\.\d{1,2}', vv_transcript['reference']):
                                            log('WARNING', 'hg19 or hg38 mapping issue for {0}-{1}'.format(vv_transcript['reference'], hgnc_current_symbol))
                        else:
                            log('WARNING', 'No transcript in {0}.json file'.format(hgnc_current_symbol))
                            continue
                db.close()

        print('.', end="", flush=True)


if __name__ == '__main__':
    main()
