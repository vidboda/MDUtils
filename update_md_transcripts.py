import os
import re
# import time
import urllib3
import certifi
import json
import psycopg2
import psycopg2.extras
import hashlib
import argparse
from precompute_spipv2 import get_db, log
from MobiDetailsApp import md_utilities


def call_vv(vv_url, http, curs, db):
    try:
        vv_data = json.loads(http.request('GET', vv_url).data.decode('utf-8'))
    except Exception:
        log('WARNING', 'No value for {0}'.format(gene['gene_symbol']))
        # disable in MD
        curs.execute(
            """
            UPDATE gene
            SET variant_creation = 'not_in_vv_json'
            WHERE gene_symbol = %s
            """,
            (gene['gene_symbol'],)
        )
        db.commit()
        vv_data['mdutils'] = 'exception'
    return vv_data

def main():
    # script meant to be croned to update transcripts in MD according to VariantValidator
    # to be ran after uta docker update for example
    # or once per 2-3 months
    # after successfull run on dev server:
    # - run check_uniprot_ids.py script
    # - rsync resources/variantvalidator/genes/*.json folder to the prod
    # - run the update_genes_from_remote.py script from the prod server
    # uses VV API genes2transcript
    # https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/NM_130786?content-type=application%2Fjson
    # vv_url_base = "https://rest.variantvalidator.org"
    # vv_url_base = "https://www608.lamp.le.ac.uk"
    vv_url_base = "http://rvv.chu-montpellier.fr"

    parser = argparse.ArgumentParser(description='Update MD gene files from VV', usage='python update_md_transcripts.py [-f path/to/dir/containing/genes/file.txt]')
    parser.add_argument('-f', '--file', default='', required=False, help='Path to the genes file to be updated')
    args = parser.parse_args()
    # get sql file list
    gene_file = None
    gene_list = ''
    if args.file and \
            os.path.isfile(args.file):
        gene_file = args.file
        for line in open(gene_file).readlines():
            gene_list += "'{0}', ".format(line.rstrip())
        gene_list = gene_list[:-2]
    elif args.file:
        log('ERROR', 'Invalid input path for gene file, please check your command')
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    if not gene_file:
        curs.execute(  # get genes - one transcript per gene (canonical) - allows update of all transcripts
            """
            SELECT gene_symbol, refseq, second_name, hgnc_name, np, chr, number_of_exons, ng, hgnc_id, strand, uniprot_id
            FROM gene
            WHERE canonical = 't'
            ORDER by gene_symbol
            """  # " ORDER BY random()"
        )
        #  WHERE canonical = 't'
    else:
        log('INFO', 'The following genes will be considered {0}'.format(gene_list))
        curs.execute(  # get genes - one transcript per gene (canonical) - allows update of all transcripts from the list
            """
            SELECT gene_symbol, refseq, second_name, hgnc_name, np, chr, number_of_exons, ng, hgnc_id, strand, uniprot_id
            FROM gene
            WHERE canonical = 't'
                AND gene_symbol IN ({0})
            ORDER by gene_symbol
            """.format(gene_list)  # " ORDER BY random()"
        )
    genes = curs.fetchall()
    count = curs.rowcount
    i = 0
    for gene in genes:
        # log('DEBUG', '{}-{}'.format(gene['gene_symbol'], i))
        i += 1
        if i % 500 == 0:
            log('INFO', '{0}/{1} genes checked'.format(i, count))

        curs.execute(
            """
            SELECT ncbi_name, genome_version
            FROM chromosomes
            WHERE name = %s
            """,
            (gene['chr'],)
        )
        ncbi_name = curs.fetchall()
        for chrom in ncbi_name:
            if chrom['genome_version'] == 'hg19':
                hg19_ncbi_chr = chrom['ncbi_name']
            elif chrom['genome_version'] == 'hg38':
                ncbi_chr = chrom['ncbi_name']
        # db.commit()
        # get VV info for the gene
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        vv_url = "{0}/VariantValidator/tools/gene2transcripts/{1}?content-type=application/json".format(vv_url_base, gene['gene_symbol'])
        # log('DEBUG', 'Calling VariantValidator gene API: {}'.format(vv_url))
        vv_data = call_vv(vv_url, http, curs, db)
        if 'mdutils' in vv_data:
            continue
        # Store json file in /mobidic_resources/variantValidator/genes/
        if 'error' in vv_data:
            log('WARNING', 'VV error for gene {0} - changed gene symbol?'.format(gene['gene_symbol']))
            # try querying by refseq instead of HGNC symbol
            vv_url = "{0}/VariantValidator/tools/gene2transcripts/{1}?content-type=application/json".format(vv_url_base, gene['refseq'])
            # log('DEBUG', 'Calling VariantValidator gene API: {}'.format(vv_url))
            vv_data = call_vv(vv_url, http, curs, db)
            if 'mdutils' in vv_data:
                continue
            if 'error' in vv_data:
                curs.execute(
                    """
                    UPDATE gene
                    SET variant_creation = 'not_in_vv_json'
                    WHERE gene_symbol = %s
                    """,
                    (gene['gene_symbol'],)
                )
                db.commit()
                continue
        if not 'current_symbol' in vv_data:
            log('DEBUG', vv_data)
        # check gene symbol
        if gene['gene_symbol'] != vv_data['current_symbol']:
            log('WARNING', 'Gene symbol to change: {0} became {1}'.format(gene['gene_symbol'], vv_data['current_symbol']))
            # update gene symbol
            curs.execute(
                """
                UPDATE gene
                SET gene_symbol = %s,
                second_name = CONCAT(second_name, ",", %s),
                hgnc_name = %s
                WHERE gene_symbol = %s
                """,
                (
                    vv_data['current_symbol'],
                    gene['gene_symbol'],
                    vv_data['current_name'],
                    gene['gene_symbol']
                )
            )
            db.commit()
            # reload vv_data with gene_symbol
            vv_url = "{0}/VariantValidator/tools/gene2transcripts/{1}?content-type=application/json".format(vv_url_base, vv_data['current_symbol'])
            vv_data = call_vv(vv_url, http, curs, db)
        # log('DEBUG', 'VV genes files path: {0}'.format(md_utilities.local_files['variant_validator']['abs_path']))
        if not os.path.isfile(
            '{0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                gene['gene_symbol']
            )
                ):
            # copy in file system
            with open(
                '{0}{1}.json'.format(
                    md_utilities.local_files['variant_validator']['abs_path'],
                    gene['gene_symbol']
                ),
                "w",
                encoding='utf-8'
            ) as vv_file:
                json.dump(
                    vv_data,
                    vv_file,
                    ensure_ascii=False,
                    indent=4
                )
            log('INFO', "VV JSON file copied for gene {0}-{1}".format(gene['gene_symbol'], gene['refseq']))
            continue
        else:
            # # md5 to check if there is an update
            new_vv_json_hash = hashlib.md5(json.dumps(vv_data, ensure_ascii=False, indent=4).encode()).hexdigest()
            # current md5
            current_vv_file_hash = hashlib.md5()
            current_vv_file = open('{0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                gene['gene_symbol']
            ), 'rb')
            current_vv_file_content = current_vv_file.read()
            current_vv_file_hash.update(current_vv_file_content)
            # log('DEBUG', '{0}-{1}'.format(new_vv_json_hash, current_vv_file_hash.hexdigest()))
            if new_vv_json_hash != current_vv_file_hash.hexdigest():
                # copy new file
                # copy in file system
                # but check before that there is sthg in the file (not Internal Server Error to avoid replacing a good file with a bad)
                if 'transcripts' in vv_data:
                    with open(
                        '{0}{1}.json'.format(
                            md_utilities.local_files['variant_validator']['abs_path'],
                            gene['gene_symbol']
                        ),
                        "w",
                        encoding='utf-8'
                    ) as vv_file:
                        json.dump(
                            vv_data,
                            vv_file,
                            ensure_ascii=False,
                            indent=4
                        )
                    log('INFO', "VV JSON file replaced for gene {0}-{1}".format(gene['gene_symbol'], gene['refseq']))
            current_vv_file.close()
        # check hgnc name
        if 'current_name' in vv_data and \
                'current_symbol' in vv_data:
            if vv_data['current_symbol'] == gene['gene_symbol'] and \
                    vv_data['current_name'] != gene['hgnc_name']:
                curs.execute(
                    """
                    UPDATE gene
                    SET hgnc_name = %s
                    WHERE gene_symbol = %s
                    """,
                    (vv_data['current_name'], gene['gene_symbol'])
                )
                db.commit()
                # log('INFO', "PROTEIN NAME UPDATE: gene {0} modified from {1} to {2}".format(gene['name'][0], gene['hgnc_name'], vv_data['current_name']))

        if 'transcripts' in vv_data:
            for transcript in vv_data['transcripts']:
                # update np
                if not gene['np']:
                    curs.execute(
                        """
                        UPDATE gene
                        SET np = %s
                        WHERE refseq = %s
                        """,
                        (transcript['translation'], gene['refseq'])
                    )
                # update nb of exons and NP acc no
                if ncbi_chr in transcript['genomic_spans']:
                    #  and \
                    #     hg19_ncbi_chr in transcript['genomic_spans']:
                    nb_exons = transcript['genomic_spans'][ncbi_chr]['total_exons']
                    # log('DEBUG', '{0}-{1}'.format(transcript['reference'], gene['refseq']))
                    if transcript['reference'] == gene['refseq']:
                        # log('DEBUG', 'Current #exons: {0} - VV #exons: {1}'.format(gene['number_of_exons'], nb_exons))
                        if int(nb_exons) != int(gene['number_of_exons']):
                            curs.execute(
                                """
                                UPDATE gene
                                SET number_of_exons = %s
                                WHERE refseq = %s
                                """,
                                (nb_exons, gene['refseq'])
                            )
                            log('INFO', "NB EXONS UPDATE: gene {0}-{1} modified from {2} to {3}".format(gene['gene_symbol'], gene['refseq'], gene['number_of_exons'], nb_exons))
                        if transcript['translation'] != gene['np']:
                            curs.execute(
                                """
                                UPDATE gene
                                SET np = %s
                                WHERE refseq = %s
                                """,
                                (transcript['translation'], gene['refseq'])
                            )
                            log('INFO', "NP UPDATE: gene {0} modified from {1} to {2}".format(gene['gene_symbol'], gene['np'], transcript['translation']))
                        db.commit()
                    # if ncbi_chr in transcript['genomic_spans'] and \
                    #         hg19_ncbi_chr in transcript['genomic_spans']:
                        # continue
                    if re.search(r'^NM_\d+\.\d{1,2}$', transcript['reference']):
                        not_new = 0
                        # transcript suitable for MD
                        match_object = re.search(r'^(N[MR]_\d+\.\d{1,2})', transcript['reference'])
                        if match_object:
                            nm_acc = match_object.group(1)
                            curs.execute(
                                """
                                SELECT gene_symbol, variant_creation
                                FROM gene
                                WHERE refseq = %s
                                """,
                                (nm_acc,)
                            )
                            res_nm = curs.fetchone()
                            if res_nm:
                                # trancript already recorded - check variant creation
                                if res_nm['variant_creation'] == 'ok':
                                    continue
                                else:
                                    not_new = 1
                            # NEED TO TEST IF THE TRANSCRIPT WORKS!!!!!
                        vv_url_var = "{0}/VariantValidator/variantvalidator/GRCh38/{1}:c.2del/all?content-type=application/json".format(vv_url_base, nm_acc)
                        # log('DEBUG', 'Calling VariantValidator API: {}'.format(vv_url_var))
                        try:
                            vv_data_var = json.loads(http.request('GET', vv_url_var).data.decode('utf-8'))
                            # log('DEBUG', vv_data)
                        except Exception:
                            log('WARNING', 'No VV result for {0}'.format(nm_acc))
                            continue
                        noupdate = None
                        if 'message' in vv_data_var and \
                                vv_data_var['message'] == 'Internal Server Error':
                            curs.execute(
                                """
                                UPDATE gene
                                SET variant_creation = 'vv_server_error'
                                WHERE refseq = %s
                                """,
                                (nm_acc,)
                            )
                            db.commit()
                            # no need to update again
                            not_new = 0
                        for first_level_key in vv_data_var:
                            if 'validation_warnings' in vv_data_var[first_level_key]:
                                for warning in vv_data_var[first_level_key]['validation_warnings']:
                                    if re.search(r'cannot be mapped directly to genome build', warning) or \
                                            re.search(r'No transcript definition for', warning) or \
                                            re.search(r'No transcripts found', warning) or \
                                            re.search(r'expected one of', warning):
                                        log('WARNING', "Cannot update gene {0} ({1}) because of {2}".format(gene['gene_symbol'], nm_acc, warning))
                                        noupdate = 1
                                        break
                        variant_creation = 'ok'
                        if hg19_ncbi_chr not in transcript['genomic_spans']:
                            # look for hg19_mapping_default which is now allowed but tagged
                            variant_creation = 'hg19_mapping_default'
                        # log('DEBUG', 'transcript: {0} - variant_creation: {1}').format()
                        if not noupdate and \
                                not_new == 1:
                            # if vv creation ok, update
                            # case of an existing gene which was identified as not_ok but now works
                            curs.execute(
                                """
                                UPDATE gene
                                SET variant_creation = %s
                                WHERE refseq = %s
                                """,
                                (variant_creation, nm_acc,)
                            )
                            db.commit()
                        elif not noupdate and \
                                not_new == 0:
                            # create new transcript and check if should become canonical
                            # we NEED
                            # name = {MD,VV} second_name = gene['second_name'], chr = gene['chr'], strand = gene['strand']
                            # number_of_exons = VV, hgnc_name = VV, prot_short = gene['prot_short'], prot_size = VV
                            # uniprot_id = NULL, ng = VV ou gene['ng'], enst = gene2ensembl, ensp = gene2ensembl
                            # canonical = MD update from current, variant_creation = ok, hgnc_id = gene['hgnc_id']
                            insert_dict = {}
                            insert_dict['second_name'] = 'NULL'
                            # if gene['name'][0] == 'ARHGEF33':
                            #     log('DEBUG', 'second_name: {0}'.format(gene['second_name']))
                            if gene['second_name']:
                                insert_dict['second_name'] = gene['second_name'].replace("'", "''")
                            insert_dict['chr'] = gene['chr']
                            insert_dict['strand'] = '+'
                            if int(transcript['genomic_spans'][ncbi_chr]['orientation']) == -1:
                                insert_dict['strand'] = '-'
                            insert_dict['np'] = transcript['translation']
                            insert_dict['hgnc_id'] = gene['hgnc_id']
                            insert_dict['variant_creation'] = variant_creation
                            insert_dict['number_of_exons'] = transcript['genomic_spans'][ncbi_chr]['total_exons']
                            insert_dict['hgnc_name'] = vv_data['current_name'].replace("'", "''")
                            insert_dict['prot_size'] = int((transcript['coding_end'] - transcript['coding_start'] + 1) / 3)-1
                            for acc in transcript:
                                ng_match = re.search(r'^(NG_\d+\.\d{1,2})$', acc)
                                if ng_match:
                                    insert_dict['ng'] = ng_match.group(1)
                                else:
                                    insert_dict['ng'] = gene['ng']
                            # gene2ensembl = subprocess.call(['/usr/bin/grep', re.split('.', transcript['reference'])[0], 'gene2ensembl'])
                            with open('gene2ensembl_hs', 'r') as f:
                                for line in f:
                                    line = line.rstrip('\n')
                                    if re.search(re.split(r'\.', transcript['reference'])[0], line):
                                        # log('DEBUG', 'NM to look for: {0}'.format(re.split(r'\.', transcript['reference'])[0]))
                                        to_ensembl = re.split('\t', line)
                                        insert_dict['enst'] = re.split(r'\.', to_ensembl[4])[0]
                                        insert_dict['ensp'] = re.split(r'\.', to_ensembl[6])[0]
                                        if not re.search(r'^ENSP\d+$', insert_dict['ensp']):
                                            insert_dict['ensp'] = 'NULL'
                                        if not re.search(r'^ENST\d+$', insert_dict['enst']):
                                            insert_dict['enst'] = 'NULL'
                                        break
                            # add transcript to MD
                            insert_dict['uniprot_id'] = 'NULL'
                            if re.split(r'\.', gene['refseq'])[0] == re.split(r'\.', transcript['reference'])[0]:
                                # uniprot id
                                if gene['uniprot_id']:
                                    insert_dict['uniprot_id'] = gene['uniprot_id']
                            s = ", "
                            t = "', '"
                            insert_dict['gene_symbol'] = gene['gene_symbol']
                            insert_dict['refseq'] = transcript['reference']
                            curs.execute(
                                """
                                INSERT INTO gene ({0})
                                VALUES ('{1}')
                                """.format(
                                    s.join(insert_dict.keys()),
                                    t.join(map(str, insert_dict.values()))
                                ).replace("'NULL'", "NULL")
                            )
                            db.commit()
                            log('INFO', 'Inserting transcript {0} for gene {1}'.format(insert_dict['refseq'], insert_dict['gene_symbol']))
                            # check MANE to possibly update the canonical field - or increased NM version 
                            if transcript['annotations']['mane_select'] is True or \
                                    (re.split(r'\.', gene['refseq'])[0] == re.split(r'\.', transcript['reference'])[0] and
                                    re.split(r'\.', gene['refseq'])[1] < re.split(r'\.', transcript['reference'])[1]):
                                    # reset MD canonical for this gene and set it for this transcript
                                    # log('INFO', 'Updating canonical for gene {0} from {1} to {2}'.format(gene['name'][0], gene['name'][1], transcript['reference']))
                                    curs.execute(
                                        """
                                        UPDATE gene
                                        SET canonical = 'f'
                                        WHERE gene_symbol = %s
                                        """,
                                        (gene['gene_symbol'],)
                                    )
                                    curs.execute(
                                        """
                                        UPDATE gene
                                        SET canonical = 't'
                                        WHERE refseq = %s
                                        """,
                                        (transcript['reference'],)
                                    )
                                    db.commit()
                else:
                    if re.search(r'NM_\d+\.\d{1,2}', transcript['reference']):
                        # log('WARNING', 'Transcript {0} from gene {1} has hg19/38 mapping issues'.format(transcript['reference'], gene['gene_symbol']))
                        default = 'hg19_mapping_default'
                        if ncbi_chr not in transcript['genomic_spans'] and \
                                hg19_ncbi_chr not in transcript['genomic_spans']:
                            default = 'mapping_default'
                        elif ncbi_chr not in transcript['genomic_spans'] and \
                                hg19_ncbi_chr in transcript['genomic_spans']:
                            default = 'hg38_mapping_default'
                        curs.execute(
                            """
                            UPDATE gene
                            SET variant_creation = %s
                            WHERE refseq = %s
                            """,
                            (default, transcript['reference'])
                        )
                        db.commit()
        # db.commit()
        print('.', end="", flush=True)
    db_pool.putconn(db)


if __name__ == '__main__':
    main()
