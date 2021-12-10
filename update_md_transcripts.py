import os
import sys
import re
import urllib3
import certifi
import json
import psycopg2
import psycopg2.extras
import time
import hashlib
from insert_genes import get_db
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config, md_utilities


def log(level, text):
    print()
    localtime = time.asctime(time.localtime(time.time()))
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


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
    vv_url_base = "https://www608.lamp.le.ac.uk"
    # vv_url_base = "http://0.0.0.0:8000/"

    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    curs.execute(  # get genes - one transcript per gene (canonical) - allows upadte of all trasncripts
        "SELECT name, second_name, prot_name, prot_short, np, chr, number_of_exons, ng, hgnc_id, strand, uniprot_id FROM gene WHERE canonical = 't' AND (name[1] like 'T%' OR name[1] like 'U%' OR name[1] like 'V%' OR name[1] like 'W%' OR name[1] like 'X%' OR name[1] like 'Y%' OR name[1] like 'Z%') ORDER by name[1]"  # " ORDER BY random()"
    )
    #  WHERE canonical = 't'
    genes = curs.fetchall()
    count = curs.rowcount
    i = 0
    for gene in genes:
        # log('DEBUG', '{}-{}'.format(gene['name'][0], i))
        i += 1
        if i % 500 == 0:
            log('INFO', '{0}/{1} genes checked'.format(i, count))

        curs.execute(
            "SELECT ncbi_name, genome_version FROM chromosomes WHERE name = %s",
            (gene['chr'],)
        )
        ncbi_name = curs.fetchall()
        for chrom in ncbi_name:
            if chrom['genome_version'] == 'hg19':
                hg19_ncbi_chr = chrom['ncbi_name']
            elif chrom['genome_version'] == 'hg38':
                ncbi_chr = chrom['ncbi_name']
        db.commit()
        # log('DEBUG', 'hg19 chr: {0} - hg38 chr: {1}'.format(hg19_ncbi_chr, ncbi_chr))
        # sys.exit()

        # # print("MD------{}".format(gene['name'][1]))
        # ncbi_chr = md_utilities.get_ncbi_chr_name(
        #             db,
        #             'chr{0}'.format(gene['chr']),
        #             'hg38'
        #         )
        # hg19_ncbi_chr = md_utilities.get_ncbi_chr_name(
        #             db,
        #             'chr{0}'.format(gene['chr']),
        #             'hg19'
        #         )
        # get VV info for the gene
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        vv_url = "{0}/VariantValidator/tools/gene2transcripts/{1}?content-type=application/json".format(vv_url_base, gene['name'][1])
        # log('DEBUG', 'Calling VariantValidator gene API: {}'.format(vv_url))
        try:
            vv_data = json.loads(http.request('GET', vv_url).data.decode('utf-8'))
        except Exception:
            log('WARNING', 'No value for {0}'.format(gene['name'][0]))
            continue
        # Store json file in /mobidic_resources/variantValidator/genes/
        if not os.path.isfile(
            '{0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                gene['name'][0]
            )
                ):
            # copy in file system
            with open(
                '{0}{1}.json'.format(
                    md_utilities.local_files['variant_validator']['abs_path'],
                    gene['name'][0]
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
            log('INFO', "VV JSON file copied for gene {0}-{1}".format(gene['name'][0], gene['name'][1]))
        else:
            # md5 to check if there is an update
            new_vv_json_hash = hashlib.md5(json.dumps(vv_data, ensure_ascii=False, indent=4).encode()).hexdigest()
            # current md5
            current_vv_file_hash = hashlib.md5()
            current_vv_file = open('{0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                gene['name'][0]
            ), 'rb')
            current_vv_file_content = current_vv_file.read()
            current_vv_file_hash.update(current_vv_file_content)
            # log('DEBUG', '{0}-{1}'.format(new_vv_json_hash, current_vv_file_hash.hexdigest()))
            if new_vv_json_hash != current_vv_file_hash.hexdigest():
                # copy new file
                # copy in file system
                with open(
                    '{0}{1}.json'.format(
                        md_utilities.local_files['variant_validator']['abs_path'],
                        gene['name'][0]
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
                log('INFO', "VV JSON file replaced for gene {0}-{1}".format(gene['name'][0], gene['name'][1]))

        # check prot name
        if 'current_name' in vv_data and \
                'current_symbol' in vv_data:
            if vv_data['current_symbol'] == gene['name'][0] and \
                    vv_data['current_name'] != gene['prot_name']:
                curs.execute(
                    "UPDATE gene SET prot_name = %s WHERE name[1] = %s",
                    (vv_data['current_name'], gene['name'][0])
                )
                db.commit()
                # log('INFO', "PROTEIN NAME UPDATE: gene {0} modified from {1} to {2}".format(gene['name'][0], gene['prot_name'], vv_data['current_name']))

        if 'transcripts' in vv_data:
            for transcript in vv_data['transcripts']:
                # update np
                if not gene['np']:
                    curs.execute(
                        "UPDATE gene SET np = %s WHERE name[2] = %s",
                        (transcript['translation'], gene['name'][1])
                    )
                # update nb of exons and NP acc no
                if ncbi_chr in transcript['genomic_spans']:
                    nb_exons = transcript['genomic_spans'][ncbi_chr]['total_exons']
                    # log('DEBUG', '{0}-{1}'.format(transcript['reference'], gene['name'][1]))
                    if transcript['reference'] == gene['name'][1]:
                        # log('DEBUG', 'Current #exons: {0} - VV #exons: {1}'.format(gene['number_of_exons'], nb_exons))
                        if int(nb_exons) != int(gene['number_of_exons']):
                            curs.execute(
                                "UPDATE gene SET number_of_exons = %s WHERE name[2] = %s",
                                (nb_exons, gene['name'][1])
                            )
                            log('INFO', "NB EXONS UPDATE: gene {0}-{1} modified from {2} to {3}".format(gene['name'][0], gene['name'][1], gene['number_of_exons'], nb_exons))
                        if transcript['translation'] != gene['np']:
                            curs.execute(
                                "UPDATE gene SET np = %s WHERE name[2] = %s",
                                (transcript['translation'], gene['name'][1])
                            )
                            log('INFO', "NP UPDATE: gene {0} modified from {1} to {2}".format(gene['name'][0], gene['np'], transcript['translation']))
                        db.commit()
                if ncbi_chr in transcript['genomic_spans'] and \
                        hg19_ncbi_chr in transcript['genomic_spans']:
                    if re.search(r'^NM_\d+\.\d{1,2}$', transcript['reference']):
                        # transcript suitable for MD
                        match_object = re.search(r'^(N[MR]_\d+\.\d{1,2})', transcript['reference'])
                        if match_object:
                            nm_acc = match_object.group(1)
                            curs.execute(
                                "SELECT name FROM gene WHERE name[2] = %s",
                                (nm_acc,)
                            )
                            res_nm = curs.fetchone()
                            if res_nm:
                                # trancript already recorded
                                continue
                            # NEED TO TEST IF THE TRANSCRIPT WORKS!!!!!
                        vv_url_var = "{0}/VariantValidator/variantvalidator/GRCh38/{1}:c.2del/all?content-type=application/json".format(vv_url_base, nm_acc)
                        log('DEBUG', 'Calling VariantValidator API: {}'.format(vv_url_var))
                        try:
                            vv_data_var = json.loads(http.request('GET', vv_url_var).data.decode('utf-8'))
                            # log('DEBUG', vv_data)
                        except Exception:
                            log('WARNING', 'No VV result for {0}'.format(nm))
                            continue
                        noupdate = None
                        for first_level_key in vv_data_var:
                            if 'validation_warnings' in vv_data_var[first_level_key]:
                                for warning in vv_data_var[first_level_key]['validation_warnings']:
                                    if re.search(r'cannot be mapped directly to genome build', warning) or \
                                            re.search(r'No transcript definition for', warning) or \
                                            re.search(r'No transcripts found', warning) or \
                                            re.search(r'expected one of', warning):
                                        log('WARNING', "Cannot update gene {0} ({1}) because of {2}".format(gene['name'][0], nm_acc, warning))
                                        noupdate = 1
                                        break
                        if not noupdate:
                            # create new transcript and check if should become canonical
                            # we NEED
                            # name = {MD,VV} second_name = gene['second_name'], chr = gene['chr'], strand = gene['strand']
                            # number_of_exons = VV, prot_name = VV, prot_short = gene['prot_short'], prot_size = VV
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
                            insert_dict['prot_short'] = gene['prot_short'].replace("'", "''")
                            insert_dict['hgnc_id'] = gene['hgnc_id']
                            insert_dict['variant_creation'] = 'ok'
                            insert_dict['number_of_exons'] = transcript['genomic_spans'][ncbi_chr]['total_exons']
                            insert_dict['prot_name'] = vv_data['current_name'].replace("'", "''")
                            insert_dict['prot_size'] = int((transcript['coding_end'] - transcript['coding_start'] + 1) / 3)
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
                            # if gene2ensembl:
                            #     to_ensembl = re.split('\t', gene2ensembl)
                            #     insert_dict['enst'] = re.split(r'\.', to_ensembl[4])[0]
                            #     insert_dict['ensp'] = re.split(r'\.', to_ensembl[6])[0]
                            # log('DEBUG', 'ENST: {0} - ENSP: {1}'.format(insert_dict['enst'], insert_dict['ensp']))
                            # add transcript to MD
                            insert_dict['uniprot_id'] = 'NULL'
                            if re.split(r'\.', gene['name'][1])[0] == re.split(r'\.', transcript['reference'])[0]:
                                # uniprot id
                                if gene['uniprot_id']:
                                    insert_dict['uniprot_id'] = gene['uniprot_id']
                            s = ", "
                            t = "', '"
                            # log('INFO', "INSERT INTO gene (name, {0}) VALUES ('{{\"{1}\",\"{2}\"}}', '{3}')".format(
                            #         s.join(insert_dict.keys()),
                            #         gene['name'][0],
                            #         transcript['reference'],
                            #         t.join(map(str, insert_dict.values()))
                            #     ).replace("'NULL'", "NULL")
                            # )
                            curs.execute(
                                "INSERT INTO gene (name, {0}) VALUES ('{{\"{1}\",\"{2}\"}}', '{3}')".format(
                                    s.join(insert_dict.keys()),
                                    gene['name'][0],
                                    transcript['reference'],
                                    t.join(map(str, insert_dict.values()))
                                ).replace("'NULL'", "NULL")
                            )
                            db.commit()
                            if re.split(r'\.', gene['name'][1])[0] == re.split(r'\.', transcript['reference'])[0]:
                                # if (transcript['annotations']['refseq_select'] or
                                #         transcript['annotations']['mane_select']) and \
                                if re.split(r'\.', gene['name'][1])[1] < re.split(r'\.', transcript['reference'])[1]:
                                    # reset MD canonical for this gene and set it for this transcript
                                    # log('INFO', 'Updating canonical for gene {0} from {1} to {2}'.format(gene['name'][0], gene['name'][1], transcript['reference']))
                                    curs.execute(
                                        "UPDATE gene SET canonical = 'f' WHERE name[1] = %s",
                                        (gene['name'][0],)
                                    )
                                    curs.execute(
                                        "UPDATE gene SET canonical = 't' WHERE name[2] = %s",
                                        (transcript['reference'],)
                                    )
                                    db.commit()
                else:
                    log('WARNING', 'Transcript {0} from gene {1} has hg19/38 mapping issues'.format(transcript['reference'], gene['name'][0]))
            # db.commit()
        # db.commit()
        print('.', end="", flush=True)


if __name__ == '__main__':
    main()
