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
    # script meant to be croned to update NM acc versions in MD according to VariantValidator
    # to be ran after uta docker update for example
    # uses VV API genes2transcript
    # https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/NM_130786?content-type=application%2Fjson
    #vv_url_base = "https://rest.variantvalidator.org"
    vv_url_base = "https://www608.lamp.le.ac.uk"
    # vv_url_base = "http://0.0.0.0:8000/"

    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)

    curs.execute(  # get genes
        "SELECT name, prot_name, nm_version, np, chr, number_of_exons FROM gene ORDER BY random() LIMIT 10"
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
        # print("MD------{}".format(gene['name'][1]))
        ncbi_chr = md_utilities.get_ncbi_chr_name(
                    db,
                    'chr{0}'.format(gene['chr']),
                    'hg38'
                )
        # get VV info for the gene
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        vv_url = "{0}/VariantValidator/tools/gene2transcripts/{1}?content-type=application/json".format(vv_url_base, gene['name'][1])
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
            #md5 to check if there is an update
            new_vv_json_hash = hashlib.md5(json.dumps(vv_data, ensure_ascii=False, indent=4).encode()).hexdigest()
            # current md5
            current_vv_file_hash = hashlib.md5()
            current_vv_file = open('{0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                gene['name'][0]
            ), 'rb')
            current_vv_file_content = current_vv_file.read()
            current_vv_file_hash.update(current_vv_file_content)
            log('DEBUG', '{0}-{1}'.format(new_vv_json_hash, current_vv_file_hash.hexdigest()))
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

            # TODO if file exists check md5 and if it's different, update the file - and all the associated variants?????
        # check gene name
        if (
            'current_name' in vv_data
            and 'current_symbol' in vv_data
            and vv_data['current_symbol'] == gene['name'][0]
            and vv_data['current_name'] != gene['prot_name']
        ):
            curs.execute(
                "UPDATE gene SET prot_name = %s WHERE name[1] = %s",
                (vv_data['current_name'], gene['name'][0])
            )
            db.commit()
            log('INFO', "NAME UPDATE: gene {0} modified from {1} to {2}".format(gene['name'][0], gene['prot_name'], vv_data['current_name']))

        if 'transcripts' in vv_data:
            # current_nm = gene['nm_version']
            ts_dict = {}
            for transcript in vv_data['transcripts']:
                # update nb of exons - to be removed later!!!!!
                if ncbi_chr[0] in transcript['genomic_spans']:
                    nb_exons = transcript['genomic_spans'][ncbi_chr[0]]['total_exons']
                    # log('DEBUG', '{0}-{1}'.format(transcript['reference'], '{0}.{1}'.format(gene['name'][1], gene['nm_version'])))
                    if transcript['reference'] == '{0}.{1}'.format(gene['name'][1], gene['nm_version']):
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
                if match_object := re.search(
                    r'^(N[MR]_\d+)\.(\d{1,2})', transcript['reference']
                ):
                    nm_acc = match_object[1]
                    # if nm_acc == gene['name'][1]:
                    nm_version = match_object[2]
                    if nm_acc not in ts_dict:
                        ts_dict[nm_acc] = [nm_version]
                    else:
                        ts_dict[nm_acc].append(nm_version)

            # do sthg with ts_dict before changing gene
            for nm in ts_dict:
                # exploring unconsistant NMs
                curs.execute(
                    "SELECT nm_version FROM gene WHERE name[2] = %s",
                    (nm,)
                )
                res_nm = curs.fetchone()
                max_vv_nm = max(ts_dict[nm])
                if not res_nm:
                    continue
                # log("DEBUG", "Gene: {0} - NM: {1} - VV Max NM: {2} - MD Current NM: {3}".format(gene['name'][0], nm, max_vv_nm, res_nm[0]))
                if int(res_nm[0]) != int(max_vv_nm):
                    # NEED TO TEST IF THE TRANSCRIPT WORKS!!!!!
                    vv_url_var = "{0}/VariantValidator/variantvalidator/GRCh38/{1}.{2}:c.1A>T/all?content-type=application/json".format(vv_url_base, nm, max_vv_nm)
                    log('DEBUG', f'Calling VariantValidator API: {vv_url_var}')
                    try:
                        vv_data = json.loads(http.request('GET', vv_url_var).data.decode('utf-8'))
                        # log('DEBUG', vv_data)
                    except Exception:
                            log('WARNING', 'No VV result for {0}.{1}'.format(nm, max_vv_nm))
                            continue
                    noupdate = None
                    for first_level_key in vv_data:
                        if 'validation_warnings' in vv_data[first_level_key]:
                            for warning in vv_data[first_level_key]['validation_warnings']:
                                if re.search(r'cannot be mapped directly to genome build', warning) or \
                                        re.search(r'No transcript definition for', warning) or \
                                        re.search(r'No transcripts found', warning) or \
                                        re.search(r'expected one of', warning):
                                    log('WARNING', "Cannot update gene {0} from {1} to {2} because of {3}".format(gene['name'][0], res_nm[0], max_vv_nm, warning))
                                    noupdate = 1
                                    break
                    if not noupdate:
                        # update nb of exons
                        new_exon_nb = gene['number_of_exons']
                        for transcript in vv_data['transcripts']:
                            if ncbi_chr[0] in transcript['genomic_spans']:
                                nb_exons = transcript['genomic_spans'][ncbi_chr[0]]['total_exons']
                                if int(nb_exons) != int(gene['number_of_exons']):
                                    new_exon_nb = nb_exons
                        curs.execute(
                            "UPDATE gene SET nm_version = %s, number_of_exons = %s WHERE name[2] = %s",
                            (max_vv_nm, new_exon_nb, nm)
                        )
                        log('INFO', "NM UPDATE: gene {0} - {1} modified from {2} to {3} - exons number from {4} to {5}".format(gene['name'][0], nm, res_nm[0], max_vv_nm, gene['number_of_exons'], new_exon_nb))
                db.commit()

        print('.', end="", flush=True)


if __name__ == '__main__':
    main()
