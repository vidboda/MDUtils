import re
import sys
import urllib3
import certifi
import json
import psycopg2
import psycopg2.extras
import time
from insert_genes import get_db
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config


def log(level, text):
    print()
    localtime = time.asctime( time.localtime(time.time()) )
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


def main():
    # script meant to be croned to update NM acc versions in MD according to VariantValidator
    # to be ran after uta docker update for example
    # uses VV API genes2transcript
    # https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/NM_130786?content-type=application%2Fjson
    vv_url_base = "https://rest.variantvalidator.org"
    # vv_url_base = "http://0.0.0.0:8000/"

    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)

    curs.execute(  # get genes
        "SELECT name, nm_version FROM gene WHERE canonical = 't' ORDER BY name"
    )
    genes = curs.fetchall()
    count = curs.rowcount
    i = 0
    for gene in genes:
        # log('DEBUG', '{}-{}'.format(gene['name'][0], i))
        i += 1
        if i % 500 == 0:
            log('INFO', '{0}/{1} genes checked'.format(i, count))
        # print("MD------{}".format(gene['name'][1]))
        # get VV info for the gene
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        vv_url = "{0}/VariantValidator/tools/gene2transcripts/{1}?content-type=application/json".format(vv_url_base, gene['name'][1])
        try:
            vv_data = json.loads(http.request('GET', vv_url).data.decode('utf-8'))
        except Exception:
           log('WARNING', 'No value for {0}'.format(gene['name'][0]))
           continue
        if 'transcripts' in vv_data:
            # current_nm = gene['nm_version']
            ts_dict = {}
            for transcript in vv_data['transcripts']:
                # print("VV------{}".format(transcript['reference']))
                match_object = re.search(r'^(N[MR]_\d+)\.(\d{1,2})', transcript['reference'])
                if match_object:
                    nm_acc = match_object.group(1)
                    # if nm_acc == gene['name'][1]:
                    nm_version = match_object.group(2)
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
                if res_nm and \
                        int(res_nm[0]) != int(max_vv_nm):
                    # NEED TO TEST IF THE TRANSCIPT WORKS!!!!!
                    vv_url_var = "{0}/VariantValidator/variantvalidator/GRCh38/{1}.{2}:c.1A>T/all?content-type=application/json".format(vv_url_base, nm, max_vv_nm)
                    log('DEBUG', 'Calling VariantValidator API: {}'.format(vv_url_var))
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
                        curs.execute(
                            "UPDATE gene SET nm_version = %s WHERE name[2] = %s",
                            (max_vv_nm, nm)
                        )
                        log('INFO', "NM UPDATE: gene {0} - {1} modified from {2} to {3}".format(gene['name'][0], nm, res_nm[0], max_vv_nm))
                db.commit()

        print('.', end="", flush=True)


if __name__ == '__main__':
    main()
