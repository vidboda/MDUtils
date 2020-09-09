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
        "SELECT name, nm_version FROM gene WHERE canonical = 't' AND name[1] ~'^R[HIM]\w+'ORDER BY name"
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
                    # if VV nm_version > current_nm
                    # => update and next
                    # if VV nm_version < current_nm
                    # 2 possibilities
                    # - update needed (MD too high)
                    # - update not needed => current_nm will arrive next
                    # so we need to keep in memory the preceeding nm_version and decide whether or not we have to update
                    # if nm_acc == gene['name'][1]:
                    #     if int(nm_version) > int(current_nm):
                    #         # print("IN--{}-{}-{}-{}-{}-".format(vv_data['current_symbol'],
                    #         #   transcript['reference'], gene['name'][1], nm_version, current_nm))
                    #         curs.execute(
                    #             "UPDATE gene SET nm_version = '{0}' WHERE name[2] = '{1}'".format(nm_version, gene['name'][1])
                    #         )
                    #         db.commit()
                    #         current_nm = nm_version
                    #         log('INFO', "NM INCREASE: gene {0} - {1} increases from {2} to {3}".format(gene['name'][0], gene['name'][1], current_nm, nm_version))
                        #  elif int(nm_version) < int(current_nm):
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
                    curs.execute(
                        "UPDATE gene SET nm_version = %s WHERE name[2] = %s",
                        (max_vv_nm, nm)
                    )
                    log('INFO', "NM UPDATE: gene {0} - {1} modified from {2} to {3}".format(gene['name'][0], nm, res_nm[0], max_vv_nm))
                db.commit()

        print('.', end="", flush=True)
if __name__ == '__main__':
    main()
