import re
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import datetime
import json
from precompute_spipv2 import get_db, log
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import md_utilities


def main():
    # script meant to update variants when the canonical form of a gene is changed
    parser = argparse.ArgumentParser(description='script meant to update variants when the canonical form of a gene is changed',
                                     usage='update_vars_when_iso_change.py -k md_api_key -g gene_hgnc')
    parser.add_argument('-k', '--api-key', default='', required=True,
                        help='Your API key visible on your profile page on the website.')
    parser.add_argument('-g', '--gene-symbol', default='', required=True,
                        help='The gene you want to update the variants from.')
    args = parser.parse_args()
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    username = None
    if len(args.api_key) != 43:
        log('ERROR', 'Invalid API key, please check it')
    else:
        api_key = args.api_key
        # user
        curs.execute(
            """
            SELECT username
            FROM mobiuser
            WHERE api_key = %s
            """,
            (api_key,)
        )
        res_user = curs.fetchone()
        if res_user is None:
            log('ERROR', 'Unknown API key')
        username = res_user['username']
        log('INFO', f'User: {username}')
    if match_obj := re.search(r'^([\w-]+)$', args.gene_symbol):
        gene_symbol = match_obj[1]
    else:
        log('ERROR', 'Invalid gene name, please check it')
    # date
    today = datetime.datetime.now()
    creation_date = '{0}-{1}-{2}'.format(
        today.strftime("%Y"), today.strftime("%m"), today.strftime("%d")
    )
    # check if gene exists and get new canonical isoform
    curs.execute(
        """
        SELECT DISTINCT(refseq)
        FROM gene
        WHERE gene_symbol = %s
            AND canonical = 't'
        """,
        (gene_symbol,)
    )
    res = curs.fetchone()
    if res is None:
        log(
            'ERROR',
            f'The gene {gene_symbol} is not present in MobiDetails, please check it',
        )
    nm = res['refseq']
    # nm_full = res['nm']
    # get all variants
    curs.execute(
        """
        SELECT a.chr,
            a.pos,
            a.pos_ref,
            a.pos_alt,
            a.g_name,
            b.c_name,
            b.id
        FROM variant a, variant_feature b
        WHERE a.feature_id = b.id
            AND b.gene_symbol = %s
            AND b.refseq != %s
            AND a.genome_version = 'hg38'
        ORDER BY a.pos
        """,
        (gene_symbol, nm)
    )
    res = curs.fetchall()
    if res is None:
        log('ERROR', 'No variant to update')
    vv_url_base = "https://rest.variantvalidator.org"
    for var in res:
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        # vv_url_base = "http://0.0.0.0:8000/"
        vv_url = "{0}/VariantValidator/variantvalidator/GRCh38/{1}-{2}-{3}-{4}/all?content-type=application/json".format(vv_url_base, var['chr'], var['pos'], var['pos_ref'], var['pos_alt'])
        log('DEBUG', f'Calling VariantValidator API: {vv_url}')
        try:
            vv_data = json.loads(http.request('GET', vv_url).data.decode('utf-8'))
            # log('DEBUG', vv_data)
        except Exception:
                log('WARNING', 'No VV result for {0}:{1}'.format(nm, var['c_name']))
                continue
        for first_level_key in vv_data:
            if match_obj := re.search(f'{nm}:c\.(.+)$', first_level_key):
                new_c_name = match_obj[1]
                log('DEBUG', 'Old c_name: {0} - New c_name: {1}'.format(var['c_name'], new_c_name))
                if new_c_name == var['c_name']:
                    curs.execute(
                        """
                        UPDATE variant_feature
                        SET refseq = %s, creation_date = %s
                        WHERE id = %s
                        """,
                        (nm, creation_date, var['id'])
                    )
                    log('INFO', f"Variant {var['c_name']} remains unchanged")
                else:
                    # likely to change are p_name, ivs_name, prot_type, start_segment_type, start_segment_number, end_segment_type, end_segment_number
                    # also need to update creation_date, creation_user
                    # get p_name
                    p_name = None
                    if 'hgvs_predicted_protein_consequence' in vv_data[first_level_key]:
                        # log('DEBUG', vv_data[first_level_key]['hgvs_predicted_protein_consequence'])
                        if 'tlr' in vv_data[first_level_key]['hgvs_predicted_protein_consequence']:
                            if match_object := re.search(
                                r'NP_\d+\.\d.*:p\.\(?(.+)\)?',
                                vv_data[first_level_key][
                                    'hgvs_predicted_protein_consequence'
                                ]['tlr'],
                            ):
                                p_name = match_object[1]
                                if re.search(r'\)$', p_name):
                                    # remove last ')'
                                    p_name = p_name[:-1]
                            else:
                                log('WARNING', 'No p_name in VV results')
                        else:
                            log('WARNING', 'No tlr in VV results')
                    else:
                        log('WARNING', 'No hgvs_predicted_protein_consequence in VV results')
                    start_segment_type = start_segment_number = end_segment_type = end_segment_number = ivs_name = None
                    # get segments type and number
                    positions = md_utilities.compute_start_end_pos(var['g_name'])
                    ncbi_chr = md_utilities.get_ncbi_chr_name(db, f"chr{var['chr']}", 'hg38')
                    start_segment_type = md_utilities.get_segment_type_from_vv(vv_data[first_level_key]['variant_exonic_positions'][ncbi_chr[0]]['start_exon'])
                    start_segment_number = md_utilities.get_segment_number_from_vv(vv_data[first_level_key]['variant_exonic_positions'][ncbi_chr[0]]['start_exon'])
                    end_segment_type = md_utilities.get_segment_type_from_vv(vv_data[first_level_key]['variant_exonic_positions'][ncbi_chr[0]]['end_exon'])
                    end_segment_number = md_utilities.get_segment_number_from_vv(vv_data[first_level_key]['variant_exonic_positions'][ncbi_chr[0]]['end_exon'])
                        # get IVS name
                    if start_segment_type == 'intron':
                        if positions[0] == positions[1]:
                            ivs_obj = re.search(r'^[\*-]?\d+([\+-]\d+)(.+)$', new_c_name)
                            ivs_name = 'IVS{0}{1}{2}'.format(start_segment_number, ivs_obj[1], ivs_obj[2])
                        elif ivs_obj := re.search(
                            r'^\d+([\+-]\d+)_\d+([\+-]\d+)(.+)$', new_c_name
                        ):
                            ivs_name = 'IVS{0}{1}_IVS{2}{3}{4}'.format(
                                start_segment_number,
                                ivs_obj[1],
                                end_segment_number,
                                ivs_obj[2],
                                ivs_obj[3],
                            )
                        elif ivs_obj := re.search(
                            r'^\d+([\+-]\d+)_(\d+)([^\+-].+)$', new_c_name
                        ):
                            ivs_name = 'IVS{0}{1}_{2}{3}'.format(
                                start_segment_number,
                                ivs_obj[1],
                                ivs_obj[2],
                                ivs_obj[3],
                            )
                        elif ivs_obj := re.search(
                            r'^(\d+)_\d+([\+-]\d+)(.+)$', new_c_name
                        ):
                            ivs_name = '{0}_IVS{1}{2}{3}'.format(
                                ivs_obj[1],
                                end_segment_number,
                                ivs_obj[2],
                                ivs_obj[3],
                            )
                    if p_name is None or \
                            start_segment_type is None or \
                            start_segment_number is None or \
                            end_segment_type is None or \
                            end_segment_number is None:
                        log('WARNING', 'A mandatory new parameter is lacking')
                        continue
                    if ivs_name is None:
                        curs.execute(
                            """
                            UPDATE variant_feature
                            SET refseq = %s,
                                c_name = %s,
                                p_name = %s,
                                start_segment_type = %s,
                                start_segment_number = %s,
                                end_segment_type = %s,
                                end_segment_number = %s,
                                creation_date = %s
                            WHERE id = %s
                            """,
                            (nm, new_c_name, p_name, start_segment_type, start_segment_number, end_segment_type, end_segment_number, creation_date, var['id'])
                        )
                    else:
                        curs.execute(
                            """
                            UPDATE variant_feature
                            SET refseq = %s,
                                c_name = %s,
                                p_name = %s,
                                ivs_name = %s,
                                start_segment_type = %s,
                                start_segment_number = %s,
                                end_segment_type = %s,
                                end_segment_number = %s,
                                creation_date = %s
                            WHERE id = %s
                            """,
                            (nm, new_c_name, p_name, ivs_name, start_segment_type, start_segment_number, end_segment_type, end_segment_number, creation_date, var['id'])
                        )
                    log('INFO', 'Variant {0} updated to {1}'.format(var['c_name'], new_c_name))

    db.commit()
    db_pool.putconn(db)


if __name__ == '__main__':
    main()
