import re
import sys
import argparse
import psycopg2
import psycopg2.extras
import urllib3
import certifi
import time
import datetime
import json
from insert_genes import get_db
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config, md_utilities


def log(level, text):
    localtime = time.asctime(time.localtime(time.time()))
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


# def compute_start_end_pos(name):
#     match_object = re.search(r'(\d+)_(\d+)[di]', name)
#     if match_object is not None:
#         return match_object.group(1), match_object.group(2)
#     else:
#         match_object = re.search(r'^(\d+)[ATGC][>=]', name)
#         if match_object is not None:
#             return match_object.group(1), match_object.group(1)
#         else:
#             # single nt del or delins
#             match_object = re.search(r'^(\d+)[d]', name)
#             if match_object is not None:
#                 return match_object.group(1), match_object.group(1)


def main():
    # script meant to update variants when the canonical form of a gene is changed
    parser = argparse.ArgumentParser(description='Define a canonical transcript per gene when several are defined',
                                     usage='python update_canonical_when_several.py -k md_api_key -g gene_hgnc')
    parser.add_argument('-k', '--api-key', default='', required=True,
                        help='Your API key visible on your profile page on the website.')
    parser.add_argument('-g', '--gene-name', default='', required=True,
                        help='The gene you want to update the variants from.')
    args = parser.parse_args()
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    username = None
    if len(args.api_key) != 43:
        log('ERROR', 'Invalid API key, please check it')
    else:
        api_key = args.api_key
        # user
        curs.execute(
            "SELECT username FROM mobiuser WHERE api_key = %s",
            (api_key,)
        )
        res_user = curs.fetchone()
        if res_user is None:
            log('ERROR', 'Unknown API key')
        username = res_user['username']
        log('INFO', 'User: {}'.format(username))
    match_obj = re.search('^([\w-]+)$', args.gene_name)
    if match_obj:
        gene_name = match_obj.group(1)
    else:
        log('ERROR', 'Invalid gene name, please check it')
    # date
    today = datetime.datetime.now()
    creation_date = '{0}-{1}-{2}'.format(
        today.strftime("%Y"), today.strftime("%m"), today.strftime("%d")
    )
    # check if gene exists and get new canonical isoform
    curs.execute(
        "SELECT DISTINCT(name[2]) as nm FROM gene WHERE name[1] = %s AND canonical = 't'",
        (gene_name,)
    )
    res = curs.fetchone()
    if res is None:
        log('ERROR', 'The gene {} is not present in MobiDetails, please check it'.format(gene_name))
    nm = res['nm']
    # nm_full = res['nm']
    # get all variants
    curs.execute(
        "SELECT a.chr, a.pos, a.pos_ref, a.pos_alt, a.g_name, b.c_name, b.id FROM variant a, variant_feature b WHERE a.feature_id = b.id AND b.gene_name[1] = %s AND b.gene_name[2] != %s AND a.genome_version = 'hg38' ORDER BY a.pos",
        (gene_name, nm)
    )
    res = curs.fetchall()
    if res is None:
        log('ERROR', 'No variant to update')
    for var in res:
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        vv_url_base = "https://rest.variantvalidator.org"
        # vv_url_base = "http://0.0.0.0:8000/"
        vv_url = "{0}/VariantValidator/variantvalidator/GRCh38/{1}-{2}-{3}-{4}/all?content-type=application/json".format(vv_url_base, var['chr'], var['pos'], var['pos_ref'], var['pos_alt'])
        log('DEBUG', 'Calling VariantValidator API: {}'.format(vv_url))
        try:
            vv_data = json.loads(http.request('GET', vv_url).data.decode('utf-8'))
            # log('DEBUG', vv_data)
        except Exception:
                log('WARNING', 'No VV result for {0}:{1}'.format(nm, var['c_name']))
                continue
        for first_level_key in vv_data:
            match_obj = re.search(r'{}:c\.(.+)$'.format(nm), first_level_key)
            if match_obj:
                new_c_name = match_obj.group(1)
                log('DEBUG', 'Old c_name: {0} - New c_name: {1}'.format(var['c_name'], new_c_name))
                if new_c_name == var['c_name']:
                    curs.execute(
                        "UPDATE variant_feature SET gene_name[2] = %s, \
                            creation_date = %s WHERE id = %s",
                        (nm, creation_date, var['id'])
                    )
                    log('INFO', 'Variant {} remains unchanged'.format(var['c_name']))
                else:
                    # likely to change are p_name, ivs_name, prot_type, start_segment_type, start_segment_number, end_segment_type, end_segment_number
                    # also need to update creation_date, creation_user
                    # get p_name
                    p_name = None
                    if 'hgvs_predicted_protein_consequence' in vv_data[first_level_key]:
                        # log('DEBUG', vv_data[first_level_key]['hgvs_predicted_protein_consequence'])
                        if 'tlr' in vv_data[first_level_key]['hgvs_predicted_protein_consequence']:
                            # log('DEBUG', vv_data[first_level_key]['hgvs_predicted_protein_consequence']['tlr'])
                            match_object = re.search(r'NP_\d+\.\d.*:p\.\(?(.+)\)?', vv_data[first_level_key]['hgvs_predicted_protein_consequence']['tlr'])
                            if match_object:
                                p_name = match_object.group(1)
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
                    ncbi_chr = md_utilities.get_ncbi_chr_name(db, 'chr{}'.format(var['chr']), 'hg38')
                    start_segment_type = md_utilities.get_segment_type_from_vv(vv_data[first_level_key]['variant_exonic_positions'][ncbi_chr[0]]['start_exon'])
                    start_segment_number = md_utilities.get_segment_number_from_vv(vv_data[first_level_key]['variant_exonic_positions'][ncbi_chr[0]]['start_exon'])
                    end_segment_type = md_utilities.get_segment_type_from_vv(vv_data[first_level_key]['variant_exonic_positions'][ncbi_chr[0]]['end_exon'])
                    end_segment_number = md_utilities.get_segment_number_from_vv(vv_data[first_level_key]['variant_exonic_positions'][ncbi_chr[0]]['end_exon'])
                    if positions[0] != positions[1]:
                        # get IVS name
                        if start_segment_type == 'intron':
                            ivs_obj = re.search(r'^\d+([\+-]\d+)_\d+([\+-]\d+)(.+)$', new_c_name)
                            if ivs_obj:
                                ivs_name = 'IVS{0}{1}_IVS{2}{3}{4}'.format(
                                    start_segment_number, ivs_obj.group(1),
                                    end_segment_number, ivs_obj.group(2), ivs_obj.group(3)
                                )
                            else:
                                ivs_obj = re.search(r'^\d+([\+-]\d+)_(\d+)([^\+-].+)$', new_c_name)
                                if ivs_obj:
                                    ivs_name = 'IVS{0}{1}_{2}{3}'.format(
                                        start_segment_number, ivs_obj.group(1),
                                        ivs_obj.group(2), ivs_obj.group(3)
                                    )
                                else:
                                    ivs_obj = re.search(r'^(\d+)_\d+([\+-]\d+)(.+)$', new_c_name)
                                    if ivs_obj:
                                        ivs_name = '{0}_IVS{1}{2}{3}'.format(
                                            ivs_obj.group(1), end_segment_number,
                                            ivs_obj.group(2), ivs_obj.group(3)
                                        )
                    else:
                        # substitutions
                        if start_segment_type == 'intron':
                            ivs_obj = re.search(r'^[\*-]?\d+([\+-]\d+)(.+)$', new_c_name)
                            ivs_name = 'IVS{0}{1}{2}'.format(
                                start_segment_number, ivs_obj.group(1), ivs_obj.group(2)
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
                            "UPDATE variant_feature SET gene_name[2] = %s, c_name = %s, p_name = %s, start_segment_type = %s, \
                            start_segment_number = %s, end_segment_type = %s, end_segment_number = %s, \
                            creation_date = %s WHERE id = %s",
                            (nm, new_c_name, p_name, start_segment_type, start_segment_number, end_segment_type, end_segment_number, creation_date, var['id'])
                        )
                    else:
                        curs.execute(
                            "UPDATE variant_feature SET gene_name[2] = %s, c_name = %s, p_name = %s, ivs_name = %s, start_segment_type = %s, \
                            start_segment_number = %s, end_segment_type = %s, end_segment_number = %s, \
                            creation_date = %s WHERE id = %s",
                            (nm, new_c_name, p_name, ivs_name, start_segment_type, start_segment_number, end_segment_type, end_segment_number, creation_date, var['id'])
                        )
                    log('INFO', 'Variant {0} updated to {1}'.format(var['c_name'], new_c_name))

    db.commit()


if __name__ == '__main__':
    main()
