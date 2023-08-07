# import os
import sys
import re
# import time
# import json
# import urllib3
import tabix
# import certifi
import psycopg2
import psycopg2.extras
from precompute_spipv2 import get_db, log
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import md_utilities, create_app

# Script to check the last and second last clinvar file against a list of variants
# returns significant changes, e.g. VUS to likely pathogenic
# can be croned after a clinvar update

result_dict = {}
clinvar_url = 'https://www.ncbi.nlm.nih.gov/clinvar'

def search_clinsig(clinvar_list):
    # retrieves CLNSIG from clinvar VCF
    match_object = re.search(r'CLNSIG=([\w\/\|]+);CLNSIGCONF=', clinvar_list[7])
    if match_object:
         return match_object.group(1)
    match_object = re.search(r'CLNSIG=([\w\/\|]+);CLNVC=', clinvar_list[7])
    if match_object:
        # log('DEBUG', clinvar_last[7])
        return match_object.group(1)
    match_object = re.search(r'CLNREVSTAT=no_interpretation_for_the_single_variant', clinvar_list[7])
    if match_object:
        return
    else:
        log('WARNING', 'Bad format for clinvar_list field: {0}'.format(clinvar_list))


def getAlleleID(clinvar_list):
    # retrieves ALLELEID from clinvar VCF
    match_object = re.search(r'ALLELEID=(\d+);CLNDISDB=', clinvar_list[7])
    if match_object:
         return match_object.group(1)


def fill_table(clinvar_last, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, generic_clinsig):
    # builds a dict of type dict{user: {variant_id: {clinvar2nd_last_version:clinsig,clinvar_last_version:clinsig}}}
    if re.search(rf'{generic_clinsig}', clinsig2nd_last) and \
            not re.search(rf'{generic_clinsig}', clinsig_last):
        vcf_str = '{0}_{1}_{2}_{3}'.format(
            var['chr'],
            var['pos'],
            var['pos_ref'],
            var['pos_alt'],
        )
        if var['mobiuser_id'] in result_dict and \
                vcf_str in result_dict[var['mobiuser_id']]:
            return
        result_dict.setdefault(var['mobiuser_id'], {'username': var['username'], 'email': var['email']})[vcf_str] = {
            clinvar2nd_last_version: clinsig2nd_last,
            clinvar_last_version: clinsig_last,
            'refseq': var['refseq'],
            'gene_symbol': var['gene_symbol'],
            'c_name': var['c_name'],
            'p_name': var['p_name'],
            'clinvar_allele_id': getAlleleID(clinvar_last),
            'mobidetails_id': var['feature_id']
        }


def get_value_from_tabix_file(tb, var):
    query = "{0}:{1}-{2}".format(var['chr'], var['pos'], var['pos'])
    try:
        records = tb.querys(query)
    except Exception as e:
        log('WARNING', 'Tabix failed for {0} with {1}'.format(query, e.args))
    for record in records:
        ref_list = re.split(',', record[3])
        alt_list = re.split(',', record[4])
        # validate nucleotides
        if var['pos_ref'] in ref_list and \
                var['pos_alt'] in alt_list:
            return record
    return None


def main():
    # get last watch versions
    last_watched_file = open("clinvar_watch_last.txt", "r")
    last_watched_version = last_watched_file.read()
    last_watched_file.close()
    log('DEBUG', 'Clinvar last watch: {0}'.format(last_watched_version))
    # get current clinvar file and version
    clinvar_last_file = md_utilities.local_files['clinvar_hg38']['abs_path']
    clinvar_last_version = md_utilities.clinvar_version
    log('DEBUG', 'Clinvar last file: {0}'.format(clinvar_last_file))
    log('DEBUG', 'Clinvar last version: {0}'.format(clinvar_last_version))
    # if no clinvar update, just leave, else print new date in watchfile
    if clinvar_last_version == last_watched_version:
        log('INFO', 'No clinvar update since last watch. Exiting.')
        sys.exit(0)
    else:
        last_watched_file = open("clinvar_watch_last.txt", "w")
        last_watched_file.write(clinvar_last_version)
        last_watched_file.close()
    # get second last clinvar file and version
    clinvar2nd_last_version = md_utilities.get_resource_current_version(
        '{0}{1}'.format(md_utilities.app_path, md_utilities.local_files['clinvar_hg38']['rel_path']),
        rf'clinvar_(\d+).vcf.gz',
        clinvar_last_version
    )
    # for testing purpose
    # clinvar2nd_last_version = '20220103'
    clinvar2nd_last_file = '{0}{1}clinvar_{2}.vcf.gz'.format(
        md_utilities.app_path,
        md_utilities.local_files['clinvar_hg38']['rel_path'],
        clinvar2nd_last_version
    )
    log('DEBUG', 'Clinvar 2nd last file: {0}'.format(clinvar2nd_last_file))
    log('DEBUG', 'Clinvar 2nd last version: {0}'.format(clinvar2nd_last_version))

    # get the list of variants
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    curs.execute(  # get transcripts
        """
        SELECT a.mobiuser_id, a.feature_id,
            b.refseq, b.c_name, b.p_name, b.gene_symbol,
            c.*, d.*
        FROM mobiuser_favourite a, variant_feature b, variant c, mobiuser d
        WHERE a.feature_id = b.id
            AND b.id = c.feature_id
            AND a.mobiuser_id = d.id
            AND d.clinvar_check = 't'
            AND a.type IN (2,3)
            AND c.genome_version = 'hg38'
        """
    )
    vars = curs.fetchall()
    log('DEBUG', '# of variants to be considered: {0}'.format(len(vars)))
    app = create_app()
    i = 0
    num_vars = curs.rowcount
    tb_last = tabix.open(clinvar_last_file)
    tb_2nd_last = tabix.open(clinvar2nd_last_file)
    for var in vars:
        print('.', end="", flush=True)
        i += 1
        if i % 1000 == 0:
            log('INFO', '{0}/{1} variants checked'.format(i, num_vars))
        # check clinvar files
        clinvar_last = get_value_from_tabix_file(tb_last, var)
        if not clinvar_last:
            # no match, next
            continue
        clinsig_last = search_clinsig(clinvar_last)
        clinvar2nd_last = get_value_from_tabix_file(tb_2nd_last, var)
        clinsig2nd_last = search_clinsig(clinvar2nd_last) if clinvar2nd_last else None
        # log('DEBUG', '{0}-{1}-{2}-{3}: last clinsig: {4} - 2nd last clinsig: {5}'.format(
        #         var['chr'],
        #         var['pos'],
        #         var['pos_ref'],
        #         var['pos_alt'],
        #         clinsig_last,
        #         clinsig2nd_last
        #     )
        # )
        # use cases:
        # {B, LB} becomes sthg else
        # {P, LP} becomes sthg else
        # VUS becomes sthg else
        # Conflicting_interpretations_of_pathogenicity becomes sthg
        # nothing becomes sthg
        if clinsig_last and \
                clinsig2nd_last:
            # {B, LB} becomes sthg else
            fill_table(clinvar_last, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, 'enign')
            # {P, LP} becomes sthg else
            fill_table(clinvar_last, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, 'athogenic')
            # VUS becomes sthg else
            fill_table(clinvar_last, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, 'Uncertain_significance')
            # Conflicting_interpretations_of_pathogenicity becomes sthg
            fill_table(clinvar_last, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, 'Conflicting_interpretations_of_pathogenicity')
        elif clinsig_last and \
                not clinsig2nd_last:
            clinsig2nd_last = 'Not recorded'
            # nothing becomes sthg
            fill_table(clinvar_last, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, 'recorded')
    db_pool.putconn(db)
    if result_dict:
        with app.app_context():
            for mobiuser in result_dict:
                email_html = """
                <p>Dear {0},</p>
                <p>Please find below a table summarizing the last changes in ClinVar concerning the variants you are watching:</p><br/>
                <table border="1" frame="hsides" rules="rows">
                    <tr>
                        <th rowspan="2" scope="col">Variant</th>
                        <th colspan="2" scope="col">ClinVar release / interpretation</th>
                    </tr>
                    <tr>
                        <th scope="col">{1}</th>
                        <th scope="col">{2}</th>
                    </tr>
                """.format(
                    result_dict[mobiuser]['username'],
                    clinvar2nd_last_version,
                    clinvar_last_version
                )
                # add variants
                for feature in result_dict[mobiuser]:
                    # search for vcf_str
                    if re.search(r'^[\dXY]{1,2}_\d+_[ATGC]+_[ATGC]+$', feature):
                        # variant
                        variant = result_dict[mobiuser][feature]
                        email_html = email_html + """
                        <tr>
                            <td><a href='https://mobidetails.iurc.montp.inserm.fr/MD/api/variant/{0}/browser/'>{1}({2}):c.{3} - {4}</a></td>
                            <td>{5}</td>
                            <td><a href='{6}?term=({7}[AlleleID])'>{8}</a></td>
                        </tr>
                        """.format(
                            variant['mobidetails_id'],
                            variant['refseq'],
                            variant['gene_symbol'],
                            variant['c_name'],
                            variant['p_name'],
                            variant[clinvar2nd_last_version],
                            clinvar_url,
                            variant['clinvar_allele_id'],
                            variant[clinvar_last_version]
                        )
                # finalize email
                email_html = email_html + """
                            </table>
                            <p>You can <a href='https://mobidetails.iurc.montp.inserm.fr/MD/auth/login' target='_blank'>connect</a> to modify your Clinvar watch settings or modify your list of followed variants.</p>
                            """
                
                md_utilities.send_email(
                    md_utilities.prepare_email_html(
                        'MobiDetails Clinvar watch',
                        email_html,
                        False
                    ),
                    '[MobiDetails - Clinvar watch]',
                    [result_dict[mobiuser]['email']]
                )
    

if __name__ == '__main__':
    main()