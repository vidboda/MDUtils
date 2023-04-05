import os
import re
import time
import json
import urllib3
import tabix
import certifi
import psycopg2
import psycopg2.extras
from precompute_spipv2 import get_db, log
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import md_utilities, create_app

# Script to check the last and second last clinvar file against a list of variants
# returns significant changes, e.g. VUS to likely pathogenic
# can be croned after a clinvar update


def search_clinsig(clinvar_list):
    match_object = re.search(r'CLNSIG=([\w\/\|]+);CLNSIGCONF=', clinvar_list[7])
    if match_object:
         return match_object.group(1)
    match_object = re.search(r'CLNSIG=([\w\/\|]+);CLNVC=', clinvar_list[7])
    if match_object:
        # log('DEBUG', clinvar_last[7])
        return match_object.group(1)
    else:
        log('WARNING', 'Bad format for clinvar_list field: {0}'.format(clinvar_list))


def trigger_alert(app, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, generic_clinsig):
    if re.search(rf'{generic_clinsig}', clinsig2nd_last) and \
            not re.search(rf'{generic_clinsig}', clinsig_last):
        with app.app_context():
            md_utilities.send_email(
                md_utilities.prepare_email_html(
                    'MobiDetails Clinvar watch',
                    """
                    <p>Dear {0},</p>
                    <p>The Clinvar interpretation of a variant that you are watching has changed with the last release:</p>
                    <ul>
                        <li>{1}({2}):c.{3} - p.{4}: in Clinvar release {5}: {6}, becomes</li>
                        <li>{1}({2}):c.{3} - p.{4}: in Clinvar release {7}: {8}</li>
                    </ul>
                    <p>Check out this variant in <a href='https://mobidetails.iurc.montp.inserm.fr/MD/api/variant/{9}/browser/' target='_blank'>MobiDetails</a> or <a href='https://mobidetails.iurc.montp.inserm.fr/MD/auth/login' target='_blank'>connect</a> to modify your Clinvar watch settings or modify your list of followed variants.</p>
                    """.format(
                        var['username'],
                        var['refseq'],
                        var['gene_symbol'],
                        var['c_name'],
                        var['p_name'],
                        clinvar2nd_last_version,
                        clinsig2nd_last,
                        clinvar_last_version,
                        clinsig_last,
                        var['feature_id']
                    ),
                    False
                ),
                '[MobiDetails - Clinvar follow up]',
                [var['email']]
        )

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
    # get current clinvar file and version
    clinvar_last_file = md_utilities.local_files['clinvar_hg38']['abs_path']
    clinvar_last_version = md_utilities.clinvar_version
    log('DEBUG', 'Clinvar last file: {0}'.format(clinvar_last_file))
    log('DEBUG', 'Clinvar last version: {0}'.format(clinvar_last_version))
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
            trigger_alert(app, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, 'enign')
            # {P, LP} becomes sthg else
            trigger_alert(app, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, 'athogenic')
            # VUS becomes sthg else
            trigger_alert(app, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, 'Uncertain_significance')
            # Conflicting_interpretations_of_pathogenicity becomes sthg
            trigger_alert(app, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, 'Conflicting_interpretations_of_pathogenicity')
        elif clinsig_last and \
                not clinsig2nd_last:
            clinsig2nd_last = 'Not recorded'
            # nothing becomes sthg
            trigger_alert(app, var, clinsig_last, clinsig2nd_last, clinvar_last_version, clinvar2nd_last_version, 'recorded')
    db_pool.putconn(db)

if __name__ == '__main__':
    main()