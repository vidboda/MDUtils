import os
import sys
import re
import json
import psycopg2
import psycopg2.extras
import time
from insert_genes import get_db
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config, md_utilities

# script meant to be used after switching the exon numbering system from MD to VV (genes json files)
# the script will "remap" all variants


def log(level, text):
    print()
    localtime = time.asctime(time.localtime(time.time()))
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


def main():
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    curs.execute(  # get all variants
        "SELECT a.c_name, a.gene_name, a.ivs_name, a.start_segment_type, a.start_segment_number, a.end_segment_type, a.end_segment_number, b.chr, b.nm_version, c.g_name FROM variant_feature a, gene b, variant c WHERE a.gene_name = b.name AND a.id = c.feature_id AND c.genome_version = 'hg38' ORDER BY random() LIMIT 150"
    )

    vars = curs.fetchall()
    for var in vars:
        positions = md_utilities.compute_start_end_pos(var['g_name'])
        ncbi_chr = md_utilities.get_ncbi_chr_name(db, f"chr{var['chr']}", 'hg38')
        # read VV json file
        try:
            json_file = open('{0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                var['gene_name'][0]
            ))
        except IOError:
            # print('file_not_found_error')
            log('WARNING','VV file_not_found_error for {0}'.format(var['gene_name'][0]))
        vv_json = json.load(json_file)
        if 'transcripts' in vv_json:
            for vv_transcript in vv_json['transcripts']:
                if (
                    vv_transcript['reference']
                    == '{0}.{1}'.format(var['gene_name'][1], var['nm_version'])
                    and ncbi_chr[0] in vv_transcript['genomic_spans']
                ):
                    # we're in
                    # we loop on exons
                    last_segment = {
                        'genomic_end': 10000000000,
                        'genomic_start': 0,
                        'exon_number': 0
                    }
                    update_start = update_end = 0
                    for exon in vv_transcript['genomic_spans'][ncbi_chr[0]]['exon_structure']:
                        # exons
                        if int(positions[0]) >= int(exon['genomic_start']) and \
                                int(positions[0]) <= int(exon['genomic_end']):
                            update_var_segment('start', 'exon', exon['exon_number'], var, curs)
                            update_start = 1
                        if int(positions[1]) >= int(exon['genomic_start']) and \
                                int(positions[1]) <= int(exon['genomic_end']):
                            update_var_segment('end', 'exon', exon['exon_number'], var, curs)
                            update_end = 1
                        # introns
                        if int(positions[0]) < int(exon['genomic_start']) and \
                                int(positions[0]) > int(last_segment['genomic_end']):
                            update_var_segment('start', 'intron', last_segment['exon_number'], var, curs)
                            update_start = 1
                        if int(positions[1]) < int(exon['genomic_start']) and \
                                int(positions[1]) > int(last_segment['genomic_end']):
                            update_var_segment('end', 'intron', last_segment['exon_number'], var, curs)
                            update_end = 1
                        #
                        if update_start == 1 and \
                                update_end == 1:
                            break
        else:
            log('ERROR', 'Error for {0}: {1}'.format(var['gene_name'][0], vv_json))


def update_var_segment(segment, segment_type, segment_number, var, curs):
    if (
        segment != 'start'
        or segment_type == var['start_segment_type']
        and segment_number == var['start_segment_number']
    ) and (
        segment != 'end'
        or segment_type == var['end_segment_type']
        and segment_number == var['end_segment_number']
    ):
        return
        # IVS name
    if segment_type == 'intron' and \
                segment == 'end':
        curs.execute(
            "SELECT start_segment_number, end_segment_number FROM variant_feature WHERE id = %s",
            (var['id'],)
        )
        if segments := curs.fetchone():
            ivs_name = None
            if ivs_obj := re.search(
                r'^[\*-]?\d+([\+-]\d+)(.+)$', var['c_name']
            ):
                ivs_name = 'IVS{0}{1}{2}'.format(
                    segments['start_segment_number'], ivs_obj[1], ivs_obj[2]
                )
            elif ivs_obj := re.search(
                r'^\d+([\+-]\d+)_\d+([\+-]\d+)(.+)$', var['c_name']
            ):
                ivs_name = 'IVS{0}{1}_IVS{2}{3}{4}'.format(
                    segments['start_segment_number'],
                    ivs_obj[1],
                    segments['end_segment_number'],
                    ivs_obj[2],
                    ivs_obj[3],
                )
            elif ivs_obj := re.search(
                r'^\d+([\+-]\d+)_(\d+)([^\+-].+)$', var['c_name']
            ):
                ivs_name = 'IVS{0}{1}_{2}{3}'.format(
                    segments['start_segment_number'],
                    ivs_obj[1],
                    ivs_obj[2],
                    ivs_obj[3],
                )
            elif ivs_obj := re.search(
                r'^(\d+)_\d+([\+-]\d+)(.+)$', var['c_name']
            ):
                ivs_name = '{0}_IVS{1}{2}{3}'.format(
                    ivs_obj[1],
                    segments['end_segment_number'],
                    ivs_obj[2],
                    ivs_obj[3],
                )
            if ivs_name:
                # update MD
                # curs.execute(
                #     "UPDATE variant_feature SET ivs_name = %s WHERE id = %s",
                #     (ivs_name, var['id'])
                # )
                log('INFO', 'New IVS name for {0}.{1}({2}):c.{3}: {4}'.format(
                    var['gene_name'][1],
                    var['nm_version'],
                    var['gene_name'][0],
                    var['c_name'],
                    ivs_name
                ))
    log('INFO', 'New segments for {0}.{1}({2}):c.{3}: {4}-{5}-{6}'.format(
        var['gene_name'][1],
        var['nm_version'],
        var['gene_name'][0],
        var['c_name'],
        segment,
        segment_type,
        segment_number,
    ))
    # update MD
    segment_type_to_update = 'start_segment_type'
    segment_number_to_update = 'start_segment_number'
    if segment == 'end':
        segment_type_to_update = 'end_segment_type'
        segment_number_to_update = 'end_segment_number'
        # curs.execute(
        #     "UPDATE variant_feature SET %s = %s, %s = %s WHERE id = %s",
        #     (segment_type_to_update, segment_type, segment_number_to_update, segment_number, var['id'])
        # )
        # db.commit()

if __name__ == '__main__':
    main()
