import os
import time
import json
import urllib3
import certifi
import psycopg2
import psycopg2.extras
from precompute_spipv2 import get_db, log
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import md_utilities

# Script to check whether the transcripts recorded in MD are present in VV json files
# If a transcript is deactivated and found in a new file, becomes ok
# can be ran with each VVTA update


def download_vv_file(gene, transcript):
    vv_json_gene_file = '{0}{1}.json'.format(
        md_utilities.local_files['variant_validator']['abs_path'],
        gene
    )
    # check if the file has already been modified today
    if ((time.time() - os.path.getmtime(vv_json_gene_file)) / 3600) > 24:
        vv_url_base = "https://rest.variantvalidator.org"
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        vv_url = "{0}/VariantValidator/tools/gene2transcripts/{1}?content-type=application/json".format(vv_url_base, transcript)
        try:
            vv_data = json.loads(http.request('GET', vv_url).data.decode('utf-8'))
            if 'transcripts' in vv_data:
                with open(
                    vv_json_gene_file,
                    "w",
                    encoding='utf-8'
                ) as vv_file:
                    json.dump(
                        vv_data,
                        vv_file,
                        ensure_ascii=False,
                        indent=4
                    )
                log('INFO', "VV JSON file copied for gene {0}-{1}".format(gene, transcript))
            else:
                log('WARNING', 'No transcript in VV JSON file for {0}'.format(gene))
        except Exception:
            log('WARNING', 'No VV JSON file for {0}'.format(gene))


def check_vv_file(gene, db, curs, ncbi_chr, ncbi_chr_hg19):
    if os.path.isfile(
        '{0}{1}.json'.format(
            md_utilities.local_files['variant_validator']['abs_path'],
            gene['name'][0]
        )
    ):
        vv_file = open('{0}{1}.json'.format(
            md_utilities.local_files['variant_validator']['abs_path'],
            gene['name'][0]
        ), 'rb')
        vv_json = json.load(vv_file)
        transcript_found = False
        try:
            transcript_found, ncbi_chr, ncbi_chr_hg19 = check_vv_transcript(transcript_found, gene, vv_json, ncbi_chr, ncbi_chr_hg19, curs, db)
            # if not transcript_found:
            #     # transcript needs to be deactivated
            #     log('WARNING', "UPDATE gene SET variant_creation = 'not_in_vv_json' WHERE  name[2] = {0}".format(gene['name'][1]))
            #     curs.execute(
            #         """
            #         UPDATE gene
            #         SET variant_creation = 'not_in_vv_json'
            #         WHERE  name[2] = %s
            #         """,
            #         (gene['name'][1],)
            #     )
            #     db.commit()
        except KeyError:
            log('WARNING', 'No transcript in file {0}{1}'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                gene['name'][0]
            ))
    return ncbi_chr, ncbi_chr_hg19


def check_vv_transcript(transcript_found, gene, vv_json, ncbi_chr, ncbi_chr_hg19, curs, db):
    for vv_transcript in vv_json['transcripts']:
        # log('DEBUG', '{0}-{1}'.format(vv_transcript['reference'], gene['name'][1]))
        if vv_transcript['reference'] == gene['refseq']:
            chrom = vv_transcript['annotations']['chromosome']
            if isinstance(chrom, list):
                chrom = chrom[0]
            if chrom not in ncbi_chr:
                curs.execute(
                    """
                    SELECT ncbi_name, genome_version
                    FROM chromosomes WHERE name = %s
                    """,
                    (chrom,)
                )
                ncbi_name = curs.fetchall()
                for chromo in ncbi_name:
                    if chromo['genome_version'] == 'hg19':
                        ncbi_chr_hg19[chrom] = chromo['ncbi_name']
                    elif chromo['genome_version'] == 'hg38':
                        ncbi_chr[chrom] = chromo['ncbi_name']
            # chek that we have hg19 and hg38
            if ncbi_chr[chrom] in vv_transcript['genomic_spans'] and \
                    ncbi_chr_hg19[chrom] in vv_transcript['genomic_spans']:
                transcript_found = True
                # check hg19 and hg38
                if gene['variant_creation'] != 'ok':
                    curs.execute(
                        """
                        UPDATE gene
                        SET variant_creation = 'ok'
                        WHERE  refseq = %s
                        """,
                        (gene['refseq'],)
                    )
                    db.commit()
                    log('INFO', 'Transcript {0} becomes ok'.format(gene['refseq']))
                return transcript_found, ncbi_chr, ncbi_chr_hg19
            else:
                # hg19/38 mapping issue
                log('WARNING', 'Transcript {0} from gene {1} has hg19/38 mapping issues'.format(gene['refseq'], gene['gene_symbol']))
                default = 'hg19_mapping_default'
                if ncbi_chr[chrom] not in vv_transcript['genomic_spans'] and \
                        ncbi_chr_hg19[chrom] not in vv_transcript['genomic_spans']:
                    default = 'mapping_default'
                elif ncbi_chr[chrom] not in vv_transcript['genomic_spans'] and \
                        ncbi_chr_hg19[chrom] in vv_transcript['genomic_spans']:
                    default = 'hg38_mapping_default'
                curs.execute(
                    """
                    UPDATE gene
                    SET variant_creation = %s
                    WHERE refseq = %s
                    """,
                    (default, gene['refseq'])
                )
                db.commit()
                return transcript_found, ncbi_chr, ncbi_chr_hg19
    # transcript not in file
    log('WARNING', 'Transcript {0} from gene {1} not found in VV json'.format(gene['refseq'], gene['gene_symbol']))
    if gene['variant_creation'] == 'ok':
        curs.execute(
            """
            UPDATE gene
            SET variant_creation = 'not_in_vv_json'
            WHERE refseq = %s
            """,
            (gene['refseq'],)
        )
        db.commit()
    return transcript_found, ncbi_chr, ncbi_chr_hg19


def main():
    # script which checks whether a transcript in MD is present in VV gene jsons
    # if not should be deactivated
    ncbi_chr = {}
    ncbi_chr_hg19 = {}
    # 1st get all transcripts
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    curs.execute(  # get transcripts
        """
        SELECT gene_symbol, refseq, variant_creation
        FROM gene
        ORDER BY name[1]
        """
    )
    genes = curs.fetchall()
    for gene in genes:
        # get VV json file
        if gene['gene_symbol'] == '':
            log('ERROR', 'A gene has no name!!!')
        if os.path.isfile(
            '{0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                gene['gene_symbol']
            )
        ):
            vv_file = open('{0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                gene['gene_symbol']
            ), 'rb')
            vv_json = json.load(vv_file)
            transcript_found = False
            try:
                transcript_found, ncbi_chr, ncbi_chr_hg19 = check_vv_transcript(transcript_found, gene, vv_json, ncbi_chr, ncbi_chr_hg19, curs, db)
                if not transcript_found:
                    download_vv_file(gene['gene_symbol'], gene['refseq'])
                    # check if downloaded file is ok?
                    # neither error nor "transcripts": []
                    # if error mark gene as unusable in MD
                    ncbi_chr, ncbi_chr_hg19 = check_vv_file(gene, db, curs, ncbi_chr, ncbi_chr_hg19)
            except KeyError:
                log('WARNING', 'No transcript in file {0}{1}'.format(
                    md_utilities.local_files['variant_validator']['abs_path'],
                    gene['gene_symbol']
                ))
                download_vv_file(gene['gene_symbol'], gene['refseq'])
                ncbi_chr, ncbi_chr_hg19 = check_vv_file(gene, db, curs, ncbi_chr, ncbi_chr_hg19)
        else:
            log('WARNING', 'File not found {0}{1}.json'.format(
                md_utilities.local_files['variant_validator']['abs_path'],
                gene['name'][0]
            ))
            download_vv_file(gene['gene_symbol'], gene['refseq'])
            ncbi_chr, ncbi_chr_hg19 = check_vv_file(gene, db, curs, ncbi_chr, ncbi_chr_hg19)
    db_pool.putconn(db)


if __name__ == '__main__':
    main()
