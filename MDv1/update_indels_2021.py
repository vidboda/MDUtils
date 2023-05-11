import re
import psycopg2
import psycopg2.extras
import twobitreader
import time
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import config, md_utilities


def log(level, text):
    localtime = time.asctime(time.localtime(time.time()))
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


def get_db():
    try:
        # read connection parameters
        params = config.mdconfig()
        db = psycopg2.connect(**params)
    except (Exception, psycopg2.DatabaseError) as error:
        log('ERROR', error)
    return db


def main():
    db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    # curs.execute(
    #     "SELECT c.pos, a.id, a.c_name, a.dna_type, a.variant_size, a.wt_seq, a.mt_seq, b.chr FROM variant_feature a, gene b, variant c WHERE a.id = c.feature_id AND a.gene_name = b.name AND b.strand = '-' AND a.dna_type IN ('deletion', 'duplication' , 'indel') AND c.genome_version = 'hg38'"
    # )
    curs.execute(
        "SELECT c.pos, a.id, a.c_name, a.dna_type, a.variant_size, a.wt_seq, a.mt_seq, b.chr, b.strand FROM variant_feature a, gene b, variant c WHERE a.id = c.feature_id AND a.gene_name = b.name AND a.dna_type ='insertion' AND c.genome_version = 'hg38'"
    )
    res = curs.fetchall()
    log('INFO', f'{len(res)} variants to be checked:')
    j = 0
    for var in res:
        pos_vcf = int(var['pos'])
        if var['dna_type'] in ['indel', 'insertion']:
            pos_vcf -= 1
        x = pos_vcf - 25
        y = pos_vcf + int(var['variant_size']) + 25

        genome = twobitreader.TwoBitFile(
            f"{md_utilities.local_files['human_genome_hg38']['abs_path']}.2bit"
        )
        current_chrom = genome[f"chr{var['chr']}"]
        seq_slice = current_chrom[x:y].upper()
        # seq2 = current_chrom[int(positions[0])+1:int(positions[0])+2]
        # return seq2
        if var['strand'] == '-':
            seq_slice = md_utilities.reverse_complement(seq_slice).upper()
        marker = 25 if var['dna_type'] != 'insertion' else 26
        begin = seq_slice[:marker]
        middle = seq_slice[marker:len(seq_slice)-marker]
        end = seq_slice[-marker:]

        new_wt_seq = None
        new_mt_seq = None
        if var['dna_type'] == 'indel':
            ins_obj = re.search(r'delins([ATGC]+)', var['c_name'])
            exp_size = abs(len(middle) - len(ins_obj[1]))
            exp = ''.join('-' for _ in range(0, exp_size))
            if len(middle) > len(ins_obj[1]):
                new_wt_seq = "{0} {1} {2}".format(begin, middle, end)
                new_mt_seq = "{0} {1}{2} {3}".format(begin, ins_obj[1], exp, end)
            else:
                new_wt_seq = "{0} {1}{2} {3}".format(begin, middle, exp, end)
                new_mt_seq = "{0} {1} {2}".format(begin, ins_obj[1], end)
        elif var['dna_type'] == 'deletion':
            new_wt_seq = "{0} {1} {2}".format(begin, middle, end)
            exp = ''.join('-' for _ in range(0, var['variant_size']))
            new_mt_seq = "{0} {1} {2}".format(begin, exp, end)
        elif var['dna_type'] == 'duplication':
            exp = ''.join('-' for _ in range(0, var['variant_size']))
            new_wt_seq = "{0} {1}{2} {3}".format(begin, middle, exp, end)
            new_mt_seq = "{0} {1}{1} {2}".format(begin, middle, end)
        elif var['dna_type'] == 'insertion':
            ins_obj = re.search(r'ins([ATGC]+)', var['c_name'])
            exp = ''.join('-' for _ in range(0, len(ins_obj[1])))
            new_wt_seq = "{0} {1} {2}".format(begin, exp, end)
            new_mt_seq = "{0} {1} {2}".format(begin, ins_obj[1], end)
        if new_wt_seq != var['wt_seq'] or \
                new_mt_seq != var['wt_seq']:
            curs.execute(
                "UPDATE variant_feature SET wt_seq = %s, mt_seq = %s WHERE id = %s",
                (new_wt_seq, new_mt_seq, var['id'])
            )
            db.commit()
            log(
                'INFO',
                f"UPDATE variant_feature SET wt_seq = '{new_wt_seq}', mt_seq = '{new_mt_seq}' WHERE id = '{var['id']}'",
            )
            log('INFO', f"variant:{var['c_name']}, WT:{new_wt_seq}, MT:{new_mt_seq}")
            j += 1
    log('INFO', f'{j} variants updated / {len(res)} checked')

if __name__ == '__main__':
    main()
