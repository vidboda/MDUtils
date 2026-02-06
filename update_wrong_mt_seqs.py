import re
import psycopg2
import psycopg2.extras
from precompute_spipv2 import get_db, log
# requires MobiDetails config module + database.ini file
from MobiDetailsApp import md_utilities


def main():
    # script meant to update mt seq that from 20251201, 
    # for variants lying in genes in minus strand, mutant sequence was wrong
    # need to be ran on prod
    # concerns all subsitutions when wt_seq = mt_seq and
    # all delins, ins and dups lying in minus strand genes
    # created after 20251201 by vcf_str method (no means to check this)
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)

    # get wrong substitutions
    curs.execute(
        """
        SELECT *
        FROM variant_feature
        WHERE wt_seq = mt_seq
        ORDER by creation_date
        """
    )
    subs_vars = curs.fetchall()
    for var in subs_vars:
        log("info", f"variant: {var['c_name']} - {var['refseq']}")
        # slice mt_seq
        mt_seq_obj = re.search(r'^([ATGC]+)\s([ATGC])\s([ATGC]+)$', var['mt_seq'])
        if mt_seq_obj:
            # log("debug", f"variant: {var['c_name']} - {var['refseq']}")
            begin = mt_seq_obj.group(1)
            middle = mt_seq_obj.group(2)
            end = mt_seq_obj.group(3)
            mt_sub_obj = re.search(r'>([ATGC])$', var['c_name'])
            if mt_sub_obj:
                if middle != mt_sub_obj.group(1):
                    # change mt_seq
                    new_mt_seq = f"{begin} {mt_sub_obj.group(1)} {end}"
                    log("info", f"UPDATE variant_feature SET mt_seq = {new_mt_seq} WHERE id = {var['id']}")
                    curs.execute(
                        """
                        UPDATE variant_feature
                        SET mt_seq = %s
                        WHERE id = %s
                        """,
                        (new_mt_seq, var['id'])
                    )
                else:
                    # change wt_seq ?
                    wt_seq_obj = re.search(r'^([ATGC]+)\s([ATGC])\s([ATGC]+)$', var['wt_seq'])
                    if wt_seq_obj:
                        begin = wt_seq_obj.group(1)
                        middle = wt_seq_obj.group(2)
                        end = wt_seq_obj.group(3)
                        wt_sub_obj = re.search(r'([ATGC])>', var['c_name'])
                        if wt_sub_obj:
                            new_wt_seq = f"{begin} {wt_sub_obj.group(1)} {end}"
                            log("info", f"UPDATE variant_feature SET wt_seq = {new_wt_seq} WHERE id = {var['id']}")
                            curs.execute(
                                """
                                UPDATE variant_feature
                                SET wt_seq = %s
                                WHERE id = %s
                                """,
                                (new_wt_seq, var['id'])
                            )
    db.commit()
    # correct ins
    curs.execute(
        """
        SELECT *
        FROM variant_feature a, gene b
        WHERE a.refseq = b.refseq
            AND a.dna_type IN ('insertion', 'indel')
            AND a.creation_date >= '2025-12-01'
            AND b.strand = '-'
        """
    )
    other_vars = curs.fetchall()
    for var in other_vars:
        log("info", f"variant: {var['c_name']} - {var['refseq']}")
        # slice mt_seq
        mt_seq_obj = re.search(r'^([ATGC]+)\s([ATGC]+)-*\s([ATGC]+)$', var['mt_seq'])
        if mt_seq_obj:
            # log("debug", f"variant: {var['c_name']} - {var['refseq']}")
            begin = mt_seq_obj.group(1)
            middle = mt_seq_obj.group(2)
            end = mt_seq_obj.group(3)
            mt_sub_obj = re.search(r'ins([ATGC]+)$', var['c_name'])
            if mt_sub_obj.group(1) != middle:
                # wrong
                new_mt_seq = f"{begin} {mt_sub_obj.group(1)} {end}"
                log("info", f"UPDATE variant_feature SET mt_seq = {new_mt_seq} WHERE id = {var['id']}")
                curs.execute(
                    """
                    UPDATE variant_feature
                    SET mt_seq = %s
                    WHERE id = %s
                    """,
                    (new_mt_seq, var['id'])
                )
    db.commit()
    db_pool.putconn(db)


if __name__ == '__main__':
    main()