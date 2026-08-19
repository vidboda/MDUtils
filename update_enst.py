import os
import sys
import psycopg2
import psycopg2.extras
import argparse
from precompute_spipv2 import get_db, log
from MobiDetailsApp import md_utilities


def main():
    # Updates Eensembl ENST according to https://ftp.ensembl.org/pub/release-116/tsv/homo_sapiens/Homo_sapiens.GRCh38.116.refseq.tsv.gz
    # and UNIPROT IDs according to https://ftp.ensembl.org/pub/release-116/tsv/homo_sapiens/Homo_sapiens.GRCh38.116.uniprot.tsv.gz
    # both are merged in resources/ens2refseq2uniprot.GRCh38.116.tsv like
    # transcript_stable_id    refseq_xref     uniprot_xref
    # transcript_stable_id    refseq_xref     uniprot_xref
    # ENST00000420522 NR_147025
    # ENST00000443270 NR_183671
    # ENST00000417651 NR_110758
    # ENST00000376030 NM_201628       Q674X7
    # ENST00000503743 NM_001370230    Q674X7
    # ENST00000503743 NM_001437721    Q674X7
    # ENST00000503743 NM_015209       Q674X7
    # ENST00000361144 NM_001018000    Q674X7
    # ENST00000400797 NM_001017999    Q674X7
    # ENST00000400798 NM_001018001    Q674X7
    # ENST00000400798 NM_001370229    Q674X7
    # ENST00000400798 NM_001370231    Q674X7
    # ENST00000404665 NR_027136
    # ENST00000489803 NR_037633
    # ENST00000454994 NR_037634       C9IYK1
    # ENST00000621281 NR_037635       A0A087WV05
    # ENST00000619609 NM_207421       Q6TGC4
    # read the file line per line
    # check whether the ENST is already in the database, if not insert it, if yes update it
    # check whether the refseq has a uniprot ID associated, if yes update it if necessary, if no insert - DEPRECATED UNIPROT IDS ARE MANAGED ELSEWHERE
    parser = argparse.ArgumentParser(description='Manage ENST data', usage='python update_enst.py [--dry-run]')
    parser.add_argument('-d', '--dry-run', default='', required=False, help='Show updates but do not commit to the database', action='store_true')
    args = parser.parse_args()
    # initiailize the database connection
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    # read the file
    i = 0
    with open("resources/ens2refseq2uniprot.GRCh38.116.tsv", "r") as f:
        for line in f:
            # stop after 1000 lines
            i += 1
            if i % 1000 == 0:
                log('INFO', '{0} transcripts checked'.format(i))
                if not args.dry_run:
                    db.commit()
                # break
            # skip the header
            if line.startswith("transcript_stable_id"):
                continue
            # log("DEBUG", f"Processing line: {line.strip()}")
            # split the line by tab
            parts = line.strip().split("\t")
            enst = parts[0]
            refseq = parts[1]
            uniprot = parts[2] if len(parts) > 2 else None
            # check whether the ENST is already in the database
            curs.execute("SELECT refseq, enst FROM gene WHERE refseq ~ CONCAT(%s, '\.\d{1,2}')", (refseq,))
            result = curs.fetchone()
            # ENST
            if not result:
                # log("WARNING", f"RefSeq {refseq} not found in database")
                continue
            if not result['enst']:
                # here we add only lacking ENSTs, we do not update existing ones
                # insert the ENST
                curs.execute("UPDATE gene set enst = %s WHERE refseq = %s", (enst, result['refseq']))
                # if result['enst'] is not None:
                log("INFO", f"Inserted ENST {enst} for RefSeq {result['refseq']} old ENST: {result['enst']}")
            else:
                if result['enst'] != enst:
                    curs.execute("UPDATE gene set enst = %s WHERE refseq = %s AND enst = %s", (enst, result['refseq'], result['enst']))
                    log("INFO", f"RefSeq {refseq} already has ENST {result['enst']} in database, new ENST {enst} inserted")
            # UNIPROT - treated in check_uniprot_ids.py
            # if result["uniprot_id"] != uniprot:
            #     # first check whether the UNIPROT is already in the database
            #     curs.execute("SELECT id FROM uniprot WHERE id = %s", (uniprot,))
            #     uniprot_result = curs.fetchone()
            #     if uniprot_result is None:
            #         # insert the UNIPROT
            #         # curs.execute("INSERT INTO uniprot (id) VALUES (%s)", (uniprot,))
            #         log("WARNING", f"UniProt ID to be inserted {uniprot}")
            #     # insert the UNIPROT
            #     # curs.execute("UPDATE gene set uniprot_id = %s WHERE refseq = %s", (uniprot, result['refseq']))
            #     if result['uniprot_id'] is not None:
            #         log("INFO", f"Inserted UniProt {uniprot} for RefSeq {result['refseq']}")
    if not args.dry_run:
        db.commit()
    db_pool.putconn(db)


if __name__ == '__main__':
    main()