import sys
import os
import argparse
import psycopg2
import psycopg2.extras
from precompute_spipv2 import get_db, log


def main():
    parser = argparse.ArgumentParser(description='Insert gnomAD data into MD',
                                     usage='python insert_gnomad.py [-d path/to/dir/containing/gnomad.v2.1.1.lof_metrics.by_gene.txt]')
    parser.add_argument('-d', '--directory', default='', required=True,
                        help='Path to the directory containing the gnomAD metrics by gene file')
    args = parser.parse_args()
    # get file
    if os.path.isfile(args.directory):
        gnomadFile = args.directory
    else:
        sys.exit('Invalid input path, please check your command')

    # get db connector and cursor
    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)
    i = 0

    for geneLine in open(gnomadFile).readlines():
        geneLineList = geneLine.split("\t")
        # print(geneLineList[0])
        curs.execute(  # exists in MD?
            """
            SELECT gene_symbol, refseq
            FROM gene
            WHERE gene_symbol = %s
                AND canonical = 't'
            """,
            (geneLineList[0],)
            # number_of_exons IN (SELECT MAX(number_of_exons) FROM gene WHERE name[1] = '{0}')".format(geneLineList[0])
        )
        mdNMFirst = curs.fetchone()
        if mdNMFirst is not None:
            # print(mdNMFirst['nm'])
            curs.execute(
                """
                SELECT DISTINCT(refseq)
                FROM gene_annotation
                WHERE gene_symbol = %s
                """,
                (geneLineList[0],)
            )  # exists in table gene_annotation? get a nm
            mdNMSecond = curs.fetchone()
            if mdNMSecond is None:
                # does not exists => creation
                i += 1
                # postGene = '{"' + mdNMFirst['gene_symbol'] + '","' + mdNMFirst['refseq'] + '"}'
                oeValues = {
                    'gene_symbol': mdNMFirst['gene_symbol'],
                    'refseq': mdNMFirst['refseq'],
                    'synoe': geneLineList[13],
                    'synlower': geneLineList[24],
                    'synupper': geneLineList[25],
                    'misoe': geneLineList[4],
                    'mislower': geneLineList[26],
                    'misupper': geneLineList[27],
                    'lofoe': geneLineList[23],
                    'loflower': geneLineList[28],
                    'lofupper': geneLineList[29]
                }
                for oeval in oeValues:
                    try:
                        oeValues[oeval] = float(oeValues[oeval])
                        oeValues[oeval] = "{:.2f}".format(oeValues[oeval])
                    except Exception:
                        next
                s = ", "
                t = "', '"
                curs.execute(
                    """
                    INSERT INTO gene_annotation ({0})
                    VALUES ('{1}')
                    """.format(
                        s.join(oeValues.keys()),
                        t.join(map(str, oeValues.values()))
                    ).replace("'NULL'", "NULL")
                )
                # curs.execute(
                #     """
                #     INSERT INTO gene_annotation
                #     VALUES('{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}')
                #     """.format(
                #         postGene,
                #         oeValues['synoe'],
                #         oeValues['synlower'],
                #         oeValues['synupper'],
                #         oeValues['misoe'],
                #         oeValues['mislower'],
                #         oeValues['misupper'],
                #         oeValues['lofoe'],
                #         oeValues['loflower'],
                #         oeValues['lofupper']
                #     )
                # )

    log('INFO', '{} annotations added'.format(i))

    db.commit()
    db_pool.putconn(db)


if __name__ == '__main__':
    main()
