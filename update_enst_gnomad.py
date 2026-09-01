import re
import gzip
import csv
import argparse
from collections import defaultdict

import psycopg2
import psycopg2.extras

from precompute_spipv2 import get_db, log

ENST_RE = re.compile(r'^ENST\d+')

REPORT_FIELDS = ['refseq', 'old_enst', 'new_enst', 'decision_source',
                 'gnomad_occurrences', 'num_candidates', 'reason']


def load_gnomad_ensts(vcf_path):
    """Single pass over the gnomAD VCF: count occurrences of each ENST
    in the INFO/BCSQ field."""
    counts = defaultdict(int)
    n_var = 0
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            n_var += 1
            if n_var % 500000 == 0:
                log('INFO', f'{n_var} variants read, {len(counts)} distinct ENSTs')
            info = line.split('\t', 8)[7]
            if 'BCSQ=' not in info:
                continue
            bcsq = info.split('BCSQ=', 1)[1].split(';', 1)[0]
            for eff in bcsq.split(','):
                fields = eff.split('|')
                # BCSQ format: effect|gene|ENST|biotype|strand|aa|dna
                if len(fields) > 2 and fields[2].startswith('ENST'):
                    # strip version suffix if present (ENST00000374708.9)
                    counts[fields[2].split('.')[0]] += 1
    log('INFO', f'VCF done: {n_var} variants, {len(counts)} distinct ENSTs')
    return counts


def is_reviewed_uniprot(acc):
    """Approximate SwissProt detection: not a TrEMBL accession (A0A...)."""
    return bool(acc) and not acc.startswith('A0A')


def choose_enst(candidates, gnomad_counts, refseq):
    """candidates: list of (enst, uniprot) tuples.
    Returns (chosen_enst, decision_source, reason)."""
    if len(candidates) == 1:
        return candidates[0][0], 'unique', 'single ENST candidate for this refseq'

    hits = [(gnomad_counts.get(enst, 0), enst, uni) for enst, uni in candidates]
    seen = [h for h in hits if h[0] > 0]

    if len(seen) == 1:
        return seen[0][1], 'gnomad', 'only one ENST candidate found in gnomAD BCSQ'
    if len(seen) > 1:
        seen.sort(key=lambda h: (-h[0], not is_reviewed_uniprot(h[2]), h[1]))
        log('DEBUG', f'{refseq}: {len(seen)} ENSTs in gnomAD, chose {seen[0][1]} '
                     f'(candidates: {[(h[1], h[0]) for h in seen]})')
        return (seen[0][1], 'gnomad',
                f'{len(seen)} ENSTs annotated in gnomAD, kept the most frequent')

    # no candidate ENST found in gnomAD -> fallback
    fallback = sorted(candidates, key=lambda c: (not is_reviewed_uniprot(c[1]), c[0]))
    log('WARNING', f'{refseq}: no candidate ENST found in gnomAD, '
                   f'falling back to {fallback[0][0]}')
    return (fallback[0][0], 'fallback',
            'no candidate ENST in gnomAD, kept reviewed UniProt / first sorted')


def main():
    parser = argparse.ArgumentParser(
        description='Update ENST choices based on gnomAD BCSQ annotations',
        usage='python update_enst_gnomad.py [--dry-run]')
    parser.add_argument('--vcf', default='gnomad_ms.fully_annotated.vcf.gz')
    parser.add_argument('--tsv', default='resources/ens2refseq2uniprot.GRCh38.116.tsv')
    parser.add_argument('-r', '--report', default='enst_update_report.csv',
                        help='CSV report file listing all changes')
    parser.add_argument('-d', '--dry-run', action='store_true',
                        help='Show updates and generate the report '
                             'but do not commit to the database')
    args = parser.parse_args()

    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)

    # 1. Load ENSTs annotated in gnomAD
    gnomad_counts = load_gnomad_ensts(args.vcf)

    # 2. Collect ENST candidates per refseq from the Ensembl TSV
    candidates = defaultdict(list)
    with open(args.tsv, 'r') as f:
        for line in f:
            if line.startswith('transcript_stable_id'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            enst, refseq = parts[0], parts[1]
            uniprot = parts[2] if len(parts) > 2 else None
            candidates[refseq].append((enst, uniprot))

    # 3. Resolve choices and update the database,
    #    recording every change in the CSV report
    n_updated = 0
    with open(args.report, 'w', newline='') as report_f:
        writer = csv.DictWriter(report_f, fieldnames=REPORT_FIELDS)
        writer.writeheader()

        for i, (refseq, cands) in enumerate(sorted(candidates.items()), start=1):
            chosen, source, reason = choose_enst(cands, gnomad_counts, refseq)
            if not chosen:
                continue

            curs.execute(
                "SELECT refseq, enst FROM gene WHERE refseq ~ CONCAT(%s, '\\.\\d{1,2}')",
                (refseq,))
            rows = curs.fetchall()
            if not rows:
                continue

            for res in rows:
                old_enst = res['enst']
                if old_enst == chosen:
                    continue
                curs.execute(
                    "UPDATE gene SET enst = %s WHERE refseq = %s",
                    (chosen, res['refseq']))
                n_updated += 1
                writer.writerow({
                    'refseq': res['refseq'],
                    'old_enst': old_enst,
                    'new_enst': chosen,
                    'decision_source': source,
                    'gnomad_occurrences': gnomad_counts.get(chosen, 0),
                    'num_candidates': len(cands),
                    'reason': reason,
                })
                log('INFO', f"{res['refseq']}: ENST {old_enst} -> {chosen} "
                            f"({'gnomAD-annotated' if source == 'gnomad' else 'fallback'})")

            # periodic commit unless dry-run
            if i % 1000 == 0 and not args.dry_run:
                db.commit()
                log('INFO', f'{i} refseq processed')

    if not args.dry_run:
        db.commit()
    log('INFO', f'Done: {n_updated} gene rows updated '
                f'{"(DRY RUN - nothing committed)" if args.dry_run else ""} '
                f'- report written to {args.report}')
    db_pool.putconn(db)


if __name__ == '__main__':
    main()