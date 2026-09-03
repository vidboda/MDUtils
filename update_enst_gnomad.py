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
                 'old_gnomad_occurrences', 'new_gnomad_occurrences',
                 'num_candidates', 'reason']

MIN_GNOMAD_DEFAULT = 2

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

def choose_enst(candidates, gnomad_counts, min_count, refseq):
    """candidates: list of (enst, uniprot) tuples.
    Returns (chosen_enst, decision_source, reason)."""
    if len(candidates) == 1:
        return candidates[0][0], 'unique', 'single ENST candidate for this refseq'

    # Filter to ENSTs meeting the minimum gnomAD occurrence threshold
    hits = [(gnomad_counts.get(enst, 0), enst, uni) for enst, uni in candidates]
    seen = [h for h in hits if h[0] >= min_count]

    if len(seen) == 1:
        return seen[0][1], 'gnomad', f'only one ENST candidate meets min-gnomad-count ({min_count})'
    if len(seen) > 1:
        seen.sort(key=lambda h: (-h[0], not is_reviewed_uniprot(h[2]), h[1]))
        log('DEBUG', f'{refseq}: {len(seen)} ENSTs meet min-gnomad-count, chose {seen[0][1]} '
                     f'(candidates: {[(h[1], h[0]) for h in seen]})')
        return (seen[0][1], 'gnomad',
                f'{len(seen)} ENSTs meet min-gnomad-count, kept the most frequent')

    # no candidate ENST meets threshold -> fallback
    fallback = sorted(candidates, key=lambda c: (not is_reviewed_uniprot(c[1]), c[0]))
    log('WARNING', f'{refseq}: no candidate ENST meets min-gnomad-count ({min_count}), '
                   f'falling back to {fallback[0][0]}')
    return (fallback[0][0], 'fallback',
            'no candidate ENST meets threshold, kept reviewed UniProt / first sorted')

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
    parser.add_argument('--min-gnomad-count', type=int, default=MIN_GNOMAD_DEFAULT,
                        help=f'Minimum gnomAD BCSQ occurrences for an ENST to be considered '
                             f'"active" (default: {MIN_GNOMAD_DEFAULT})')
    parser.add_argument('--test-limit', type=int, default=None,
                        help='Test mode: process only N refseqs for validation (e.g., 10000)')
    args = parser.parse_args()

    db_pool, db = get_db()
    curs = db.cursor(cursor_factory=psycopg2.extras.DictCursor)

    # TEST MODE WARNING
    if args.test_limit:
        log('INFO', f'TEST MODE: will process only {args.test_limit} refseqs for validation')

    # 1. Load ENSTs annotated in gnomAD
    gnomad_counts = load_gnomad_ensts(args.vcf)

    # 2. Collect ENST candidates per refseq from the Ensembl TSV
    candidates = defaultdict(list)
    total_refseq = 0
    with open(args.tsv, 'r') as f:
        for line in f:
            if line.startswith('transcript_stable_id'):
                continue
            total_refseq += 1
            if args.test_limit and total_refseq > args.test_limit:
                log('INFO', f'Test limit reached at {total_refseq} refseqs')
                break
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            enst, refseq = parts[0], parts[1]
            uniprot = parts[2] if len(parts) > 2 else None
            candidates[refseq].append((enst, uniprot))

    actual_processed = sum(1 for _ in candidates)
    log('INFO', f'Loaded {actual_processed} refseq candidates '
                f'({"limited by test-limit" if args.test_limit else ""})')

    # 3. Resolve choices and update the database,
    #    recording every change in the CSV report
    n_updated = 0
    n_kept = 0
    n_skipped_no_change = 0
    with open(args.report, 'w', newline='') as report_f:
        writer = csv.DictWriter(report_f, fieldnames=REPORT_FIELDS)
        writer.writeheader()

        for i, (refseq, cands) in enumerate(sorted(candidates.items()), start=1):
            chosen, source, reason = choose_enst(cands, gnomad_counts,
                                                args.min_gnomad_count, refseq)
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
                    n_skipped_no_change += 1
                    continue

                # SAFEGUARD: keep the old ENST if it is annotated in gnomAD
                # while the newly chosen ENST is not or less annotated.
                # Switching would lose gnomAD annotation coverage.
                old_is_active = gnomad_counts.get(old_enst, 0) >= args.min_gnomad_count
                # new_is_active = gnomad_counts.get(chosen, 0) >= gnomad_counts.get(old_enst, 0)
                new_is_active = gnomad_counts.get(chosen, 0) >= args.min_gnomad_count

                if old_is_active and not new_is_active:
                    n_kept += 1
                    writer.writerow({
                        'refseq': res['refseq'],
                        'old_enst': old_enst,
                        'new_enst': '',
                        'decision_source': 'kept_old_gnomad',
                        'old_gnomad_occurrences': gnomad_counts.get(old_enst, 0),
                        'new_gnomad_occurrences': gnomad_counts.get(chosen, 0),
                        'num_candidates': len(cands),
                        'reason': 'old ENST is gnomAD-annotated (>=threshold), new ENST is not',
                    })
                    log('INFO', f"{res['refseq']}: kept ENST {old_enst} "
                                f"(gnomAD-annotated, {gnomad_counts.get(old_enst, 0)} >= {gnomad_counts.get(chosen, 0)}), "
                                f"not switching to {chosen} ({gnomad_counts.get(chosen, 0)})")
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
                    'old_gnomad_occurrences': gnomad_counts.get(old_enst, 0),
                    'new_gnomad_occurrences': gnomad_counts.get(chosen, 0),
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

    summary = (f'{n_updated} gene rows updated, '
               f'{n_kept} kept unchanged (old ENST gnomAD-protected), '
               f'{n_skipped_no_change} skipped (no change needed) '
               f'{"(DRY RUN - nothing committed)" if args.dry_run else ""} '
               f'{"(TEST MODE)" if args.test_limit else ""} '
               f'- report written to {args.report}')
    log('INFO', f'Done: {summary}')
    db_pool.putconn(db)

if __name__ == '__main__':
    main()