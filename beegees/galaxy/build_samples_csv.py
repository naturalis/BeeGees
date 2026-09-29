"""Resolve a user-supplied samples CSV for BeeGees.

The user provides a CSV with columns ID, forward, [reverse], and optional
taxonomy columns where forward/reverse are bare filenames (not full paths).
This script maps those filenames to absolute paths using the Galaxy collection,
then writes the resolved CSV that BeeGees expects.

Usage:
    python build_samples_csv.py
        --samples-csv  /path/to/user_samples.csv
        --output       samples.csv
        --collection-names  name1 name2 ...
        --collection-paths  /abs/path1 /abs/path2 ...
        [--reads-dir   reads]
"""

import argparse
import csv
import os
import sys

REQUIRED_COLS = {
    'ID':      ['ID'],
    'forward': ['forward'],
}
OPTIONAL_READ_COL = {
    'reverse': ['reverse'],
}
TAX_COLS = ['phylum', 'class', 'order', 'family', 'genus', 'species']


def err(stage, msg):
    print(f'[build_samples_csv] ERROR ({stage}): {msg}', file=sys.stderr)
    sys.exit(1)


def resolve_col(fieldnames, aliases):
    for a in aliases:
        if a in fieldnames:
            return a
    return None


def is_gzipped(path):
    with open(path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def link_read(filename, target, reads_dir):
    """Symlink a Galaxy dataset (.dat) under its real filename so tools can detect gzip by extension."""
    name = os.path.basename(filename)
    if name in ('', '.', '..'):
        err('file resolution', f'Invalid filename: "{filename}"')
    gz = is_gzipped(target)
    if gz and not name.endswith('.gz'):
        name += '.gz'
    elif not gz and name.endswith('.gz'):
        name = name[:-3]
    link = os.path.abspath(os.path.join(reads_dir, name))
    if os.path.lexists(link):
        if os.path.realpath(link) != os.path.realpath(target):
            err('file resolution', f'Two different files map to the same name: {name}')
        return link
    os.symlink(os.path.abspath(target), link)
    return link


FASTQ_EXTS = ('.fastqsanger', '.fastq', '.fq')


def strip_fastq_ext(name):
    if name.endswith('.gz'):
        name = name[:-3]
    for ext in FASTQ_EXTS:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


def build_path_map(names, paths):
    """Map collection element names, and their extension-less stems, to dataset paths."""
    exact = dict(zip(names, paths))
    stems = {}
    for name, path in exact.items():
        stems.setdefault(strip_fastq_ext(name), []).append(path)
    return exact, stems


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples-csv',       required=True)
    parser.add_argument('--output',            required=True)
    parser.add_argument('--collection-names',  nargs='*', default=[])
    parser.add_argument('--collection-paths',  nargs='*', default=[])
    parser.add_argument('--reads-dir',         default='reads')
    args = parser.parse_args()

    # ── Stage 1: input file exists ──
    if not os.path.isfile(args.samples_csv):
        err('samples CSV', f'File not found: {args.samples_csv}')

    # ── Stage 2: collection consistency ──
    if len(args.collection_names) != len(args.collection_paths):
        err('collection',
            f'--collection-names ({len(args.collection_names)}) and '
            f'--collection-paths ({len(args.collection_paths)}) have different counts')
    path_map, stem_map = build_path_map(args.collection_names, args.collection_paths)

    # ── Stage 3: parse and validate CSV columns ──
    with open(args.samples_csv, newline='') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []

        col = {}
        missing = []
        for key, aliases in REQUIRED_COLS.items():
            found = resolve_col(fieldnames, aliases)
            if not found:
                missing.append(f'{key} (accepted: {", ".join(aliases)})')
            else:
                col[key] = found
        if missing:
            err('columns', 'Missing required column(s):\n  ' + '\n  '.join(missing)
                + f'\n  Found: {", ".join(fieldnames)}')

        rev_col = None
        for aliases in OPTIONAL_READ_COL.values():
            found = resolve_col(fieldnames, aliases)
            if found:
                rev_col = found
                break

        tax_present = [c for c in TAX_COLS if c in fieldnames]
        rows = list(reader)

    # ── Stage 4: sample ID validation ──
    ids = [r[col['ID']] for r in rows]
    bad_ids = [i for i in ids if '_' in i]
    if bad_ids:
        err('sample IDs',
            'Sample IDs must not contain underscores. '
            f'Offending: {", ".join(bad_ids)}. Use hyphens instead.')
    dupes = [i for i in ids if ids.count(i) > 1]
    if dupes:
        err('sample IDs', f'Duplicate IDs: {", ".join(set(dupes))}')

    # ── Stage 5: resolve filenames to symlinks of collection datasets ──
    os.makedirs(args.reads_dir, exist_ok=True)

    def resolve(filename, label):
        if not filename:
            return ''
        if filename in path_map:
            return link_read(filename, path_map[filename], args.reads_dir)
        # Galaxy collections often drop extensions (test_R1 vs test_R1.fastq.gz)
        candidates = stem_map.get(strip_fastq_ext(filename), [])
        if len(candidates) == 1:
            return link_read(filename, candidates[0], args.reads_dir)
        if len(candidates) > 1:
            err('file resolution',
                f'{label}: "{filename}" matches more than one collection element; '
                'make the collection element names unique.')
        err('file resolution',
            f'{label}: "{filename}" not found in the collection.\n'
            f'  Available: {", ".join(sorted(path_map))}')

    out_rows = []
    for r in rows:
        out = {
            'ID':      r[col['ID']],
            'forward': resolve(r[col['forward']], f"sample {r[col['ID']]} forward"),
        }
        if rev_col:
            out['reverse'] = resolve(r.get(rev_col, ''), f"sample {r[col['ID']]} reverse")
        for tc in tax_present:
            out[tc] = r.get(tc, '')
        out_rows.append(out)

    # ── Stage 6: write output ──
    out_cols = ['ID', 'forward'] + (['reverse'] if rev_col else []) + tax_present
    try:
        with open(args.output, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=out_cols)
            w.writeheader()
            w.writerows(out_rows)
    except OSError as e:
        err('output', f'Could not write {args.output}: {e}')

    print(f'[build_samples_csv] resolved {len(out_rows)} sample(s) to {args.output}')


if __name__ == '__main__':
    main()
