"""Rebuild ``src/tcrkit/data/imgt_gdomains.csv`` from IMGT's G-DOMAIN protein display.

    python -m tcrkit._build.fetch_imgt_gdomains            # downloads the page
    python -m tcrkit._build.fetch_imgt_gdomains page.html  # or parses a saved copy

The page aligns MH G domains under a ruler of IMGT positions: ``|`` and ``.``
are positions (one more each), letters are insertions (``7A``, ``54B``), and the
digits after 49 are positions 49.1-49.7. Positions before 1, printed as digits
10..1, become ``-10`` .. ``-1``. Position numbers are counted from 1, and every
number printed above a ``|`` is checked against the count.

Only the MH1 / MH2 sections are kept (not the G-LIKE domains of other proteins),
and only human, mouse and rat entries. Rows whose residues fall off the ruler
are skipped and reported. Data: IMGT, CC BY 4.0 (see src/tcrkit/data/README.md).

The CSV is shipped with tcrkit, so users never need to run this; run it only to
refresh the data from IMGT.
"""

from __future__ import annotations

import csv
import html
import re
import sys
import urllib.request
from pathlib import Path

__all__ = ['parse', 'main']

URL = 'https://www.imgt.org/IMGTrepertoireMH/Proteins/protein/G-DOMAIN/Gdomains.html'
OUT = Path(__file__).resolve().parent.parent / 'data' / 'imgt_gdomains.csv'
SECTIONS = {'G-ALPHA1 [D1] (MH1)': 'G-ALPHA1', 'G-ALPHA2 [D2] (MH1)': 'G-ALPHA2',
            'G-ALPHA [D1] (MH2)': 'G-ALPHA', 'G-BETA [D1] (MH2)': 'G-BETA'}
SPECIES = ('Homsap', 'Musmus', 'Ratnor')


def _ruler_labels(tick: str, ruler: str) -> list[str | None]:
    """The IMGT label of every ruler column that is a position, None elsewhere."""
    labels: list[str | None] = [None] * len(ruler)
    first = ruler.index('|')
    lead = [i for i in range(first) if ruler[i].isdigit()]
    for k, i in enumerate(lead):
        labels[i] = str(-(len(lead) - k))
    num = 0
    for i in range(first, len(ruler)):
        c = ruler[i]
        if c in '|.':
            num += 1
            labels[i] = str(num)
            if c == '|':
                printed = {m.group() for m in re.finditer(r'\d+', tick) if m.start() <= i < m.end()}
                if printed != {str(num)}:
                    raise ValueError(f'ruler column {i}: counted {num}, page prints {printed}')
        elif c.isalpha():
            labels[i] = f'{num}{c}'
        elif c.isdigit():
            labels[i] = f'{num}.{c}'
    return labels


def parse(text: str) -> list[dict]:
    """One row per (domain, entry): sequence and space-separated IMGT labels."""
    lines = html.unescape(re.sub(r'<[^>]+>', '', text)).split('\n')
    rows, seen = [], set()
    for i, line in enumerate(lines):
        domain = SECTIONS.get(line.strip())
        if domain is None:
            continue
        labels = _ruler_labels(lines[i + 4], lines[i + 5])
        j = i + 6
        while lines[j].strip():
            entry = lines[j]
            j += 1
            m = re.match(r'(\S+) (\S+)\s+(?:\(([^)]*)\)|(\S+))', entry)
            if m is None or m.group(1) not in SPECIES:
                continue
            seq, labs = [], []
            for col in range(m.end(), len(entry)):
                c = entry[col]
                if c in ' .()':
                    continue
                if col >= len(labels) or labels[col] is None:
                    print(f'skipped {domain} {m.group(2)}: residue {c!r} at column {col} '
                          'is off the ruler', file=sys.stderr)
                    seq = None
                    break
                seq.append(c.upper())
                labs.append(labels[col])
            if seq is None:
                continue
            row = {'domain': domain, 'species': m.group(1), 'gene': m.group(2),
                   'allele': m.group(3) or m.group(4), 'sequence': ''.join(seq),
                   'imgt': ' '.join(labs)}
            key = (domain, row['gene'], row['allele'], row['sequence'])
            if key not in seen:
                seen.add(key)
                rows.append(row)
    return rows


def main(argv: list[str]) -> None:
    if argv:
        text = Path(argv[0]).read_text()
    else:
        with urllib.request.urlopen(URL, timeout=60) as r:
            text = r.read().decode()
    rows = parse(text)
    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f'{len(rows)} rows -> {OUT}')


if __name__ == '__main__':
    main(sys.argv[1:])
