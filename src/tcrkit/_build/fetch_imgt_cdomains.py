"""Rebuild ``src/tcrkit/data/imgt_cdomains.csv`` from IMGT's C-DOMAIN protein display.

    python -m tcrkit._build.fetch_imgt_cdomains            # downloads the page
    python -m tcrkit._build.fetch_imgt_cdomains page.html  # or parses a saved copy

The page aligns C domains (IG / TR constant domains) and C-LIKE domains (B2M, MHC
alpha3 / alpha2 / beta2, ...) under one ruler of IMGT positions. ``|`` and ``.`` are
positions; the number printed above a ``|`` sets the count (the C-DOMAIN scheme skips
46-76), otherwise each position is one more than the last. Digit runs are insertions:
counting up after position ``n`` they are ``n.1``, ``n.2`` ...; counting down they
belong to the next position, ``m.7`` ... ``m.1`` (IMGT's ``84.1-84.7, 85.7-85.1``, and
``1.8`` ... ``1.1`` before position 1).

Only human, mouse and rat entries are kept. Data: IMGT, CC BY 4.0 (see
src/tcrkit/data/README.md). The CSV ships with tcrkit; run this only to refresh it.
"""

from __future__ import annotations

import csv
import html
import re
import sys
import urllib.request
from pathlib import Path

__all__ = ['parse', 'main']

URL = 'https://www.imgt.org/IMGTrepertoireMH/Proteins/protein/C-DOMAIN/CLdomains.html'
OUT = Path(__file__).resolve().parent.parent / 'data' / 'imgt_cdomains.csv'
SPECIES = ('Homsap', 'Musmus', 'Ratnor')


def _ruler_labels(tick: str, ruler: str) -> list[str | None]:
    """The IMGT label of every ruler column that is a position, None elsewhere."""
    labels: list[str | None] = [None] * len(ruler)
    tokens = [(m.start(), m.end(), m.group()) for m in re.finditer(r'\d+', tick)]
    num, pending = 0, []                 # pending: columns of a run counting down
    i = 0
    while i < len(ruler):
        c = ruler[i]
        if c in '|.':
            num += 1
            if c == '|':
                # the tick printed over this bar; neighbours can run together ('2326'
                # is 23 and 26), so a tick containing the count confirms it, and a
                # larger one is a jump in the scheme (45 -> 77)
                over = [t for s, e, t in tokens if s <= i < e]
                if over and str(num) not in over[0]:
                    if int(over[0]) < num:
                        raise ValueError(f'ruler column {i}: counted {num}, page prints {over[0]}')
                    num = int(over[0])
            labels[i] = str(num)
            for col, d in pending:
                labels[col] = f'{num}.{d}'
            pending = []
        elif c.isdigit():
            j = i
            while j < len(ruler) and ruler[j].isdigit():
                j += 1
            run = ruler[i:j]
            up = 0
            while up < len(run) and int(run[up]) == up + 1:
                up += 1
            for k in range(up):
                labels[i + k] = f'{num}.{run[k]}'
            pending = [(i + k, run[k]) for k in range(up, len(run))]
            i = j
            continue
        i += 1
    return labels


def parse(text: str) -> list[dict]:
    """One row per entry: species, gene, accession, chain domain, sequence, IMGT labels."""
    lines = html.unescape(re.sub(r'<[^>]+>', '', text)).split('\n')
    head = next(i for i, line in enumerate(lines) if re.match(r'\s*\d+\s+10\s+15', line))
    tick, ruler = lines[head], lines[head + 1]
    labels = _ruler_labels(tick, ruler)
    first = next(i for i, lab in enumerate(labels) if lab is not None)
    rows = []
    for line in lines[head + 2:]:
        m = re.match(r'(\w+)\s*,(\S+)', line)
        if m is None:
            if line.strip().startswith('Dots indicate'):
                break
            continue
        if m.group(1) not in SPECIES:
            continue
        fields = line[m.end():first].split('(')[0].split()   # '(R)' starts the residues
        seq, labs = [], []
        for col in range(first, len(line)):
            c = line[col]
            if not c.isalpha():
                continue
            if col >= len(labels) or labels[col] is None:
                print(f'skipped {m.group(2)}: residue {c!r} at column {col} is off the ruler',
                      file=sys.stderr)
                seq = None
                break
            seq.append(c.upper())
            labs.append(labels[col])
        if not seq:
            continue
        gene = re.sub(r'\(\d+\)$', '', m.group(2))
        number = fields.pop() if fields and re.fullmatch(r'\d|C[HL]\d', fields[-1]) else ''
        rows.append({'species': m.group(1), 'gene': gene,
                     'accession': fields[-1] if fields else '', 'domain': number,
                     'sequence': ''.join(seq), 'imgt': ' '.join(labs)})
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
