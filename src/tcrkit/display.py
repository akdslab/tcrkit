"""Print a sequence with its annotation underneath it, one line per annotation::

    print_residues(annotate_mhc(seq, per_residue=True))

    A0201  365 residues
              1         11        21        31        41        51
    sequence  MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRF
    topology  SSSSSSSSSSSSSSSSSSSSSSSSOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    domain                            111111111111111111111111111111111111
    ...

Works on the per-residue tables of :func:`tcrkit.annotate_mhc` and :func:`tcrkit.anarci`.
"""

from __future__ import annotations

import sys

import pandas as pd

__all__ = ['print_residues', 'SYMBOLS', 'COLORS']

#: One character per value, for each column that can be printed. Anything not listed
#: prints as its first character; empty and ``'-'`` print as a blank.
SYMBOLS = {
    'topology': {'signal_peptide': 'S', 'extracellular': 'O', 'transmembrane': 'M',
                 'cytoplasmic': 'I'},
    'domain': {'G-ALPHA1': '1', 'G-ALPHA2': '2', 'G-ALPHA': 'a', 'G-BETA': 'b'},
    'region': {'A': 'A', 'AB': '.', 'B': 'B', 'BC': '.', 'C': 'C', 'CD': '.', 'D': 'D',
               'helix': 'H',
               'FR1': 'f', 'FR2': 'f', 'FR3': 'f', 'FR4': 'f',
               'CDR1': '1', 'CDR2': '2', 'CDR3': '3'},
    'helix': {True: 'H', False: ' '},
}

_MEANING = {
    'topology': {'S': 'signal peptide', 'O': 'extracellular', 'M': 'transmembrane',
                 'I': 'cytoplasmic'},
    'domain': {'1': 'G-ALPHA1', '2': 'G-ALPHA2', 'a': 'G-ALPHA', 'b': 'G-BETA'},
    'region': {'A': 'strand A', 'B': 'strand B', 'C': 'strand C', 'D': 'strand D',
               '.': 'turn (AB, BC, CD)', 'H': 'helix region (IMGT 50-92)',
               'f': 'framework', '1': 'CDR1', '2': 'CDR2', '3': 'CDR3'},
    'helix': {'H': 'groove-helix residue'},
}

_DEFAULT = ('topology', 'domain', 'region', 'helix')

#: Terminal colours (ANSI codes) per printed character, for ``color=``. Characters not
#: listed print uncoloured.
COLORS = {
    'region': {'A': '34', 'B': '34', 'C': '34', 'D': '34', '.': '90', 'H': '31',
               '1': '33', '2': '32', '3': '35'},
    'topology': {'S': '33', 'M': '35', 'I': '36'},
    'domain': {'1': '34', '2': '32', 'a': '34', 'b': '32'},
    'helix': {'H': '31'},
}


def _paint(text: str, key: str, palette: dict) -> str:
    """`text` with each character coloured by the matching character of `key`."""
    out, current = [], None
    for ch, k in zip(text, key):
        code = palette.get(k)
        if code != current:
            out.append('\033[0m' if code is None else f'\033[{code}m')
            current = code
        out.append(ch)
    if current is not None:
        out.append('\033[0m')
    return ''.join(out)


def _symbol(column: str, value) -> str:
    if isinstance(value, float) and pd.isna(value):
        return ' '
    table = SYMBOLS.get(column, {})
    if value in table:
        return table[value]
    text = str(value)
    return ' ' if text in ('', '-') else text[0]


def _ruler(start: int, length: int) -> str:
    """Position numbers every 10 residues: ``'1         11        21'``."""
    line = [' '] * (length + 8)
    for offset in range(0, length, 10):
        label = str(start + offset)
        line[offset:offset + len(label)] = label
    return ''.join(line).rstrip()


def print_residues(residues: pd.DataFrame, columns=None, width: int = 60,
                   file=None, color=False) -> None:
    """Print each sequence of a per-residue table with its annotation underneath.

    Args:
        residues: the per-residue table of :func:`tcrkit.annotate_mhc` or
            :func:`tcrkit.anarci` (``per_residue=True``): one row per residue, with
            ``position`` and ``aa``, indexed by sequence.
        columns: which annotation lines to print, top to bottom; default: whichever of
            ``topology``, ``domain``, ``region`` and ``helix`` the table has. Any column
            works; see :data:`SYMBOLS` for how values become characters.
        width: residues per line.
        file: where to print (default: stdout).
        color: colour the residues and one annotation line by that line's values, with
            terminal colour codes (shown by terminals and Jupyter). ``True`` colours by
            ``region``; or name the line, e.g. ``'topology'``. See :data:`COLORS`.
    """
    out = file or sys.stdout
    if not {'position', 'aa'} <= set(residues.columns):
        raise ValueError('print_residues needs a per-residue table: '
                         'annotate_mhc(..., per_residue=True) or anarci(..., per_residue=True)')
    columns = [c for c in (columns or _DEFAULT) if c in residues.columns]
    by = ('region' if color is True else color) or None
    if by is not None and by not in columns:
        raise ValueError(f'color={color!r}: no {by!r} line to colour by; lines: {columns}')
    palette = COLORS.get(by, {}) if by else {}
    pad = max(len('sequence'), *(len(c) for c in columns)) + 2
    used = {c: set() for c in columns}
    for name in dict.fromkeys(residues.index):
        res = residues.loc[[name]]
        seq = ''.join(res['aa'])
        tracks = {c: ''.join(_symbol(c, v) for v in res[c]) for c in columns}
        for c, t in tracks.items():
            used[c] |= set(t) - {' '}
        print(f'{name}  {len(seq)} residues', file=out)
        for i in range(0, len(seq), width):
            print(' ' * pad + _ruler(int(res['position'].iloc[i]), len(seq[i:i + width])),
                  file=out)
            chunk = seq[i:i + width]
            key = tracks[by][i:i + width] if by else ''
            print(f"{'sequence':<{pad}}{_paint(chunk, key, palette) if by else chunk}", file=out)
            for c, t in tracks.items():
                line = t[i:i + width].rstrip()
                if c == by:
                    line = _paint(line, key, palette)
                print(f'{c:<{pad}}{line}'.rstrip(), file=out)
            print(file=out)
    legend = []
    for c in columns:
        meaning = _MEANING.get(c, {})
        items = [f'{_paint(ch, ch, palette) if c == by else ch} {meaning.get(ch, "")}'.rstrip()
                 for ch in sorted(used[c])]
        if items:
            legend.append(f"{c}: {', '.join(items)}")
    if legend:
        print('\n'.join(legend), file=out)
