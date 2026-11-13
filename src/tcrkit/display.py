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

__all__ = ['print_residues', 'mark_string', 'SYMBOLS', 'COLORS']

#: One character per value, for each column that can be printed. Anything not listed
#: prints as its first character; empty and ``'-'`` print as a blank.
SYMBOLS = {
    'topology': {'signal_peptide': 'S', 'extracellular': 'O', 'transmembrane': 'M',
                 'cytoplasmic': 'I'},
    'domain': {'G-ALPHA1': '1', 'G-ALPHA2': '2', 'G-ALPHA': 'a', 'G-BETA': 'b',
               'C-LIKE': 'c', 'V': 'V', 'C': 'C'},
    'gene_segment': {'L': 'L', 'V': 'V', 'J': 'J', 'N/D': 'n', 'C': 'C'},
    'region': {'A': 'A', 'AB': '.', 'B': 'B', 'BC': '.', 'C': 'C', 'CD': '.', 'D': 'D',
               'DE': ':', 'E': 'E', 'EF': ':', 'F': 'F', 'FG': ':', 'G': 'G',
               'helix': 'H',
               'FR1': 'f', 'FR2': 'f', 'FR3': 'f', 'FR4': 'f',
               'CDR1': '1', 'CDR2': '2', 'CDR3': '3'},
    'helix': {True: 'H', False: ' '},
}

_MEANING = {
    'topology': {'S': 'signal peptide', 'O': 'extracellular', 'M': 'transmembrane',
                 'I': 'cytoplasmic'},
    'domain': {'1': 'G-ALPHA1', '2': 'G-ALPHA2', 'a': 'G-ALPHA', 'b': 'G-BETA',
               'c': 'C-LIKE', 'V': 'V domain', 'C': 'constant domain'},
    'gene_segment': {'L': 'leader (V gene)', 'V': 'V gene', 'J': 'J gene',
                     'n': 'N/D (neither V nor J)',
                     'C': 'C gene'},
    'region': {'A': 'strand A', 'B': 'strand B', 'C': 'strand C', 'D': 'strand D',
               'E': 'strand E', 'F': 'strand F', 'G': 'strand G',
               '.': 'turn (AB, BC, CD)', ':': 'turn (DE, EF, FG)', 'H': 'helix region (IMGT 50-92)',
               'f': 'framework', '1': 'CDR1', '2': 'CDR2', '3': 'CDR3'},
    'helix': {'H': 'groove-helix residue'},
}

_DEFAULT = ('topology', 'domain', 'gene_segment', 'region', 'helix')

#: Terminal colours (ANSI codes) per printed character, for ``color=``. Characters not
#: listed print uncoloured.
COLORS = {
    'region': {'A': '34', 'B': '34', 'C': '34', 'D': '34', 'E': '34', 'F': '34', 'G': '34',
               '.': '90', ':': '90', 'H': '31',
               '1': '33', '2': '32', '3': '35'},
    'topology': {'S': '33', 'M': '35', 'I': '36'},
    'domain': {'1': '34', '2': '32', 'a': '34', 'b': '32', 'V': '34', 'C': '32', 'c': '32'},
    'gene_segment': {'L': '90', 'V': '34', 'J': '33', 'n': '35', 'C': '32'},
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


def show(residues: pd.DataFrame, option) -> None:
    """What ``show=`` does in :func:`tcrkit.anarci` / :func:`tcrkit.annotate_mhc`: print
    `residues` in colour if `option` is ``True``, or with `option` as print_residues
    options if it is a dict (colour on unless it says ``'color': False``); nothing if it
    is false. Returns None."""
    if option:
        print_residues(residues, **{'color': True,
                                    **(option if isinstance(option, dict) else {})})


def mark_string(values, column: str) -> str:
    """One character per value, as :func:`print_residues` prints `column`, with ``'-'``
    where it prints a blank: ``mark_string(['FR1', 'CDR1'], 'region') == 'f1'``."""
    return ''.join(_symbol(column, v) for v in values).replace(' ', '-')


def _ruler(start: int, length: int) -> str:
    """Position numbers every 10 residues: ``'1         11        21'``."""
    line = [' '] * (length + 8)
    for offset in range(0, length, 10):
        label = str(start + offset)
        line[offset:offset + len(label)] = label
    return ''.join(line).rstrip()


def print_residues(residues: pd.DataFrame, columns=None, width: int | None = 60,
                   file=None, color=False) -> None:
    """Print each sequence of a per-residue table with its annotation underneath.

    Args:
        residues: the per-residue table of :func:`tcrkit.annotate_mhc` or
            :func:`tcrkit.anarci` (``per_residue=True``): one row per residue, with
            ``position`` and ``aa``, indexed by sequence.
        columns: which annotation lines to print, top to bottom; default: whichever of
            ``topology``, ``domain``, ``gene_segment``, ``region`` and ``helix`` the
            table has (``domain`` is left out when ``gene_segment`` is there: on a TCR
            it says the same). Any column works; see :data:`SYMBOLS` for how values
            become characters.
        width: residues per line; ``None`` puts each whole sequence on one line.
        file: where to print (default: stdout).
        color: colour the residues and one annotation line by that line's values, with
            terminal colour codes (shown by terminals and Jupyter). ``True`` colours by
            ``region``; or name the line, e.g. ``'topology'``. See :data:`COLORS`.
    """
    out = file or sys.stdout
    if not {'position', 'aa'} <= set(residues.columns):
        raise ValueError('print_residues needs a per-residue table: '
                         'annotate_mhc(..., per_residue=True) or anarci(..., per_residue=True)')
    if columns is None:
        columns = [c for c in _DEFAULT if c in residues.columns]
        if 'gene_segment' in columns and 'domain' in columns:
            columns.remove('domain')
    columns = [c for c in columns if c in residues.columns]
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
        step = width or max(len(seq), 1)
        for i in range(0, len(seq), step):
            print(' ' * pad + _ruler(int(res['position'].iloc[i]), len(seq[i:i + step])),
                  file=out)
            chunk = seq[i:i + step]
            key = tracks[by][i:i + step] if by else ''
            print(f"{'sequence':<{pad}}{_paint(chunk, key, palette) if by else chunk}", file=out)
            for c, t in tracks.items():
                line = t[i:i + step].rstrip()
                if c == by:
                    line = _paint(line, key, palette)
                print(f'{c:<{pad}}{line}'.rstrip(), file=out)
            print(file=out)
    legend = []
    for c in columns:
        for line in _legend(c, used[c], palette if c == by else {}, pad):
            legend.append(line)
    if legend:
        print('\n'.join(legend), file=out)


def _legend(column: str, chars, palette: dict, pad: int, width: int = 80) -> list[str]:
    """The legend of one line, wrapped at `width`: ``region  A B C D strands, H helix``."""
    meaning = _MEANING.get(column, {})
    order = list(meaning)
    chars = sorted(chars, key=lambda ch: (order.index(ch) if ch in order else len(order), ch))
    strands = [ch for ch in chars if meaning.get(ch, '').startswith('strand ')]
    items = []                                       # (plain text, painted text)
    if strands:
        items.append((' '.join(strands) + (' strands' if len(strands) > 1 else ' strand'),
                      ' '.join(_paint(s, s, palette) for s in strands)
                      + (' strands' if len(strands) > 1 else ' strand')))
    for ch in chars:
        if ch not in strands:
            text = f' {meaning[ch]}' if meaning.get(ch) else ''
            items.append((ch + text, _paint(ch, ch, palette) + text))
    lines, plain, painted = [], '', ''
    for p, q in items:
        if plain and pad + len(plain) + 2 + len(p) > width:
            lines.append(painted + ',')
            plain = painted = ''
        plain, painted = (f'{plain}, {p}', f'{painted}, {q}') if plain else (p, q)
    if plain:
        lines.append(painted)
    return [f"{column if i == 0 else '':<{pad}}{line}" for i, line in enumerate(lines)]
