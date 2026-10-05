"""Command-line interface: one subcommand per capability.

    tcrkit stitch         build a full TCR nt/aa sequence from V/J calls + CDR3
    tcrkit anarci         annotate V/J/CDR regions of a protein sequence
    tcrkit normalize-tcr  standardize TCR gene/allele names to IMGT
    tcrkit normalize-mhc  standardize HLA/MHC allele names
    tcrkit add-cdr3-junctions   put the IMGT junction residues back on a CDR3

These are separate steps that compose: a CDR3 reported as IMGT 105-117 will not stitch,
so repair it with ``add-cdr3-junctions`` first and feed the result to ``stitch``.

Every subcommand works on values given on the command line, and each also has a table
mode (``--in-csv`` plus the relevant ``--*-col`` flags) that reads a CSV/TSV and writes
one back. Output goes to stdout as a readable table unless ``--out-csv`` or ``--json``
is given.
"""

from __future__ import annotations

import argparse
import json
import sys

__all__ = ['main']


# ---------------------------------------------------------------------------
# Shared I/O helpers
# ---------------------------------------------------------------------------

def _read_table(path, sep=None):
    """Read a CSV/TSV, inferring the separator from the extension unless given."""
    import pandas as pd
    if sep is None:
        sep = '\t' if str(path).lower().endswith(('.tsv', '.tab', '.txt')) else ','
    return pd.read_csv(path, sep=sep)


def _emit(df, args, index=True):
    """Write a result frame to --out-csv, or to stdout as JSON or a readable table."""
    import pandas as pd

    if getattr(args, 'out_csv', None):
        sep = '\t' if str(args.out_csv).lower().endswith(('.tsv', '.tab')) else ','
        df.to_csv(args.out_csv, sep=sep, index=index)
        print(f"wrote {len(df):,} rows to {args.out_csv}", file=sys.stderr)
        return
    if getattr(args, 'json', False):
        payload = df.to_dict(orient='records' if not index else 'index')
        print(json.dumps(payload, indent=2, default=str))
        return
    with pd.option_context('display.max_columns', None, 'display.width', 0,
                           'display.max_colwidth', 60):
        print(df.to_string(index=index))


def _emit_scalars(pairs, args):
    """Write ``[(input, result), ...]`` as JSON or as aligned 'input -> result' lines."""
    if getattr(args, 'json', False):
        print(json.dumps([{'input': i, 'result': r} for i, r in pairs], indent=2,
                         default=str))
        return
    width = max((len(str(i)) for i, _ in pairs), default=0)
    for i, r in pairs:
        print(f"{str(i):<{width}}  ->  {'-' if r is None else r}")


def _require(args, *names):
    """Fail cleanly when a subcommand was given neither a value nor a table."""
    if not any(getattr(args, n, None) for n in names):
        flags = ' / '.join('--' + n.replace('_', '-') for n in names)
        raise SystemExit(f"error: give either {flags}. See --help.")


# ---------------------------------------------------------------------------
# stitch
# ---------------------------------------------------------------------------

def _cmd_stitch(args):
    from .stitch import stitch

    if args.in_csv:
        df = _read_table(args.in_csv, args.sep)
        mapping = {'v': args.v_col, 'j': args.j_col, 'cdr3': args.cdr3_col}
        if args.species_col:
            mapping['species'] = args.species_col
        out = stitch(df, col_mapping=mapping, chain=args.chain, species=args.species,
                     no_leader=args.no_leader, constant=args.constant,
                     engine=args.engine, progress_bar=not args.quiet,
                     run_anarci=args.run_anarci)
        if args.join:
            out = df.join(out.drop(columns=['v', 'j', 'cdr3'], errors='ignore'),
                          rsuffix='_stitched')
        return _emit(out, args)

    _require(args, 'v')
    out = stitch(args.v, args.j, args.cdr3, chain=args.chain, species=args.species,
                 no_leader=args.no_leader, constant=args.constant,
                 engine=args.engine, run_anarci=args.run_anarci, progress_bar=False)
    if args.json:
        print(json.dumps(out.iloc[0].to_dict(), indent=2, default=str))
    else:
        for k, v in out.iloc[0].items():
            print(f"{k}: {v}")


def _add_stitch(sub):
    p = sub.add_parser(
        'stitch', help='build a full TCR nt/aa sequence from V/J calls + CDR3',
        description='Rebuild a full-length TCR sequence with Stitchr. The CDR3 must '
                    'already carry its IMGT junction residues - run '
                    "'tcrkit fix-junction' first if it does not.",
        epilog="examples:\n"
               "  tcrkit stitch --v TRBV19 --j TRBJ2-7 --cdr3 CASSIRSSYEQYF\n"
               "  tcrkit stitch --in-csv tcrs.csv --v-col v --j-col j --cdr3-col cdr3 \\\n"
               "                --chain beta --out-csv stitched.csv\n"
               "\n"
               "a CDR3 reported as IMGT 105-117 needs repairing first:\n"
               "  tcrkit add-cdr3-junctions --in-csv tcrs.csv --cdr3-col cdr3 \\\n"
               "                            --j-col j --out-col cdr3 --out-csv fixed.csv\n"
               "  tcrkit stitch --in-csv fixed.csv --cdr3-col cdr3 --v-col v --j-col j",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--v', help='V gene / allele call')
    p.add_argument('--j', help='J gene / allele call')
    p.add_argument('--cdr3', help='CDR3 amino-acid sequence')
    p.add_argument('--chain', default='beta',
                   choices=['alpha', 'beta', 'gamma', 'delta'], help='(default: beta)')
    p.add_argument('--species', default='human',
                   help="'human' or 'mouse'; spellings like 'Mus musculus' are "
                        "understood (default: human)")
    p.add_argument('--constant',
                   help="constant-region allele, or 'auto' to let Stitchr choose "
                        "with its J-gene-aware rule (default: tcrkit's fixed "
                        'DEFAULT_CONSTANT)')
    p.add_argument('--no-leader', action='store_true',
                   help="omit the V gene's leader sequence")
    p.add_argument('--run-anarci', action='store_true',
                   help='re-read the stitched sequence with ANARCI')
    p.add_argument('--engine', default='stitchr', choices=['stitchr', 'thimble'],
                   help="'thimble' runs Stitchr's CLI once for the whole input "
                        '(default: stitchr)')
    p.add_argument('--in-csv', help='CSV/TSV to stitch row-wise instead')
    p.add_argument('--v-col', default='v_call', help='(default: v_call)')
    p.add_argument('--j-col', default='j_call', help='(default: j_call)')
    p.add_argument('--cdr3-col', default='cdr3_aa', help='(default: cdr3_aa)')
    p.add_argument('--species-col', help='read species per row from this column')
    p.add_argument('--join', action='store_true',
                   help='join the results onto the input columns')
    p.set_defaults(func=_cmd_stitch)
    return p


# ---------------------------------------------------------------------------
# anarci
# ---------------------------------------------------------------------------

def _cmd_anarci(args):
    from .anarci import anarci, read_sequences

    if args.in_csv:
        df = _read_table(args.in_csv, args.sep)
        if args.seq_col not in df.columns:
            raise SystemExit(f"error: --seq-col {args.seq_col!r} is not a column in "
                             f"{args.in_csv}. Columns: {list(df.columns)}")
        seqs = df[args.seq_col]
    elif args.file:
        seqs = read_sequences(args.file)
    else:
        _require(args, 'seq')
        seqs = list(args.seq)

    # `allow` must be omitted entirely, not passed as None: its mere presence turns on
    # ANARCI's chain-type restriction.
    extra = {'allow': set(args.allow)} if args.allow else {}
    out = anarci(seqs, progress_bar=not args.quiet,
                 include_imgt_index=args.include_imgt, **extra)
    if args.in_csv and args.join:
        out = df.join(out.drop(columns=['sequence'], errors='ignore'),
                      rsuffix='_anarci')
    _emit(out, args)


def _add_anarci(sub):
    p = sub.add_parser(
        'anarci', help='annotate V/J/CDR regions of a protein sequence',
        description='IMGT-number protein sequences with ANARCI. Uses an in-process '
                    'pyhmmer backend when no hmmscan binary is available.',
        epilog="examples:\n"
               "  tcrkit anarci --seq MNAGVTQTPKFRVLKTGQSMTLLCAQDMNHDY...\n"
               "  tcrkit anarci --file seqs.fasta --out-csv numbered.csv\n"
               "  tcrkit anarci --in-csv tcrs.csv --seq-col sequence_aa --allow B",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--seq', nargs='+', help='one or more amino-acid sequences')
    p.add_argument('--file', help='plain-text (one per line) or FASTA file')
    p.add_argument('--in-csv', help='CSV/TSV with a column of sequences')
    p.add_argument('--seq-col', default='sequence_aa', help='(default: sequence_aa)')
    p.add_argument('--allow', nargs='+', metavar='CHAIN',
                   help='restrict to these ANARCI chain types. Pass "A B" for '
                        'alpha/beta TCRs - without it an alpha chain can be reported '
                        'as delta with a TRDJ gene')
    p.add_argument('--include-imgt', action='store_true',
                   help='add an imgt_dict column of position -> residue')
    p.add_argument('--join', action='store_true',
                   help='join the results onto the input columns (--in-csv only)')
    p.set_defaults(func=_cmd_anarci)
    return p


# ---------------------------------------------------------------------------
# normalize-tcr
# ---------------------------------------------------------------------------

def _cmd_normalize_tcr(args):
    from .normalize import normalize_tcr

    if args.in_csv:
        df = _read_table(args.in_csv, args.sep)
        if not args.cols:
            raise SystemExit('error: --in-csv needs --cols naming the columns to '
                             'standardize.')
        out = df.copy()
        for col in args.cols:
            mapping = {'symbol': col}
            if args.species_col:
                mapping['species'] = args.species_col
            res = normalize_tcr(df, col_mapping=mapping, species=args.species,
                                precision=args.precision)
            if args.keep_original:
                out[f'{col}_original'] = df[col]
            out[f'{col}{args.suffix}' if args.suffix else col] = res['standardized']
        return _emit(out, args, index=False)

    _require(args, 'symbols')
    out = normalize_tcr(list(args.symbols), species=args.species,
                        precision=args.precision)
    _emit_scalars(list(zip(out['symbol'], out['standardized'])), args)


def _add_normalize_tcr(sub):
    p = sub.add_parser(
        'normalize-tcr', help='standardize TCR gene/allele names to IMGT',
        description='Standardize TCR gene / allele symbols to IMGT nomenclature '
                    'with tidytcells.',
        epilog="examples:\n"
               "  tcrkit normalize-tcr TCRBV20S1 TRAV23/6\n"
               "  tcrkit normalize-tcr TRBV20-1*01 --precision gene\n"
               "  tcrkit normalize-tcr --in-csv tcrs.csv --cols v_call j_call \\\n"
               "                       --species human --out-csv clean.csv",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('symbols', nargs='*', help='gene / allele symbols to standardize')
    p.add_argument('--species', default='any',
                   help="'human', 'mouse', or 'any' to try both (default: any)")
    p.add_argument('--precision', default='allele', choices=['allele', 'gene'],
                   help='(default: allele)')
    p.add_argument('--in-csv', help='CSV/TSV to standardize columns of')
    p.add_argument('--cols', nargs='+', help='columns to standardize (--in-csv)')
    p.add_argument('--species-col', help='read species per row from this column')
    p.add_argument('--suffix', help="write to '<col><suffix>' instead of replacing")
    p.add_argument('--keep-original', action='store_true',
                   help="also keep the untouched values in '<col>_original'")
    p.set_defaults(func=_cmd_normalize_tcr)
    return p


# ---------------------------------------------------------------------------
# normalize-mhc
# ---------------------------------------------------------------------------

def _cmd_normalize_mhc(args):
    from .normalize import normalize_mhc

    if args.in_csv:
        df = _read_table(args.in_csv, args.sep)
        if args.allele_col not in df.columns:
            raise SystemExit(f"error: --allele-col {args.allele_col!r} is not a column "
                             f"in {args.in_csv}. Columns: {list(df.columns)}")
        alleles = df[args.allele_col]
    else:
        _require(args, 'alleles')
        alleles = list(args.alleles)

    out = normalize_mhc(alleles, paired=args.paired, chain=args.chain,
                        errors=args.errors, unified=args.unified,
                        drop_mutations=args.drop_mutations,
                        restrict_allele_fields=args.fields)
    _emit(out, args, index=bool(args.in_csv))


def _add_normalize_mhc(sub):
    p = sub.add_parser(
        'normalize-mhc', help='standardize HLA/MHC allele names',
        description='Parse and standardize MHC/HLA allele names with mhcgnomes.',
        epilog="examples:\n"
               "  tcrkit normalize-mhc 'A*0201' HLA-B57:01 'H2-Kb'\n"
               "  tcrkit normalize-mhc 'DRB1*04:01' --paired\n"
               "  tcrkit normalize-mhc --in-csv pmhc.csv --allele-col mhc \\\n"
               "                       --errors coerce --out-csv clean.csv",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('alleles', nargs='*', help="allele names, or 'alpha__beta' pairs")
    p.add_argument('--paired', action='store_true',
                   help='return alpha/beta pairs, inferring the partner chain')
    p.add_argument('--chain', choices=['alpha', 'beta'],
                   help='pick one chain of a native pair')
    p.add_argument('--unified', action='store_true',
                   help='express single alleles in the paired layout')
    p.add_argument('--drop-mutations', action='store_true',
                   help='strip mutations into an mhc_mutations column')
    p.add_argument('--fields', type=int, default=2, metavar='N',
                   help='allele fields to keep, 2 -> A*02:01 (default: 2)')
    p.add_argument('--errors', default='raise', choices=['raise', 'coerce', 'ignore'],
                   help="'coerce' records failures in the error column "
                        "(default: raise)")
    p.add_argument('--in-csv', help='CSV/TSV with a column of alleles')
    p.add_argument('--allele-col', default='mhc', help='(default: mhc)')
    p.set_defaults(func=_cmd_normalize_mhc)
    return p


# ---------------------------------------------------------------------------
# fix-junction
# ---------------------------------------------------------------------------

def _cmd_fix_junction(args):
    from .cdr3_junction import add_cdr3_junctions, load_junction_table

    table = load_junction_table(args.table) if args.table else None

    if args.in_csv:
        df = _read_table(args.in_csv, args.sep)
        mapping = {'cdr3': args.cdr3_col, 'j': args.j_col}
        if args.v_col:
            mapping['v'] = args.v_col
        if args.species_col:
            mapping['species'] = args.species_col
        res = add_cdr3_junctions(df, col_mapping=mapping, table=table,
                             keep_unknown=not args.drop_unknown)
        out = df.copy()
        out[args.out_col or args.cdr3_col] = res['cdr3_fixed']
        return _emit(out, args, index=False)

    _require(args, 'cdr3')
    if not args.j:
        raise SystemExit('error: --j is required (it determines the 118 residue).')
    res = add_cdr3_junctions(list(args.cdr3), args.j, args.v, args.species,
                             table=table)
    _emit_scalars(list(zip(res['cdr3'], res['cdr3_fixed'])), args)


def _add_fix_junction(sub):
    p = sub.add_parser(
        'add-cdr3-junctions', help='put the IMGT junction residues back on a CDR3',
        description='Put the IMGT junction residues (104 C ... 118 F/W/C/L/V) back '
                    'on CDR3s reported as 105-117 only.',
        epilog="examples:\n"
               "  tcrkit add-cdr3-junctions --cdr3 ASSYSGNTEAF --j TRBJ1-1\n"
               "  tcrkit add-cdr3-junctions --cdr3 AVMDSNYQLI --j TRAJ33*01 --species human\n"
               "  tcrkit add-cdr3-junctions --in-csv tcrs.csv --cdr3-col cdr3 \\\n"
               "                            --j-col j_call --out-col cdr3_fixed \\\n"
               "                            --out-csv fixed.csv",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--cdr3', nargs='+', help='CDR3 amino-acid sequence(s)')
    p.add_argument('--j', help='J gene / allele call (required)')
    p.add_argument('--v', help='V gene / allele call (optional, more specific)')
    p.add_argument('--species', help="'human' / 'mouse' (detected from the gene if omitted)")
    p.add_argument('--table', help='a junction lookup CSV to use instead of the bundled one')
    p.add_argument('--in-csv', help='CSV/TSV to repair a column of')
    p.add_argument('--cdr3-col', default='cdr3_aa', help='(default: cdr3_aa)')
    p.add_argument('--j-col', default='j_call', help='(default: j_call)')
    p.add_argument('--v-col', help='V call column (optional)')
    p.add_argument('--species-col', help='species column (optional)')
    p.add_argument('--out-col', help='write to this column (default: overwrite --cdr3-col)')
    p.add_argument('--drop-unknown', action='store_true',
                   help='leave rows whose J gene is unknown empty instead of unchanged')
    p.set_defaults(func=_cmd_fix_junction)
    return p


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def build_parser():
    from . import __version__

    parser = argparse.ArgumentParser(
        prog='tcrkit',
        description='TCR-processing toolkit: stitch, ANARCI numbering, gene/allele '
                    'normalization and CDR3 junction repair.',
        epilog="Run 'tcrkit <command> --help' for a command's options and examples.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version=f'tcrkit {__version__}')

    sub = parser.add_subparsers(dest='command', metavar='<command>')
    for add in (_add_stitch, _add_anarci, _add_normalize_tcr, _add_normalize_mhc,
                _add_fix_junction):
        p = add(sub)
        # Options every subcommand shares, added last so they sort after the specific ones.
        p.add_argument('--out-csv', help='write results to this CSV/TSV')
        p.add_argument('--json', action='store_true', help='print results as JSON')
        p.add_argument('--sep', help='input separator (inferred from the extension)')
        p.add_argument('-q', '--quiet', action='store_true',
                       help='suppress progress bars and summaries')
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, 'command', None):
        parser.print_help()
        return 1
    try:
        args.func(args)
    except BrokenPipeError:          # e.g. piped into `head`
        return 0
    except KeyboardInterrupt:
        return 130
    return 0


if __name__ == '__main__':
    sys.exit(main())
