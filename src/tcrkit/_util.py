"""Internal helpers shared by the capability modules. Not part of the public API.

Three jobs: detecting which shape a caller passed their input in (:func:`collect`),
shaping the result on the way out (:func:`rows`, :func:`squeeze`), and keeping noisy
third-party libraries quiet.
"""

from __future__ import annotations

import contextlib
import io
import logging
import re
import warnings

import pandas as pd

__all__ = [
    'mute_stdout', 'mute_stderr', 'mute_warnings',
    'progress', 'collect', 'rows', 'squeeze', 'per_row', 'drop_unset',
    'species_of', 'SPECIES', 'fmt_warnings', 'is_empty', 'ERROR_COL',
]

#: The one name for "what went wrong" in every frame tcrkit returns.
ERROR_COL = 'error'

#: The species names tcrkit speaks, in and out. Backends are handed whatever spelling
#: they want internally; the caller only ever sees these.
SPECIES = ('human', 'mouse')

# Everything seen in the wild that means one of the above. Matched after stripping
# anything that is not a letter and lowercasing, so 'Mus musculus', 'MUS_MUSCULUS',
# 'Mus Musculus' and 'mus.musculus' all collapse to the same key.
_SPECIES_ALIASES = {
    'human': 'human', 'homosapiens': 'human', 'hsapiens': 'human', 'hs': 'human',
    'hsap': 'human', 'homosapien': 'human',
    'mouse': 'mouse', 'musmusculus': 'mouse', 'mmusculus': 'mouse', 'mm': 'mouse',
    'mmus': 'mouse', 'murine': 'mouse', 'mus': 'mouse',
}


def species_of(value, default='human', allow=()):
    """Any spelling of a species -> ``'human'`` or ``'mouse'``.

    The single place species names are interpreted, so every public function accepts the
    same spellings and reports the same two back. Case, spacing and punctuation are
    ignored, which is what makes a free-text organism column usable::

        species_of('Mus musculus')  -> 'mouse'
        species_of('HOMO SAPIENS')  -> 'human'
        species_of(None)            -> 'human'   (the default)

    Args:
        value: the name to interpret. Empty / NaN / None falls back to `default`.
        default: what an absent value means. Pass ``None`` to get ``None`` back instead.
        allow: extra accepted values passed through unchanged, e.g. ``('any',)`` for
            :func:`tcrkit.normalize_tcr`, which can search both species.

    Raises:
        ValueError: on a name that is neither recognised nor in `allow` - a typo in a
            keyword should not quietly become the default.
    """
    if is_empty(value):
        return default
    raw = str(value).strip()
    if raw in allow or raw.lower() in allow:
        return raw.lower()
    key = re.sub(r'[^a-z]', '', raw.lower())
    if key in _SPECIES_ALIASES:
        return _SPECIES_ALIASES[key]
    raise ValueError(
        f"Unrecognised species {value!r}. Use one of {list(SPECIES) + list(allow)} - "
        f"spellings like 'Mus musculus' or 'homoSapiens' are understood too."
    )


def species_or_default(value, default='human', allow=()):
    """:func:`species_of`, but an unrecognised name falls back instead of raising.

    For species read out of the caller's own data, where free text is expected and one
    odd row should not abort the call. An explicit keyword still uses
    :func:`species_of`, so a typo there is reported.
    """
    try:
        return species_of(value, default=default, allow=allow)
    except ValueError:
        return default


# ---------------------------------------------------------------------------
# Silence helpers
#
# Stitchr prints to stdout, ANARCI prints to stdout/stderr, and pint (pulled in
# by mhcgnomes' dependencies) logs on import. None of that is useful when these
# are called from a library, so each wrapper can mute its backend.
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def mute_stdout(enabled: bool = True):
    if not enabled:
        yield
    else:
        with io.StringIO() as buf, contextlib.redirect_stdout(buf):
            yield


@contextlib.contextmanager
def mute_stderr(enabled: bool = True):
    if not enabled:
        yield
    else:
        with io.StringIO() as buf, contextlib.redirect_stderr(buf):
            yield


@contextlib.contextmanager
def mute_warnings(enabled: bool = True):
    if not enabled:
        yield
    else:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _log = logging.getLogger('pint')
            _prev = _log.getEffectiveLevel()
            _log.setLevel(logging.ERROR)
            try:
                yield
            finally:
                _log.setLevel(_prev)


# ---------------------------------------------------------------------------
# Progress bars
#
# tqdm is a declared dependency, but tcrkit is also usable by dropping src/ on
# sys.path in an environment that may not have it, so degrade to a no-op rather
# than failing at import time.
# ---------------------------------------------------------------------------

def progress(iterable, desc: str = '', total=None, enabled: bool = True):
    """`tqdm` over `iterable` when available and `enabled`, else the iterable itself."""
    if not enabled:
        return iterable
    try:
        from tqdm import tqdm
    except ImportError:
        return iterable
    return tqdm(iterable, desc=desc, total=total)


# ---------------------------------------------------------------------------
# Input normalisation
# ---------------------------------------------------------------------------

def _is_listy(value) -> bool:
    """True for an ordered collection of values, but not for a string."""
    if isinstance(value, (str, bytes, dict)):
        return False
    return isinstance(value, (list, tuple, pd.Series, pd.Index)) or hasattr(value, '__array__')


def collect(fields, args=(), col_mapping=None, optional=(),
            required=None) -> pd.DataFrame:
    """Normalise any supported input shape into one input frame. Internal.

    This is what lets every public function take a single value or a whole table
    through the same argument. `fields` are the required input field names (e.g.
    ``('v', 'j', 'cdr3')``) and `args` the positional values they were given.

    Two calling forms are recognised:

    **Parallel form** - more than the first positional was given, or the first is a
    plain scalar. Each argument becomes a column, scalars broadcasting to the longest
    one::

        collect(('v', 'j'), ('TRBV19', 'TRBJ2-7'))            # 1 row
        collect(('v', 'j'), (['a', 'b'], 'TRBJ2-7'))          # 2 rows, j broadcast

    **Container form** - only the first positional was given and it carries everything:

    =========================  =================================================
    input                      interpretation
    =========================  =================================================
    DataFrame                  columns named after the fields, or via `col_mapping`
    Series                     the values of the first field, index preserved
    ``{'v': [...], ...}``      dict of columns
    ``{'v': .., 'j': ..}``     one record
    ``[{...}, {...}]``         list of records
    ``[(id, value), ...]``     ids and values, for a single-field function
    ``[v1, v2, ...]``          values of the first field
    ``'TRBV19'``               one value
    =========================  =================================================

    `optional` names extra fields that are used when the container happens to carry
    them (e.g. a per-row ``species`` column) and simply left out otherwise, so the
    caller can fall back to its own keyword default.

    `required` narrows which of `fields` a container must actually supply; it defaults
    to all of them. ``add_cdr3_junctions`` takes ``cdr3``/``j``/``v``/``species`` positionally
    but only needs the first two, so it passes ``required=('cdr3', 'j')`` and a frame
    with no V column is fine.

    Returns:
        A DataFrame with one column per field in `fields`, plus any of `optional` that
        were present. The index is carried over from the input where it has one, so
        results can be assigned straight back onto the caller's frame.
    """
    fields = list(fields)
    optional = list(optional)
    required = list(fields if required is None else required)
    args = list(args) + [None] * (len(fields) - len(args))
    first = args[0] if args else None
    rest_given = any(a is not None for a in args[1:])

    def frame(cols: dict, index=None) -> pd.DataFrame:
        out = pd.DataFrame(cols, index=index)
        # keep the declared field order, then whichever optionals turned up
        order = [f for f in fields if f in out.columns]
        order += [f for f in optional if f in out.columns]
        return out[order]

    # ---------------- parallel / single-value form ----------------
    if rest_given or not (isinstance(first, (dict, pd.DataFrame)) or _is_listy(first)):
        n, index = 1, None
        for a in args:
            if _is_listy(a):
                n = max(n, len(a))
                if index is None and isinstance(a, pd.Series):
                    index = a.index
        cols = {}
        for f, a in zip(fields, args):
            if _is_listy(a):
                values = list(a.values) if isinstance(a, pd.Series) else list(a)
                if len(values) not in (n, 1):
                    raise ValueError(
                        f"{f!r} has {len(values)} values but another argument has {n}; "
                        f"parallel arguments must be the same length (or scalar)."
                    )
                cols[f] = values * n if len(values) == 1 else values
            else:
                cols[f] = [a] * n
        return frame(cols, index)

    data = first
    wanted = fields + optional

    # ---------------- DataFrame ----------------
    if isinstance(data, pd.DataFrame):
        m = {**{f: f for f in wanted}, **(col_mapping or {})}
        missing = [m[f] for f in required if m.get(f) not in data.columns]
        if missing:
            raise ValueError(
                f"Missing columns in DataFrame: {missing}. Pass col_mapping to say where "
                f"{required} live, e.g. col_mapping={{{required[0]!r}: 'your_column'}}."
            )
        cols = {f: data[m[f]].values for f in wanted if m.get(f) in data.columns}
        # a non-required field the frame does not carry is simply 'not given'
        cols.update({f: [None] * len(data) for f in fields if f not in cols})
        return frame(cols, data.index)

    # ---------------- Series ----------------
    if isinstance(data, pd.Series):
        return frame({fields[0]: data.values,
                      **{f: [None] * len(data) for f in fields[1:]}}, data.index)

    # ---------------- dict ----------------
    if isinstance(data, dict):
        keyed = {(col_mapping or {}).get(f, f) for f in wanted} & set(data)
        if keyed and all(_is_listy(data[k]) for k in keyed):
            m = {**{f: f for f in wanted}, **(col_mapping or {})}
            n = max(len(data[m[f]]) for f in wanted if m.get(f) in data)
            cols = {}
            for f in wanted:
                key = m.get(f)
                if key in data:
                    values = data[key]
                    cols[f] = list(values.values) if isinstance(values, pd.Series) \
                        else list(values)
                elif f in fields:
                    cols[f] = [None] * n
            return frame(cols)
        if keyed:                                   # one record
            m = {**{f: f for f in wanted}, **(col_mapping or {})}
            return frame({f: [data.get(m.get(f))] for f in wanted if m.get(f) in data
                          or f in fields})
        # id -> value mapping, single-field functions
        return frame({fields[0]: list(data.values()),
                      **{f: [None] * len(data) for f in fields[1:]}},
                     list(data.keys()))

    # ---------------- list / tuple / array ----------------
    items = list(data.values) if isinstance(data, pd.Index) else list(data)
    if not items:
        return frame({f: [] for f in fields})

    if all(isinstance(x, dict) for x in items):                     # records
        m = {**{f: f for f in wanted}, **(col_mapping or {})}
        present = [f for f in wanted
                   if f in fields or any(m.get(f) in r for r in items)]
        return frame({f: [r.get(m.get(f)) for r in items] for f in present})

    tuples = all(isinstance(x, (tuple, list)) for x in items)
    if tuples and len(fields) > 1 and all(len(x) == len(fields) for x in items):
        return frame({f: [x[i] for x in items] for i, f in enumerate(fields)})
    if tuples and len(fields) == 1 and all(len(x) == 2 for x in items):
        return frame({fields[0]: [v for _, v in items]}, [k for k, _ in items])

    return frame({fields[0]: items,                                 # plain values
                  **{f: [None] * len(items) for f in fields[1:]}})


def rows(inp: pd.DataFrame, fields=None):
    """Iterate an input frame as plain dicts, preserving each column's values.

    ``DataFrame.iterrows()`` builds a Series per row and so finds a common dtype for the
    whole row: a column of ``None`` next to a string column comes back as float ``nan``,
    which then reads as the string ``'nan'`` downstream instead of 'not given'. Zipping
    the columns keeps every value exactly as it was stored.
    """
    fields = list(inp.columns) if fields is None else [f for f in fields
                                                      if f in inp.columns]
    for values in zip(*(inp[f] for f in fields)):
        yield dict(zip(fields, values))


def per_row(value, inp: pd.DataFrame, field: str, default, validate=None):
    """Resolve a keyword that may be a scalar, a per-row sequence, or an input column.

    The precedence every public function uses for settings like ``species`` / ``chain``:

    * ``value`` given as a sequence  -> one entry per row
    * ``value`` given as a scalar    -> broadcast to every row (wins over the input)
    * ``value`` left as ``None``     -> the input's own `field` column, if it has one
    * otherwise                      -> `default`

    An explicit keyword beats a column of the same name, so a caller can always override
    what their frame happens to carry.

    `validate` is applied to explicitly-given values only. A typo in a keyword is a bug
    and should fail loudly; the same junk sitting in the caller's own column is free-text
    data and stays the caller's problem to handle leniently.
    """
    if value is None:
        if field in inp.columns:
            return [default if is_empty(v) else v for v in inp[field]]
        return [default] * len(inp)
    if _is_listy(value):
        values = list(value.values) if isinstance(value, pd.Series) else list(value)
        if len(values) != len(inp):
            raise ValueError(
                f"{field}= has {len(values)} entries but the input has {len(inp)} rows."
            )
    else:
        values = [value] * len(inp)
    if validate is not None:
        for v in dict.fromkeys(values):
            validate(v)
    return values


def drop_unset(out: pd.DataFrame, fields) -> pd.DataFrame:
    """Drop echoed input columns the caller never actually supplied.

    A function declares its input fields up front, so :func:`collect` creates a column
    for each - including the optional ones left out. Echoing a column of nothing back is
    just noise, so an all-empty optional field is dropped from the result.
    """
    empty = [f for f in fields if f in out.columns and out[f].isna().all()]
    return out.drop(columns=empty)


def squeeze(out: pd.DataFrame, echo=()):
    """Drop the input echo and collapse a result frame to its plainest useful form.

    What every public function does for ``squeeze=True``. One rule: keep only the result
    columns, then

    * one row, one result column   -> the scalar
    * one row, several columns     -> a dict
    * many rows, one column        -> a Series
    * otherwise                    -> the DataFrame unchanged

    So ``normalize_tcr('TCRBV20S1', squeeze=True)`` gives ``'TRBV20-1'`` and
    ``normalize_tcr([...], squeeze=True)`` gives a Series of standardized symbols.
    """
    res = out.drop(columns=[c for c in echo if c in out.columns])
    if len(res) == 1:
        row = res.iloc[0]
        return row.iloc[0] if res.shape[1] == 1 else row.to_dict()
    if res.shape[1] == 1:
        return res.iloc[:, 0]
    return res


# Values that mean 'not reported' in the sources these wrappers read.
_EMPTY_TOKENS = ('', '-', 'nan', 'none', 'na', 'n/a', '<na>', 'null')


def is_empty(value) -> bool:
    """True for None / NaN / '' / '-' and the other placeholders sources use."""
    try:
        if value is None or (not isinstance(value, (list, tuple)) and pd.isna(value)):
            return True
    except (TypeError, ValueError):
        pass
    return str(value).strip().lower() in _EMPTY_TOKENS


def fmt_warnings(caught) -> str:
    """Recorded warnings as one ' | '-joined string, '-' if there were none.

    DeprecationWarnings are dropped - they are about the libraries, not about
    the sequence being processed.
    """
    msgs = [str(w.message).strip() for w in caught
            if not issubclass(w.category, DeprecationWarning)]
    return ' | '.join(dict.fromkeys(m for m in msgs if m)) or '-'
