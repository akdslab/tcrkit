---
name: add-capability
description: Rules for changing tcrkit's code - new or changed functions, options, CLI subcommands, docs. Use for any development work on the tcrkit package.
---

# Developing tcrkit

**Hard rules**
- Change nothing outside the tcrkit directory.
- Backward compatible: never change existing functions, arguments, defaults or outputs.
  If a break is needed, explain it and get approval first. New arguments go last.

**Conventions** (every public function)
- One function per capability; takes one value or any collection (`_util.collect`).
- Returns a DataFrame with the input echoed, one `error` column (report, never raise),
  `col_mapping=` for column names (never guess), `squeeze=True` for the bare value.
- Species are `human` / `mouse`. Steps stay separate (e.g. repair CDR3, then stitch).

**Where things go**: `src/` layout (hatchling, uv); export in `src/tcrkit/__init__.py`;
CLI in `cli.py`; run-once data scripts in `src/tcrkit/_build/`; tests in `tests/`.

**Every change**: update `docs/api.md` (all options), `README.md` only for headline
features, a notebook example if useful, and the matching skill in `.claude/skills/`.

```bash
uv sync && uv run pytest        # or .venv/bin/python -m pytest
```
