---
name: release
description: Prepare a tcrkit version bump, changelog entry and commit, and check the repo is ready to push.
---

# Release / commit checklist

1. Version lives only in `src/tcrkit/__init__.py` (`__version__`; `pyproject.toml` reads
   it). Current series `0.1`; every change bumps the patch: `0.1.4` -> `0.1.5`.
2. Every bump gets an entry at the top of `CHANGELOG.md`: one short line per change,
   breaking ones in **bold**. Keep the file uncrowded.
3. `uv run pytest`; re-run both notebooks
   (`jupyter nbconvert --to notebook --execute --inplace notebooks/<nb>.ipynb`).
4. Grep `README.md`, `docs/api.md`, `src/tcrkit/data/README.md` and docstrings for stale
   names, files and counts. Check docs and skills match the change.
5. Any behaviour change to existing code needs the user's approval.
6. Commit only what the user asked; never `sandbox/` or scratch notebooks, and
   `CLAUDE.md`, `CHANGELOG.md`, `.claude/` only when told. One commit unless asked.
