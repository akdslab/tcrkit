"""The documentation and the Claude Code skills that come with the package::

    tcrkit.docs()                        # {name: path} of every page
    tcrkit.docs('api')                   # the path of docs/api.md
    tcrkit.install_skills()              # copy the skills into ~/.claude/skills
    tcrkit.install_skills(project=True)  # ... or into ./.claude/skills

The pages are the repository's ``docs/`` folder, and the skills those of its
``.claude/skills/`` listed in :data:`SKILLS`. Both are put in the wheel
(``tcrkit/docs``, ``tcrkit/skills``), so they match the installed version; in a source
checkout they are read from the repository itself. On the command line:
``tcrkit docs [name]`` and ``tcrkit skills install``.
"""

from __future__ import annotations

import filecmp
import shutil
from pathlib import Path

__all__ = ['SKILLS', 'docs', 'skills', 'install_skills']

#: The skills for users of the package (``.claude/skills/<name>/SKILL.md``).
SKILLS = ['tcrkit-tcr', 'tcrkit-mhc', 'tcrkit-metrics']

_HERE = Path(__file__).resolve().parent


def _folder(installed: str, checkout: str) -> Path | None:
    """`installed` inside the installed package, else `checkout` in the source checkout."""
    for path in (_HERE / installed, _HERE.parents[1] / checkout):
        if path.is_dir():
            return path
    return None


def docs(name: str | None = None):
    """The documentation pages: ``{name: path}``, or the path of one page.

    Args:
        name: a page, e.g. ``'api'`` (with or without ``.md``).

    Raises:
        FileNotFoundError: no such page (the message lists them), or no pages at all.
    """
    folder = _folder('docs', 'docs')
    pages = {} if folder is None else {p.stem: p for p in sorted(folder.glob('*.md'))
                                         if p.stem != 'README'}
    if not pages:
        raise FileNotFoundError('no documentation found: it comes with the tcrkit wheel '
                                '(tcrkit/docs) or a source checkout (docs/)')
    if name is None:
        return pages
    key = name.removesuffix('.md')
    if key not in pages:
        raise FileNotFoundError(f'no page {name!r}; the pages are {", ".join(pages)}')
    return pages[key]


def skills() -> dict[str, Path]:
    """``{name: folder}`` of the skills in :data:`SKILLS` (each holds a ``SKILL.md``)."""
    root = _folder('skills', '.claude/skills')
    if root is None:
        raise FileNotFoundError('no skills found: they come with the tcrkit wheel '
                                '(tcrkit/skills) or a source checkout (.claude/skills)')
    return {name: root / name for name in SKILLS if (root / name / 'SKILL.md').is_file()}


def install_skills(target: str | Path | None = None, project: bool = False,
                   force: bool = False) -> dict[str, str]:
    """Copy the skills where Claude Code finds them.

    Args:
        target: the skills folder to copy into. Default: ``~/.claude/skills`` (every
            project), or ``./.claude/skills`` with `project`.
        project: install into the current directory's ``.claude/skills``.
        force: replace a skill of the same name that differs (by default it is kept).

    Returns:
        ``{skill: status}``: ``'installed'``, ``'up to date'``, ``'updated'`` or
        ``'kept: differs (force=True to replace)'``.
    """
    if target is None:
        target = (Path.cwd() if project else Path.home()) / '.claude' / 'skills'
    target = Path(target)
    status = {}
    for name, source in skills().items():
        dest = target / name
        if not dest.exists():
            shutil.copytree(source, dest)
            status[name] = 'installed'
        elif _same_tree(source, dest):
            status[name] = 'up to date'
        elif force:
            shutil.rmtree(dest)
            shutil.copytree(source, dest)
            status[name] = 'updated'
        else:
            status[name] = 'kept: differs (force=True to replace)'
    return status


def _same_tree(a: Path, b: Path) -> bool:
    cmp = filecmp.dircmp(a, b)
    if cmp.left_only or cmp.right_only or cmp.funny_files:
        return False
    _, mismatch, errors = filecmp.cmpfiles(a, b, cmp.common_files, shallow=False)
    return not mismatch and not errors and all(
        _same_tree(a / d, b / d) for d in cmp.common_dirs)
