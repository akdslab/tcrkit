"""The documentation and skills that come with the package: docs, install_skills, CLI."""

import tomllib
from pathlib import Path

import pytest

import tcrkit
from tcrkit.cli import main
from tcrkit.resources import SKILLS, skills

ROOT = Path(__file__).parents[1]


def test_docs():
    pages = tcrkit.docs()
    assert 'api' in pages
    assert tcrkit.docs('api.md') == tcrkit.docs('api') == pages['api']
    with pytest.raises(FileNotFoundError, match='the pages are .*api'):
        tcrkit.docs('nowhere')


def test_skills_match_the_wheel():
    found = skills()
    assert list(found) == SKILLS
    for name, folder in found.items():
        assert f'name: {name}\n' in (folder / 'SKILL.md').read_text()
    config = tomllib.loads((ROOT / 'pyproject.toml').read_text())
    packaged = config['tool']['hatch']['build']['targets']['wheel']['force-include']
    assert packaged.pop('docs') == 'tcrkit/docs'
    assert packaged == {f'.claude/skills/{s}': f'tcrkit/skills/{s}' for s in SKILLS}


def test_install_skills(tmp_path):
    assert set(tcrkit.install_skills(tmp_path).values()) == {'installed'}
    assert set(tcrkit.install_skills(tmp_path).values()) == {'up to date'}
    edited = tmp_path / SKILLS[0] / 'SKILL.md'
    edited.write_text('my own version')
    assert tcrkit.install_skills(tmp_path)[SKILLS[0]].startswith('kept')
    assert tcrkit.install_skills(tmp_path, force=True)[SKILLS[0]] == 'updated'
    assert edited.read_text() == (skills()[SKILLS[0]] / 'SKILL.md').read_text()


def test_cli(tmp_path, capsys):
    assert main(['docs']) == 0
    assert 'api' in capsys.readouterr().out
    assert main(['docs', 'api', '--path']) == 0
    assert capsys.readouterr().out.strip() == str(tcrkit.docs('api'))
    with pytest.raises(SystemExit, match='no page'):
        main(['docs', 'nowhere'])
    assert main(['skills', 'install', '--target', str(tmp_path)]) == 0
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(SKILLS)
