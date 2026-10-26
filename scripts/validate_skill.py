#!/usr/bin/env python3
"""Validate this skill against the Agent Skills spec and report host compatibility.

Stdlib only, no PyYAML: a skill that needs `pip install` before it can be checked
is not portable, and portability is the whole claim being tested here.

Checks:
  1. YAML frontmatter parses and carries the required `name` and `description`.
  2. `name` obeys the spec grammar and matches the parent directory.
  3. `description` is 1-1024 chars.
  4. Frontmatter uses ONLY the six spec fields (warn on host-specific extras, which
     some hosts hard-error on: "Unexpected key(s) in SKILL.md frontmatter").
  5. Every script imports from the standard library only.
  6. SKILL.md stays under the recommended 500 lines.

Usage:
  python3 scripts/validate_skill.py            # validate + compatibility table
  python3 scripts/validate_skill.py --quiet    # exit code only
Exit code 0 = valid, 1 = at least one error.
"""
import argparse
import ast
import pathlib
import re
import sys

SPEC_FIELDS = {"name", "description", "license", "compatibility", "metadata", "allowed-tools"}

# Fields other hosts add. Legal in Claude Code, a hard error when packaging for
# claude.ai / the Skills API, and meaningless elsewhere.
HOST_ONLY_FIELDS = {
    "when_to_use", "argument-hint", "arguments", "disable-model-invocation",
    "user-invocable", "disallowed-tools", "model", "effort", "context", "agent",
    "background", "hooks", "paths", "shell",
}

NAME_RE = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
HOSTS = [
    # host, where it looks, how it resolves the skill
    ("Hermes", "~/.hermes/profiles/<profile>/skills/<category>/<name>/", "auto-registers /<name>"),
    ("Claude Code", "~/.claude/skills/<name>/  or  .claude/skills/<name>/", "/<name>"),
    ("Codex CLI", "~/.agents/skills/<name>/  or  <repo>/.agents/skills/<name>/", "/skills or $<name>"),
    ("Antigravity", "agy plugin import claude  (skill inside a plugin bundle)", "/skills in-session"),
]


def parse_frontmatter(text):
    """Minimal YAML subset: flat `key: value` plus one level of nested map."""
    m = re.match(r"^---\n(.*?)\n---\n", text, re.S)
    if not m:
        return None, "no YAML frontmatter block at the top of SKILL.md"
    data, parent = {}, None
    for raw in m.group(1).splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        if raw.startswith(("  ", "\t")):          # nested key
            if parent is None:
                return None, f"nested value without a parent key: {raw!r}"
            k, _, v = raw.strip().partition(":")
            data[parent][k.strip()] = v.strip().strip('"').strip("'")
            continue
        k, _, v = raw.partition(":")
        k, v = k.strip(), v.strip()
        if not v:                                  # opens a nested map
            parent = k
            data[k] = {}
        else:
            parent = None
            data[k] = v.strip('"').strip("'")
    return data, None


def stdlib_imports(paths):
    """Return {script: [non-stdlib module names]}."""
    std = set(sys.stdlib_module_names)
    bad = {}
    for p in paths:
        mods = set()
        for n in ast.walk(ast.parse(p.read_text(encoding="utf-8"))):
            if isinstance(n, ast.Import):
                mods |= {a.name.split(".")[0] for a in n.names}
            elif isinstance(n, ast.ImportFrom) and n.level == 0 and n.module:
                mods.add(n.module.split(".")[0])
        extra = sorted(mods - std)
        if extra:
            bad[p.name] = extra
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    root = pathlib.Path(__file__).resolve().parent.parent
    errors, warnings, notes = [], [], []

    skill = root / "SKILL.md"
    if not skill.exists():
        print(f"FAIL  SKILL.md missing in {root}", file=sys.stderr)
        return 1
    text = skill.read_text(encoding="utf-8")

    fm, err = parse_frontmatter(text)
    if err:
        errors.append(err)
        fm = {}

    # 1-3. required fields and their constraints
    name, desc = fm.get("name", ""), fm.get("description", "")
    if not name:
        errors.append("required field `name` is missing")
    else:
        if not NAME_RE.fullmatch(name):
            errors.append(f"name {name!r} must be lowercase alphanumerics and single hyphens")
        if len(name) > 64:
            errors.append(f"name is {len(name)} chars, spec max is 64")
        if name != root.name:
            errors.append(f"name {name!r} must match the directory name {root.name!r}")
    if not desc:
        errors.append("required field `description` is missing")
    elif not 1 <= len(desc) <= 1024:
        errors.append(f"description is {len(desc)} chars, spec allows 1-1024")

    # 4. field surface
    for k in fm:
        if k in HOST_ONLY_FIELDS:
            errors.append(
                f"frontmatter field {k!r} is a host-specific extension, not in the spec; "
                "portable skills use only the six spec fields")
        elif k not in SPEC_FIELDS:
            warnings.append(f"unknown frontmatter field {k!r}")
    if fm.get("license") in ("", None):
        notes.append("no `license` field - see LICENSE for the actual grant")
    if not fm.get("compatibility"):
        notes.append("no `compatibility` field: hosts cannot warn about network/python needs")

    # 5. dependency-free scripts
    py = sorted((root / "scripts").glob("*.py")) if (root / "scripts").is_dir() else []
    if not py:
        notes.append("no scripts/ directory: instruction-only skill")
    for script, mods in stdlib_imports(py).items():
        errors.append(f"{script} imports non-stdlib modules {mods}; the skill must run "
                      "without a pip install")

    # 6. size budget
    lines = text.count("\n") + 1
    if lines > 500:
        warnings.append(f"SKILL.md is {lines} lines; spec recommends under 500")

    if not a.quiet:
        n_scripts = len(py)
        print(f"skill      : {name or '(unnamed)'}")
        print(f"version    : {fm.get('metadata', {}).get('version', '(unset)')}")
        print(f"dir        : {root}")
        print(f"description: {len(desc)} chars   body: {lines} lines   scripts: {n_scripts} (stdlib only)")
        print(f"fields     : {', '.join(fm) or '(none)'}")
        print()
        print("Host compatibility")
        print("-" * 78)
        for host, where, how in HOSTS:
            verdict = "OK " if not any("host-specific extension" in e for e in errors) else "HOST-ONLY"
            print(f"{verdict:10} {host:13} {how:26} {where}")
        print()
        for n in notes:
            print(f"note   {n}")
        for w in warnings:
            print(f"warn   {w}")
        for e in errors:
            print(f"FAIL   {e}")
        print()
        print("VALID" if not errors else f"INVALID ({len(errors)} error(s))")
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
