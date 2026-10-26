#!/usr/bin/env python3
"""Offline tests for the ACMG code engine. Stdlib only, no network.

These assert the two behaviours that matter and that are easy to regress:
  1. PM2 does NOT fire when the variant is present in gnomAD (AC > 0).
  2. PM2 does NOT fire when the identifier consistency check failed, even if
     gnomAD reported the variant absent - a typo must never become evidence.

Run: python3 tests/test_engine.py
"""
import json
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
ENGINE = ROOT / "scripts" / "acmg_codes.py"
FIXTURE = ROOT / "examples" / "g245s.json"

failures = []


def run(payload) -> str:
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(payload, f)
        tmp = f.name
    r = subprocess.run([sys.executable, str(ENGINE), tmp, "--scheme", "count3"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        failures.append(f"engine exited {r.returncode}: {r.stderr[:300]}")
    return r.stdout


def check(label, cond, detail=""):
    if cond:
        print(f"  ok   {label}")
    else:
        failures.append(f"{label} {detail}")
        print(f"  FAIL {label} {detail}")


print("case 1: real payload, variant present in gnomAD")
base = json.loads(FIXTURE.read_text())
out = run(base)
check("PP3 applied", "PP3" in out)
check("PM2 withheld (AC=5, not absent)", "PM2 NOT applied" in out,
      f"gnomad ac={base['gnomad']['joint']['ac']}")

print("case 2: gnomAD says absent, but identifiers disagree")
bad = json.loads(json.dumps(base))
bad["gnomad"] = {"variant_id": "17-9999999-C-T", "dataset": "gnomad_r4",
                 "absent": True, "af": None}
bad["consistency"] = {"ok": False, "problems": ["synthetic mismatch for the test"],
                      "identifiers": {}}
out = run(bad)
check("PM2 suppressed on unverified id", "PM2 NOT applied" in out)
check("suppression names the mismatch", "identifiers disagree" in out)

print("case 3: absent AND consistent -> PM2 must fire")
good_absent = json.loads(json.dumps(bad))
good_absent["consistency"] = {"ok": True, "problems": [], "identifiers": {}}
out = run(good_absent)
check("PM2 applied", "PM2" in out and "PM2 NOT applied" not in out)

print()
if failures:
    print(f"FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("all engine tests passed")
