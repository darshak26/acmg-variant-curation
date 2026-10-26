#!/usr/bin/env python3
"""Offline tests for the ACMG code engine and the identifier guard. Stdlib only, no network.

Asserts the behaviours that are easy to regress and expensive to get wrong:
  1. PM2 does NOT fire when the variant is present in gnomAD (AC > 0).
  2. PM2 does NOT fire when the identifiers were not cross-checked.
  3. The lean uses ClinGen SVI points, not a code count.
  4. A genomic (g.) ClinVar record is cross-checked, not silently skipped.
  5. A guard that could not compare anything reports NOT ok.
  6. CNVs and structural variants are rejected, not half-curated.

The ceiling is tested deliberately: the seven public-data codes max out at 6
pathogenic points (PM2 1 + PP3 1 + PS3 4), below the 10-point threshold. This
tool therefore CANNOT reach a pathogenic classification on its own - by design,
and the tests pin that so nobody "fixes" it into over-calling.

Run: python3 tests/test_engine.py
"""
import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
ENGINE = ROOT / "scripts" / "acmg_codes.py"
FETCH = ROOT / "scripts" / "fetch_evidence.py"
FIXTURE = ROOT / "examples" / "g245s.json"

failures = []


def _load_fetch():
    spec = importlib.util.spec_from_file_location("fetch_evidence", FETCH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # module level is definitions only, no network
    return mod


def run(payload, *extra) -> str:
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(payload, f)
        tmp = f.name
    r = subprocess.run([sys.executable, str(ENGINE), tmp, *extra],
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


def payload(**over):
    """Synthetic payload: absent from gnomAD, no predictors, identifiers verified."""
    d = {"fetched_at": "synthetic",
         "clinvar": {"title": "NM_000546.6(TP53):c.733G>A (p.Gly245Ser)",
                     "classification": "Pathogenic", "review_status": "x"},
         "gnomad": {"variant_id": "17-7674230-C-T", "dataset": "gnomad_r4",
                    "absent": True, "af": None},
         "vep": {"scores": {}},
         "consistency": {"ok": True, "problems": [], "checks": ["synthetic"],
                         "identifiers": {}}}
    d.update(over)
    return d


print("case 1: real payload - variant present in gnomAD")
base = json.loads(FIXTURE.read_text())
out = run(base)
check("PP3 applied", "PP3" in out)
check("PM2 withheld (AC>0)", "PM2 NOT applied" in out,
      f"gnomad ac={base['gnomad']['joint']['ac']}")

print("case 2: gnomAD absent, identifiers disagree -> PM2 suppressed")
bad = json.loads(json.dumps(base))
bad["gnomad"] = {"variant_id": "17-9999999-C-T", "dataset": "gnomad_r4",
                 "absent": True, "af": None}
bad["consistency"] = {"ok": False, "problems": ["synthetic mismatch"], "identifiers": {}}
out = run(bad)
check("PM2 suppressed on unverified id", "PM2 NOT applied" in out)
check("suppression names the mismatch", "identifiers disagree" in out)

print("case 3: absent AND verified -> PM2 fires")
out = run(payload())
check("PM2 applied", "PM2" in out and "PM2 NOT applied" not in out)

print("case 4: lean is SVI points, not a vote count")
out = run(payload())
check("one Supporting code is VUS, not pathogenic-leaning",
      "Verdict: VUS" in out and "Pathogenic-leaning" not in out)
check("points line printed", "Points (ClinGen SVI): pathogenic 1 [PM2 1]" in out)
check("1 point explained as insufficient", "below the 10-point threshold" in out)

ceiling = payload(vep={"scores": {"revel": 0.95}})
out = run(ceiling, "--functional", "damaging", "--scheme", "revel")
check("PM2+PP3+PS3 = 6 points still VUS (public-data ceiling)",
      "Verdict: VUS" in out and "pathogenic 6" in out, out.split("Points")[-1][:90])

freq = payload(gnomad={"variant_id": "x", "dataset": "gnomad_r4", "absent": False,
                       "af": 0.2, "faf95_popmax": 0.2})
out = run(freq)
check("high frequency -> BA1 -> BENIGN", "Verdict: BENIGN" in out and "BA1" in out)

print("case 5: ClinVar conflict is flagged, not hidden")
out = run(payload(clinvar={"title": "x", "classification": "Benign",
                           "review_status": "criteria provided"}))
check("conflict note printed", "CONFLICT WITH CLINVAR" in out)
check("verdict stays VUS", "Verdict: VUS" in out)

print("case 6: identifier guard (offline, no network)")
fe = _load_fetch()
cv_g = {"title": "NC_000017.11:g.7579817G>C", "classification": "Benign"}
vep_g = {"start": 7579817, "seq_region_name": "17", "allele_string": "G/C", "strand": 1}
r = fe.check_identity(cv_g, {"variant_id": "17-7579817-G-C"}, vep_g,
                      "NC_000017.11:g.7579817G>C")
check("genomic HGVS is actually cross-checked", r["ok"] and len(r["checks"]) >= 2,
      str(r["checks"]))
r = fe.check_identity(cv_g, {"variant_id": "17-9999999-G-C"}, vep_g,
                      "NC_000017.11:g.7579817G>C")
check("genomic record at the wrong locus fails", not r["ok"])
r = fe.check_identity(None, {"variant_id": "17-7579817-G-C"}, None, None)
check("nothing to compare -> not ok",
      (not r["ok"]) and "no cross-check" in r["problems"][0])
r = fe.check_identity({"title": "NM_000546.6(TP53):c.733G>A (p.Gly245Ser)"},
                      {"variant_id": "17-7674230-C-T"},
                      # live VEP returns G/A (transcript orientation) for this coding
                      # HGVS on minus-strand TP53; complement -> C/T matches gnomAD
                      {"start": 7674230, "seq_region_name": "17",
                       "allele_string": "G/A", "strand": -1},
                      "NM_000546.6:c.733G>A")
check("accession prefix does not cause a false mismatch", r["ok"], str(r["problems"]))
r = fe.check_identity({"title": "NM_000546.6(TP53):c.710T>A (p.Met237Lys)"},
                      None, None, "NM_000546.6:c.733G>A")
check("wrong coding variant still fails", not r["ok"])

print("case 7: out-of-scope variant type is rejected (canned ClinVar response)")
fe._get = lambda url, **kw: json.dumps({"result": {"150702": {
    "title": "GRCh38/hg38 17p13.1(chr17:7544902-7670581)x3",
    "obj_type": "copy number gain",
    "germline_classification": {"description": "Pathogenic"}}}})
try:
    fe.fetch_clinvar("150702")
    check("CNV raises OutOfScope", False, "no exception raised")
except fe.OutOfScope as e:
    check("CNV raises OutOfScope", "out of scope" in str(e).lower())
    check("rejection cites Riggs 2020", "Riggs" in str(e))

print()
if failures:
    print(f"FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("all engine tests passed")
