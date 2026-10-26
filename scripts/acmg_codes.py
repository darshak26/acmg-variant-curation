#!/usr/bin/env python3
"""Apply the 7-code ACMG/AMP 2015 public-database subset to one fetched variant.

Input is the JSON from fetch_evidence.py. Output is an auditable table: every code
carries the number that triggered it, so a reader can recheck the call by hand.

Two in-silico schemes, both defensible, and they can disagree:
  count3  - PP3 when >=3 of the AVAILABLE predictors call damaging (project default).
            NA is excluded from numerator AND denominator; 'indeterminate' scores
            vote neither way.
  revel   - PP3 at REVEL >= 0.75, BP4 at REVEL <= 0.15 (single calibrated predictor).
Declare which scheme was used in Methods. Do not mix them.

PS3/BS3 are never inferred by this script. Pass --functional from a paper you read.

Usage:
  python3 acmg_codes.py D393V.json --scheme count3
  python3 acmg_codes.py D393V.json --special-bs1 0.001 --functional none
"""
import argparse, json, sys

DAMAGING = "pathogenic"
BENIGN = "benign"


def _polyphen_dir(label):
    if not label:
        return None
    l = label.lower()
    if any(k in l for k in ("damaging", "deleterious")):
        return "damaging"
    if "benign" in l or "tolerated" in l:
        return "tolerated"
    return "indeterminate"


def _sift_dir(label):
    if not label:
        return None
    l = label.lower()
    if "deleter" in l:
        return "damaging"
    if "tolerat" in l:
        return "tolerated"
    return "indeterminate"


def classify_predictors(scores):
    """Return per-tool direction and the three buckets. Missing tool -> None, not 'tolerated'."""
    rows, damaging, tolerated, indeterminate = [], [], [], []

    def add(tool, value, direction):
        rows.append({"tool": tool, "value": value, "call": direction})
        if direction == "damaging":
            damaging.append(tool)
        elif direction == "tolerated":
            tolerated.append(tool)
        elif direction == "indeterminate":
            indeterminate.append(tool)

    v = scores.get("revel")
    if v is not None:
        add("REVEL", v, "damaging" if v >= 0.5 else "tolerated")
    v = scores.get("cadd_phred")
    if v is not None:
        add("CADD", v, "damaging" if v >= 20 else "tolerated")
    v = scores.get("bayesdel_addaf")
    if v is not None:
        add("BayesDel", v, "damaging" if v >= 0.13 else ("tolerated" if v <= -0.36 else "indeterminate"))
    d = _polyphen_dir(scores.get("polyphen"))
    if d:
        add("PolyPhen", scores.get("polyphen"), d)
    d = _sift_dir(scores.get("sift"))
    if d:
        add("SIFT", scores.get("sift"), d)
    return rows, damaging, tolerated, indeterminate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jsonfile")
    ap.add_argument("--scheme", choices=["count3", "revel"], default="count3")
    ap.add_argument("--bs1", type=float, default=0.01,
                    help="BS1 frequency threshold; override with gene-specific prevalence")
    ap.add_argument("--ba1", type=float, default=0.05)
    ap.add_argument("--functional", choices=["none", "damaging", "benign"], default="none",
                    help="from a functional assay you read, never inferred")
    ap.add_argument("--pm2-max-af", type=float, default=0.0,
                    help="PM2 'absent' threshold. 0.0 = strictly absent. Raising it is a "
                         "methodological choice and must be declared.")
    a = ap.parse_args()
    d = json.load(open(a.jsonfile, encoding="utf-8"))

    g = d.get("gnomad") or {}
    vep = (d.get("vep") or {}).get("scores") or {}
    af = g.get("af")
    faf95 = g.get("faf95_popmax")
    freq_for_ba1 = faf95 if faf95 is not None else af

    codes, notes = [], []

    # --- frequency codes
    if freq_for_ba1 is not None and freq_for_ba1 > a.ba1:
        codes.append({"code": "BA1", "direction": BENIGN, "strength": "standalone",
                      "basis": f"grpmax AF {freq_for_ba1:.3g} > {a.ba1}",
                      "dataset": g.get("dataset")})
    if af is not None and af > a.bs1:
        codes.append({"code": "BS1", "direction": BENIGN, "strength": "strong",
                      "basis": f"AF {af:.3g} > {a.bs1}", "dataset": g.get("dataset")})
    # gnomAD answers "Variant not found" identically for a typo and for a genuinely
    # ultra-rare allele. PM2 is a pathogenic code, so it may not rest on an
    # identifier nobody verified: require consistency.ok when the payload carries it.
    cons = d.get("consistency") or {}
    identity_ok = cons.get("ok", True)
    if (g.get("absent") or (af is not None and af <= a.pm2_max_af)) and not identity_ok:
        notes.append(
            "PM2 NOT applied: the gnomAD/ClinVar/VEP identifiers disagree "
            f"({'; '.join(cons.get('problems') or [])}). An unverified id returning "
            "'Variant not found' is indistinguishable from a typo, and a typo must "
            "never become pathogenic evidence. Re-run fetch with the correct id.")
    elif g.get("absent") or (af is not None and af <= a.pm2_max_af):
        basis = ("not found in " + str(g.get("dataset"))) if g.get("absent") else \
                f"AF {af:.3g} <= {a.pm2_max_af}"
        codes.append({"code": "PM2", "direction": DAMAGING, "strength": "supporting",
                      "basis": basis, "dataset": g.get("dataset")})
    if af is not None and af > a.pm2_max_af and not g.get("absent"):
        notes.append(
            f"PM2 NOT applied: variant is present in {g.get('dataset')} "
            f"(ac={g.get('joint', {}).get('ac') if isinstance(g.get('joint'), dict) else '?'}, "
            f"AF {af:.3g}) and --pm2-max-af is {a.pm2_max_af}. Under the strict "
            f"'absent from gnomAD' reading PM2 cannot fire. If the project instead uses "
            f"'ultra-rare' as the PM2 condition, declare that threshold in Methods.")

    # --- in-silico codes
    rows, damaging, tolerated, indet = classify_predictors(vep)
    if a.scheme == "count3":
        available = damaging + tolerated + indet
        if len(available) >= 2 and len(damaging) >= 3 and not tolerated:
            codes.append({"code": "PP3", "direction": DAMAGING, "strength": "supporting",
                          "basis": f"{len(damaging)}/{len(available)} available tools damaging "
                                   f"({', '.join(damaging)})"})
        elif len(available) >= 2 and len(tolerated) >= 3 and not damaging:
            codes.append({"code": "BP4", "direction": BENIGN, "strength": "supporting",
                          "basis": f"{len(tolerated)}/{len(available)} available tools tolerated "
                                   f"({', '.join(tolerated)})"})
        else:
            notes.append(
                f"No PP3/BP4 on the count3 scheme: {len(damaging)} damaging, "
                f"{len(tolerated)} tolerated, {len(indet)} indeterminate, "
                f"{len(rows) - len(available)} missing.")
        if len(available) < 2:
            notes.append("Fewer than 2 predictors had data: ClinGen SVI says do not assign "
                         "PP3/BP4 on one tool.")
    else:
        r = vep.get("revel")
        if r is None:
            notes.append("revel scheme chosen but REVEL has no value for this position.")
        elif r >= 0.75:
            codes.append({"code": "PP3", "direction": DAMAGING, "strength": "supporting",
                          "basis": f"REVEL {r} >= 0.75"})
        elif r <= 0.15:
            codes.append({"code": "BP4", "direction": BENIGN, "strength": "supporting",
                          "basis": f"REVEL {r} <= 0.15"})
        else:
            notes.append(f"REVEL {r} sits in the indeterminate band (0.15-0.75): "
                         f"no PP3/BP4 on the calibrated scheme.")

    # --- functional codes (human-entered only)
    if a.functional == "damaging":
        codes.append({"code": "PS3", "direction": DAMAGING, "strength": "strong",
                      "basis": "functional assay reported a damaging effect (supplied by curator)"})
    elif a.functional == "benign":
        codes.append({"code": "BS3", "direction": BENIGN, "strength": "strong",
                      "basis": "functional assay reported no effect (supplied by curator)"})

    # --- lean
    path = [c for c in codes if c["direction"] == DAMAGING]
    ben = [c for c in codes if c["direction"] == BENIGN]
    if any(c["code"] == "BA1" for c in codes):
        lean = "Benign (standalone BA1) - curation complete, no further evidence needed"
    elif len(ben) > len(path):
        lean = "Benign-leaning"
    elif len(path) > len(ben):
        lean = ("Pathogenic-leaning, untested (no functional assay)" if a.functional == "none"
                else "Pathogenic-leaning")
    elif path and ben:
        lean = "Conflicting evidence - flag for review, do not force a call"
    else:
        lean = "VUS - insufficient evidence"

    print(f"Variant: {(d.get('clinvar') or {}).get('title')}")
    print(f"ClinVar: {(d.get('clinvar') or {}).get('classification')} "
          f"({(d.get('clinvar') or {}).get('review_status')})")
    print(f"gnomAD {g.get('dataset')}: AF={af} homozygotes={g.get('homozygotes')} "
          f"faf95_popmax={faf95} ({g.get('popmax_population')})")
    print("\nPredictor calls:")
    for r in rows:
        print(f"  {r['tool']:<9} {str(r['value']):<12} {r['call']}")
    print("\nCodes:")
    for c in codes:
        print(f"  {c['code']:<4} {c['direction']:<9} {c['strength']:<11} {c['basis']}")
    if not codes:
        print("  (none derivable)")
    print(f"\nLean: {lean}")
    if notes:
        print("\nNotes and caveats:")
        for n in notes:
            print(f"  - {n}")
    print(f"\nScheme: {a.scheme} | BS1>{a.bs1} | BA1>{a.ba1} | PM2 max AF {a.pm2_max_af} | "
          f"functional={a.functional}")
    print("Thresholds and dataset versions must be reported in Methods.")


if __name__ == "__main__":
    main()
