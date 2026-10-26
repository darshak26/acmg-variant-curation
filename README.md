# acmg-variant-curation

An [Agent Skills](https://agentskills.io/specification) package that derives ACMG/AMP 2015 codes
for a germline SNV/indel from four public sources — **ClinVar, gnomAD, Ensembl VEP, PubMed** — and
prints the number that triggered every code.

Built for a VUS-curation project; gene-agnostic. Works as a skill in Claude Code, Codex CLI,
Hermes, and Antigravity.

**This is evidence assembly, not a clinical classification.** See [Scope](#scope) for exactly what
it does and does not cover.

## Why it exists

An LLM asked "what is the ACMG classification of this variant" answers from training data. The
number sounds precise and cannot be checked. This skill converts that into a question you *can*
check: **did every number come from a live query in this session, and does it clear the threshold
printed beside it?**

Two things make that enforceable rather than aspirational:

- The scripts refuse to print a code without its basis line (`PP3 · 3/3 available tools damaging`).
- A **missing** predictor is `NA`, which is excluded from numerator *and* denominator. It never
  votes damaging, never votes tolerated. NA is not zero.

## Quick start

Requires Python **3.10+** (stdlib only, nothing to `pip install`) and outbound HTTPS.
Nothing else.

```bash
# 0. confirm all four endpoints answer before trusting a batch
python3 scripts/fetch_evidence.py --selftest

# 1-3. ClinVar + gnomAD + VEP for one variant (three identifiers, same variant)
python3 scripts/fetch_evidence.py \
    --clinvar 12365 \
    --gnomad 17-7674230-C-T \
    --hgvs "NM_000546.6:c.733G>A" --transcript ENST00000269305 \
    --pubmed-terms "p.Gly245Ser,G245S" \
    --out G245S.json

# 4. read PubMed by hand - the payload's action_required says why

# 4b. trust nothing until this says ok: true
python3 -c "import json;c=json.load(open('G245S.json'))['consistency'];print(c['ok'],c['problems'])"

# 5. derive codes
python3 scripts/acmg_codes.py G245S.json --scheme count3
```

Real output for TP53 p.Gly245Ser (in [`examples/g245s.json`](examples/g245s.json)):

```
Variant: NM_000546.6(TP53):c.733G>A (p.Gly245Ser)
ClinVar: Pathogenic (reviewed by expert panel)
gnomAD gnomad_r4: AF=3.098419434278186e-06 homozygotes=0 faf95_popmax=6.8e-07 (nfe)

Predictor calls:
  REVEL     0.947        damaging
  CADD      27.5         damaging
  BayesDel  0.550935     damaging
  SIFT      deleterious  damaging

Codes:
  PP3  pathogenic supporting  4/4 available tools damaging (REVEL, CADD, BayesDel, SIFT)

Lean: Pathogenic-leaning, untested (no functional assay)
```

Two things are worth reading closely. **PM2 did not fire** — the variant is in gnomAD (AC=5), so
under a strict "absent from gnomAD" reading it cannot, and the tool says so rather than quietly
rounding 5 alleles down to zero. And the count is **4/4, not 3/4** — BayesDel came from the second
VEP call, which is why the payload records `vep.bayesdel_merge: "ok"`. A failed merge there would
shrink the denominator and could flip the call.

## The four sources, in order

Order matters — population frequency can end curation early (BA1), and functional evidence
outranks every in-silico predictor.

| Step | Source | What it yields |
|---|---|---|
| 1 | ClinVar (NCBI E-utilities) | classification, review status, last evaluated |
| 2 | gnomAD v4 GraphQL | AF, grpmax faf95, homozygotes, absence |
| 3 | Ensembl VEP REST | REVEL, CADD, BayesDel, PolyPhen, SIFT, AlphaMissense |
| 4 | PubMed E-utilities | whether a functional assay exists at all |
| 5 | `acmg_codes.py` | the codes, each with its basis |

## Scope

**Covered** — seven codes derivable from public databases:

| Code | Direction | Strength | Source | Fires when |
|---|---|---|---|---|
| BA1 | Benign | Standalone | gnomAD | grpmax AF > 0.05 — curation ends |
| BS1 | Benign | Strong | gnomAD | AF > prevalence (default 0.01) |
| PM2 | Pathogenic | Supporting | gnomAD | absent, **and the identifier verified** |
| PP3 | Pathogenic | Supporting | VEP | ≥3 available predictors damaging |
| BP4 | Benign | Supporting | VEP | ≥3 available predictors tolerated |
| PS3 | Pathogenic | Strong | literature | assay shows loss of function (human-read) |
| BS3 | Benign | Strong | literature | assay shows no effect (human-read) |

**Not covered, deliberately** — PS4, PM6, PP1, PS2 need case-level data this pipeline cannot see.
CNVs, structural variants and somatic curation are out of scope: the PP3/BP4 machinery means
nothing there. Clinical use needs a ClinGen expert panel; a curation from this tool is a proposal.

## Two in-silico schemes

PP3 has two defensible implementations that disagree on real variants. Pick one, declare it in
Methods, never mix.

- `count3` (default) — PP3 when ≥3 of the **available** predictors call damaging. Works without a
  calibrated meta-predictor.
- `revel` — single calibrated predictor: PP3 at REVEL ≥ 0.75, BP4 at ≤ 0.15, middle is
  indeterminate. The stricter ClinGen SVI 2022 reading.

## The identifier guard

The failure this tool was built to prevent is not a wrong score — it is a **right score attached to
the wrong variant**. ClinVar will happily return a pathogenic classification for variant A while
gnomAD and VEP answer for variant B, and every code still fires.

`fetch_evidence.py` therefore cross-checks all three identifiers and emits a `consistency` block:

```json
{
  "ok": false,
  "problems": ["gnomAD id is at 17:7676154 but VEP resolved 17:7674230 -> DIFFERENT POSITION"],
  "identifiers": { "requested_hgvs": "c.733G>A", "clinvar_hgvs": "c.733G>A", "gnomad_id": "17-7676154-C-T", "vep_locus": "17:7674230 G/A" }
}
```

`acmg_codes.py` then **refuses PM2** on that payload. This matters because gnomAD answers
`Variant not found` identically for a typo and for a genuinely ultra-rare allele — so absence alone
cannot prove anything about a variant nobody verified.

## Installation

The repository root *is* the skill. Copy the whole directory to your host's skill location.

| Host | Location | Invoke |
|---|---|---|
| Claude Code | `~/.claude/skills/acmg-variant-curation/` or `.claude/skills/acmg-variant-curation/` in a repo | `/acmg-variant-curation` |
| Codex CLI | `~/.agents/skills/acmg-variant-curation/` (not `~/.codex/skills`) or `<repo>/.agents/skills/` | `/skills` or `$acmg-variant-curation` |
| Hermes | `~/.hermes/profiles/<profile>/skills/research/acmg-variant-curation/` | `/acmg-variant-curation` |
| Antigravity | inside a plugin bundle: `agy plugin import claude` | `/skills` in-session |

The Python floor above is tested on 3.10 and 3.12 in CI, not assumed.

```bash
git clone https://github.com/darshak26/acmg-variant-curation.git
cp -r acmg-variant-curation ~/.claude/skills/
```

**Codex CLI caveat:** its default sandbox is `workspace-write` with network access **off**. Every
step here is a live HTTPS call, so the skill does nothing until you enable network access for the
sandbox.

## Verify it yourself

```bash
python3 scripts/validate_skill.py   # spec conformance + host compatibility table
python3 tests/test_engine.py        # offline: does PM2 still refuse a bogus ID?
python3 scripts/fetch_evidence.py --selftest   # network: are all four endpoints up?
```

`validate_skill.py` checks frontmatter against the spec, rejects host-specific fields (which some
packaging paths hard-error on), and proves every script imports from the standard library only.

## Versioning

Version lives in `metadata.version` in `SKILL.md`; see [CHANGELOG.md](CHANGELOG.md). Any change to
thresholds, schemes, or the guard's behaviour is a breaking change — a curation's codes must be
reproducible from the version recorded in Methods.

Pin the version in your write-up alongside the dataset and access date:

> Variants were annotated against ClinVar, gnomAD (dataset `gnomad_r4`, accessed YYYY-MM-DD) and
> Ensembl VEP (transcript NM_000546.6/ENST00000269305). In-silico evidence was assessed on
> algorithms with annotation for that position; algorithms without data were excluded and recorded
> as NA. PP3/BP4 was applied on the `count3` scheme. Tool: acmg-variant-curation v1.0.0.

## Honest limitations

- **Not validated against expert calls.** It has been exercised on a handful of TP53 variants. No
  concordance study against a ClinGen panel-classified set has been run. Until that exists, treat
  the output as a draft for review, not a result.
- **PS3/BS3 are never inferred.** PubMed returns candidate papers; you read them. A title is not a
  functional assay, and a gene-only query returns thousands of hits of noise.
- **Ensembl VEP flaps.** Identical requests return 200, 500 and 503 within seconds — roughly one
  success in two, measured 2026-10-26. The client retries 5xx six times with jitter, so a run takes
  1–2 minutes. A 500 is not a bad parameter.
- **Public data only.** No somatic, no CNV/SV, no phenotype, no segregation, no case-level data.

## Citation

> Padmane S. *acmg-variant-curation* (v1.0.0): ACMG code derivation from public variant databases
> as an agent skill. https://github.com/darshak26/acmg-variant-curation

## License

MIT — see [LICENSE](LICENSE). ACMG/AMP 2015 criteria: Richards et al., *Genet Med* 17:405–424.
