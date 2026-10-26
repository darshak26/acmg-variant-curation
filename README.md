# acmg-variant-curation

**Derive ACMG/AMP 2015 codes for a germline variant from ClinVar, gnomAD, Ensembl VEP and
PubMed — with the number that triggered every code printed beside it.**

[![CI](https://github.com/darshak26/acmg-variant-curation/actions/workflows/validate.yml/badge.svg)](https://github.com/darshak26/acmg-variant-curation/actions/workflows/validate.yml)
[![Release](https://img.shields.io/github/v/release/darshak26/acmg-variant-curation)](https://github.com/darshak26/acmg-variant-curation/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Agent Skills](https://img.shields.io/badge/Agent%20Skills-spec%20compliant-6f42c1.svg)](https://agentskills.io/specification)

An [Agent Skills](https://agentskills.io/specification) package: works as a skill in Claude Code,
Codex CLI, Hermes and Antigravity. Stdlib-only Python, no install step, no API keys.

> **This is evidence assembly, not a clinical classification.** It has not been validated against
> a ClinGen expert-panel set. See [Limitations](#limitations) first.

---

An LLM asked "what is the ACMG classification of this variant" answers from training data. The
number sounds precise and cannot be checked. This tool converts that into a question you *can*
check: **did every number come from a live query, and does it clear the threshold printed beside
it?**

Two properties make that enforceable rather than aspirational:

- Every code prints its basis — `PP3 · 4/4 available tools damaging (REVEL, CADD, BayesDel, SIFT)`.
- A **missing** predictor is `NA`, excluded from numerator *and* denominator. It never votes
  damaging, never votes tolerated. NA is not zero.

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [CLI reference](#cli-reference)
- [How it works](#how-it-works)
- [ACMG codes](#acmg-codes)
- [The lean](#the-lean)
- [The identifier guard](#the-identifier-guard)
- [Compatibility](#compatibility)
- [Project structure](#project-structure)
- [Development](#development)
- [Contributing](#contributing)
- [Limitations](#limitations)
- [Versioning](#versioning)
- [Citation](#citation)
- [License](#license)

## Installation

Requires **Python 3.10+** and outbound HTTPS. Nothing else — no `pip install`, no API key, no
config file.

The repository root *is* the skill. Copy the directory to your host's skill location:

```bash
git clone https://github.com/darshak26/acmg-variant-curation.git
cp -r acmg-variant-curation ~/.claude/skills/        # Claude Code
```

| Host | Skill location | Invoke |
|---|---|---|
| Claude Code | `~/.claude/skills/acmg-variant-curation/` or `.claude/skills/…` in a repo | `/acmg-variant-curation` |
| Codex CLI | `~/.agents/skills/acmg-variant-curation/` — **not** `~/.codex/skills` | `/skills` or `$acmg-variant-curation` |
| Hermes | `~/.hermes/profiles/<profile>/skills/research/acmg-variant-curation/` | `/acmg-variant-curation` |
| Antigravity | inside a plugin bundle: `agy plugin import claude` | `/skills` in-session |

Not using an agent? Run the scripts directly — they are standalone CLI tools.

> **Codex CLI note:** its default sandbox is `workspace-write` with network access **off**. Every
> step here is a live HTTPS call, so the skill does nothing until you enable network for the
> sandbox.

## Quick start

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

# 4. read PubMed by hand — the payload's action_required explains why

# 5. derive codes
python3 scripts/acmg_codes.py G245S.json --scheme count3
```

Real output for TP53 p.Gly245Ser, from [`examples/g245s.json`](examples/g245s.json):

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

Points (ClinGen SVI): pathogenic 1 [PP3 1]  |  benign 0
Lean: VUS - insufficient evidence (1 pathogenic / 0 benign points).
      Direction pathogenic, which is not a classification
Verdict: VUS
```

Read those last three lines first. **This variant is Pathogenic by ClinGen expert panel and the tool
still says VUS** — deliberately. Richards 2015 has no combination rule under which one Supporting
code classifies anything. Calling it "pathogenic-leaning" because one code outnumbered zero would
be flattering, not accurate. The tool reports direction, states plainly that direction is not a
classification, and names what would close the gap: PS4 case data or a read PS3 assay, neither of
which is derivable from these databases.

Two details also worth reading. **PM2 did not fire** — the variant is in gnomAD (AC=5), so under a
strict "absent from gnomAD" reading it cannot, and the tool says so instead of quietly rounding five
alleles down to zero. And the count is **4/4, not 3/4** — BayesDel comes from a second VEP call, and
the payload records `vep.bayesdel_merge: "ok"`, because a failed merge there would shrink the
denominator and could flip the call.

## CLI reference

### `scripts/fetch_evidence.py`

| Flag | Required | Meaning |
|---|---|---|
| `--clinvar` | one of these | ClinVar variation ID, e.g. `12365` |
| `--gnomad` | one of these | gnomAD v4 ID, GRCh38: `chr-pos-ref-alt`, e.g. `17-7674230-C-T` |
| `--hgvs` | one of these | HGVS with accession — `NM_000546.6:c.733G>A` or `NC_000017.11:g.7674230C>T` |
| `--transcript` | no | MANE Select ENST, filters VEP to one isoform |
| `--pubmed-terms` | no | Comma-separated protein-change spellings, e.g. `"p.Gly245Ser,G245S"` |
| `--out` | no | JSON output path (default: stdout) |
| `--selftest` | no | Check all four endpoints and exit |

Exit codes: `0` success, `2` out of scope (CNV/SV).

### `scripts/acmg_codes.py`

| Flag | Default | Meaning |
|---|---|---|
| `jsonfile` | — | Payload from `fetch_evidence.py` |
| `--scheme` | `count3` | `count3` or `revel` — pick one, declare it, never mix |
| `--bs1` | `0.01` | BS1 frequency threshold; override with gene-specific prevalence |
| `--ba1` | `0.05` | BA1 standalone-frequency threshold |
| `--pm2-max-af` | `0.0` | PM2 "absent" threshold. Raising it is a methodological choice — declare it |
| `--functional` | `none` | `damaging` / `benign`, from an assay **you read**. Never inferred |

## How it works

Order matters — population frequency can end curation early (BA1), and functional evidence outranks
every in-silico predictor.

| Step | Source | Yields |
|---|---|---|
| 1 | ClinVar (NCBI E-utilities) | classification, review status, last evaluated, variant type |
| 2 | gnomAD v4 (GraphQL) | AF, grpmax faf95, homozygotes, absence |
| 3 | Ensembl VEP (REST) | REVEL, CADD, BayesDel, PolyPhen, SIFT, AlphaMissense |
| 4 | PubMed (E-utilities) | whether a functional assay exists at all |
| 5 | `acmg_codes.py` | the codes, each with its basis |

**Binaries, not just variants.** `PS3`/`BS3` are never inferred. PubMed returns candidate papers and
you read them; a title is not a functional assay.

## ACMG codes

Seven codes are derivable from public data:

| Code | Direction | Strength | Source | Fires when |
|---|---|---|---|---|
| BA1 | Benign | Standalone | gnomAD | grpmax AF > 0.05 — curation ends |
| BS1 | Benign | Strong | gnomAD | AF > prevalence (default 0.01) |
| PM2 | Pathogenic | Supporting | gnomAD | absent, **and the identifier verified** |
| PP3 | Pathogenic | Supporting | VEP | ≥3 available predictors damaging |
| BP4 | Benign | Supporting | VEP | ≥3 available predictors tolerated |
| PS3 | Pathogenic | Strong | literature | assay shows loss of function (human-read) |
| BS3 | Benign | Strong | literature | assay shows no effect (human-read) |

Predictor thresholds: REVEL ≥ 0.5 · CADD PHRED ≥ 20 · BayesDel addAF ≥ 0.13 (≤ −0.36 tolerated,
between = indeterminate, votes neither) · PolyPhen D/P · SIFT deleterious.

**Not covered, deliberately:** PS4, PM6, PP1 and PS2 need case-level data this pipeline cannot see.
CNVs, structural variants and somatic curation are out of scope — the PP3/BP4 machinery means
nothing there. Clinical use needs a ClinGen expert panel.

Out-of-scope input is **rejected, not half-curated** — exit code 2, the variant type and size, and a
pointer to the right standard:

```
[OUT OF SCOPE] Variant 150702 is a copy number gain, ~126 kb: GRCh38/hg38 17p13.1(chr17:7544902-7670581)x3
acmg-variant-curation curates germline SNVs and small indels under Richards 2015.
For CNVs use the ACMG/ClinGen 2020 technical standards (Riggs et al., Genet Med 22:245-257) ...
```

### Two in-silico schemes

PP3 has two defensible implementations that disagree on real variants. Pick one, declare it in
Methods, never mix.

- **`count3`** (default) — PP3 when ≥3 of the *available* predictors call damaging. Works without a
  calibrated meta-predictor.
- **`revel`** — single calibrated predictor: PP3 at REVEL ≥ 0.75, BP4 at ≤ 0.15, middle is
  indeterminate. The stricter ClinGen SVI 2022 reading.

## The lean

Codes carry direction *and strength*. Strength converts to ClinGen SVI points (Tavtigian et al.
2020): Very Strong 8, Strong 4, Moderate 2, Supporting 1. **Pathogenic at ≥ 10, benign at ≥ 7**,
otherwise VUS. BA1 stays standalone.

**The ceiling is a design property, not a bug.** The seven public-data codes max out at **6 points**
(PM2 1 + PP3 1 + PS3 4), so this tool cannot reach a pathogenic classification on its own. A test
pins that so it can't be "fixed" into over-calling. Expect `VUS` with a direction on almost
everything — that is the honest answer for public data.

This replaced a code-count lean that called any variant with more pathogenic codes than benign ones
"pathogenic-leaning", which meant a variant merely absent from gnomAD and nothing else read as
pathogenic — synonymous, intronic and UTR variants alike. When the codes and ClinVar disagree, the
conflict is printed as `CONFLICT WITH CLINVAR` rather than resolved silently in the tool's favour.

## The identifier guard

The failure this tool exists to prevent is not a wrong score — it is a **right score attached to the
wrong variant**. ClinVar will happily return a pathogenic classification for variant A while gnomAD
and VEP answer for variant B, and every code still fires.

`fetch_evidence.py` cross-checks the identifiers and emits a `consistency` block:

```json
{
  "ok": false,
  "problems": ["gnomAD id is at 17:7676154 but VEP resolved 17:7674230 -> DIFFERENT POSITION"],
  "checks": ["requested HGVS == ClinVar HGVS (c.733G>A)"],
  "identifiers": { "gnomad_id": "17-7676154-C-T", "vep_locus": "17:7674230 G/A" }
}
```

`acmg_codes.py` then **refuses PM2** on that payload. This matters because gnomAD answers
`Variant not found` identically for a typo and for a genuinely ultra-rare allele — absence alone
proves nothing about a variant nobody verified.

A guard that compares nothing must not report success. If no cross-check is possible, `ok` is
`false` and PM2 is withheld; `consistency.checks` lists the comparisons actually performed.

## Compatibility

| Host | Loads | Notes |
|---|---|---|
| Claude Code | ✅ | Follows the Agent Skills spec; loads from personal, project or plugin level |
| Codex CLI | ✅ | `~/.agents/skills/`; needs sandbox network access enabled |
| Hermes | ✅ | Auto-registers `/<name>` as a slash command |
| Antigravity | ✅ | Via a Claude-format plugin bundle |

Frontmatter uses only the six spec fields (`name`, `description`, `license`, `compatibility`,
`metadata`, `allowed-tools`). No host-specific extensions — those hard-error on some packaging
paths. Paths inside `SKILL.md` are relative to the skill root.

Run `python3 scripts/validate_skill.py` to check this yourself.

## Project structure

```
acmg-variant-curation/
├── SKILL.md                  skill definition: rules, codes, pitfalls, deliverable
├── README.md                 this file
├── CHANGELOG.md              version history; threshold/guard changes are breaking
├── CITATION.cff              machine-readable citation
├── LICENSE · NOTICE          MIT + scope disclaimer
├── scripts/
│   ├── fetch_evidence.py     ClinVar + gnomAD + VEP + PubMed → one payload
│   ├── acmg_codes.py         payload → codes, points, lean
│   └── validate_skill.py     spec conformance + host compatibility
├── references/
│   └── thresholds.md         full tables, scheme comparison, worked examples
├── examples/
│   └── g245s.json            real fetched payload, used as test fixture
├── tests/
│   └── test_engine.py        offline assertions, no network
└── .github/workflows/
    └── validate.yml          CI: Python 3.10 + 3.12
```

## Development

```bash
python3 scripts/validate_skill.py              # spec conformance + compatibility table
python3 tests/test_engine.py                   # 19 offline assertions, no network
python3 scripts/fetch_evidence.py --selftest   # network: are all four endpoints up?
```

CI runs the first two on Python 3.10 and 3.12. The live-endpoint check is manual
(`workflow_dispatch`) because Ensembl VEP flaps roughly half the time and would make CI red for
reasons unrelated to the code.

**After editing either script, prove them with three payloads**, not one — a happy path exercises
no guard:

| Case | Expect |
|---|---|
| all identifiers correct | `consistency.ok = true`, codes print normally |
| wrong `--clinvar` | `ok = false`, `WRONG VARIATION ID` on stderr |
| wrong `--gnomad` position | `ok = false`, `DIFFERENT POSITION`, PM2 withheld |

Allow 1–2 minutes per payload: the elapsed time is VEP retries, not a hang.

## Contributing

Issues and pull requests are welcome. Please:

1. **Add a test** for any change to the guard or the lean — those are the behaviours that matter.
2. **Keep the scripts stdlib-only.** CI fails on a non-stdlib import; `pip install` defeats the
   point of a portable skill.
3. **Never widen a threshold without a citation** in `references/thresholds.md`.
4. Update `CHANGELOG.md` and bump `metadata.version`. A change to thresholds, schemes or guard
   behaviour is a **breaking** change.

## Limitations

- **Not validated against expert calls.** Exercised on a handful of TP53 variants. No concordance
  study against a ClinGen panel-classified set has been run. Treat output as a draft for review, not
  a result.
- **Cannot reach a pathogenic call.** The seven codes max out at 6 of the 10 points required. This
  is deliberate; see [The lean](#the-lean).
- **PS3/BS3 are never inferred.** You read the papers.
- **Ensembl VEP flaps.** Identical requests return 200, 500 and 503 within seconds — roughly one
  success in two, measured 2026-10-26. The client retries 5xx six times with jitter, so a run takes
  1–2 minutes. A 500 is not a bad parameter.
- **Public data only.** No somatic, no CNV/SV, no phenotype, no segregation, no case-level data.
- **Codes can change between versions.** The lean and guard are part of the contract; a v1.0.0 and
  a v2.0.0 run on the same variant may disagree.

## Versioning

[Semantic versioning](https://semver.org/). The version lives in `metadata.version` in `SKILL.md`
and in [CHANGELOG.md](CHANGELOG.md). Pin it in your write-up alongside the dataset and access date:

> Variants were annotated against ClinVar, gnomAD (dataset `gnomad_r4`, accessed YYYY-MM-DD) and
> Ensembl VEP (transcript NM_000546.6/ENST00000269305). In-silico evidence was assessed on
> algorithms with annotation for that position; algorithms without data were excluded and recorded
> as NA. PP3/BP4 was applied on the `count3` scheme. Lean on ClinGen SVI points (Tavtigian et al.
> 2020). Tool: acmg-variant-curation v2.0.0.

## Citation

> Padmane S. *acmg-variant-curation* (v2.0.0): ACMG code derivation from public variant databases
> as an agent skill. https://github.com/darshak26/acmg-variant-curation

Machine-readable metadata in [`CITATION.cff`](CITATION.cff).

## License

[MIT](LICENSE) — see also [NOTICE](NOTICE) for the scope disclaimer.

## Acknowledgements

ACMG/AMP 2015 criteria: Richards et al., *Genet Med* 17:405–424. Points weighting: Tavtigian et al.,
*Genet Med* 22:1653–1657 (2020). CNV standards cited in out-of-scope messages: Riggs et al.,
*Genet Med* 22:245–257 (2020).

Data from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), [gnomAD](https://gnomad.broadinstitute.org/),
[Ensembl VEP](https://www.ensembl.org/Tools/VEP) and [PubMed](https://pubmed.ncbi.nlm.nih.gov/).
Not affiliated with or endorsed by any of them.
