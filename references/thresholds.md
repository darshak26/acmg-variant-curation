# Thresholds, schemes, and a worked example

Every number here was fetched live from this box; each block carries its own access date.
Do not carry a number forward from this file into a new curation — re-fetch it.

## Frequency thresholds (gnomAD)

| Code | Direction | Strength | Threshold | Note |
|---|---|---|---|---|
| BA1 | Benign | Standalone | grpmax AF > 0.05 | Fires alone; curation stops |
| BS1 | Benign | Strong | AF > gene-specific prevalence | Default 0.01; TP53 panel threshold overrides |
| PM2 | Pathogenic | Supporting | Absent (AF = 0) | "Ultra-rare" needs a declared cutoff |

Use `faf95.popmax` for the grpmax comparison; use `ac/an` for the plain AF. gnomAD v4's
`joint` block exposes `ac`, `an`, `homozygote_count`, `faf95{popmax, popmax_population}` —
**not** `af`.

## Predictor thresholds

| Tool | Damaging | Tolerated | Indeterminate | VEP field |
|---|---|---|---|---|
| REVEL | ≥ 0.5 | < 0.5 | — | `revel` (with `REVEL=1`) |
| CADD PHRED | ≥ 20 | < 20 | — | `cadd_phred` (with `CADD=1`) |
| BayesDel addAF | ≥ 0.13 | ≤ −0.36 | −0.36 … 0.13 | `bayesdel_addaf_score` (dbNSFP param) |
| PolyPhen | D / P | B | — | `polyphen_prediction` (with `PolyPhen=1`) |
| SIFT | deleterious | tolerated | — | `sift_prediction` (default, no flag) |
| AlphaMissense | likely_pathogenic | likely_benign | ambiguous | `alphamissense.am_class` |

BayesDel and AlphaMissense are not ACMG codes on their own — they are extra voices. The
count3 scheme counts BayesDel; AlphaMissense is reported but not counted by either scheme.

## The two PP3/BP4 schemes

| | count3 (project default) | revel (ClinGen SVI 2022 style) |
|---|---|---|
| Condition | ≥3 available tools damaging, none tolerated | REVEL ≥ 0.75 |
| Denominator | only tools with data; NA excluded | REVEL only |
| Minimum data | ≥2 tools | REVEL present |
| Effect on D393V | PP3 fires (4/4 damaging) | PP3 does not fire (0.69 is mid-band) |

That divergence is the point of declaring a scheme. Same variant, same scores, two honest
answers. A paper that switches schemes between variants is not reproducible.

## Worked example — TP53 p.Asp393Val, NM_000546.6:c.1178A>T

Fetched values (gnomad_r4, accessed 2026-10-24):

| Source | Value |
|---|---|
| ClinVar 3634421 | VCV003634421 · Uncertain significance · criteria provided, single submitter |
| gnomAD joint | ac 2 / an 1,613,494 → AF 1.24e-06 · homozygotes 0 · faf95 popmax 3.65e-06 (sas) |
| VEP | REVEL 0.69 · CADD 23.4 · BayesDel addAF 0.312 · SIFT deleterious · AlphaMissense likely_pathogenic |
| VEP | PolyPhen absent (no coverage at this position — verified, not a bug) |
| PubMed | 0 hits for the exact protein change under `[tiab]` phrase search |

Derived:

- **PP3 fires** under count3 (4 of 4 available damaging). Does not fire under revel.
- **PM2 does not fire.** Two alleles exist, so the variant is not absent. Under AF = 0 the
  condition is false. It would need a declared ultra-rare cutoff (e.g. AF < 1e-5) to apply —
  and then the Methods must say so.
- PS3/BS3: no published assay for this change → no functional code.

Lean: pathogenic-leaning, untested. Two versions of the write-up are defensible; what is
not defensible is claiming PM2 while reporting two carriers in the same table.

## Worked example — TP53 p.Gly245Ser, NM_000546.6:c.733G>A

The same variant run twice: once with every identifier correct, once with a plausible but
wrong gnomAD id. The scores are identical because VEP was never the problem; only the
frequency source changed, and that is enough to invent a pathogenic code.

Correct identifiers (gnomad_r4, accessed 2026-10-26): ClinVar 12365, gnomAD
`17-7674230-C-T`, ClinVar SPDI `NC_000017.11:7674229:C:T` (0-based → 1-based 7674230).

| Source | Value |
|---|---|
| ClinVar 12365 | `c.733G>A (p.Gly245Ser)` · Pathogenic · reviewed by expert panel |
| gnomAD r4 joint | ac 5 · AF 3.10e-06 · homozygotes 0 · faf95 popmax 6.8e-07 (nfe) |
| VEP | REVEL 0.947 · CADD 27.5 · SIFT deleterious · AlphaMissense likely_pathogenic · PolyPhen absent |
| PubMed | 54 hits for `p.Gly245Ser`/`G245S` — titles only, none read |

Derived: **PP3 fires under both schemes** (REVEL 0.947 ≥ 0.75; 3–4 of available tools
damaging). **PM2 does not fire** — five carriers exist, so "absent" is false, and no
ultra-rare cutoff was declared. PS3/BS3: no assay read, so none assigned. Lean:
pathogenic-leaning, untested. ClinVar's expert-panel call is context, not a code.

### The trap: same variant, wrong gnomAD id

Substituting `17-7676154-C-T` — a real TP53 position, 1,924 bp away — produced
`absent: true` and a PM2. Nothing else in the pipeline objected. gnomAD answers
`Variant not found` identically for a wrong position and for a genuine absence, so the
payload is internally consistent and reads as "ultra-rare pathogenic variant" when the
variant is not even the one being curated. The consistency guard now blocks this by
comparing the id against VEP's resolved locus; PM2 additionally requires `ok: true`.

Lesson worth carrying past this gene: frequency databases treat "not found" as a normal
answer, so an unverified identifier converts a lookup failure into evidence.

## Sources for the thresholds

- ACMG/AMP 2015 standards (Richards et al.) — code definitions and strengths
- ClinGen SVI recommendations — PP3/BP4 calibration, PM2 downgrade, NA handling
- gnomAD v4 GraphQL, Ensembl VEP REST (REST endpoints verified live, this box)
