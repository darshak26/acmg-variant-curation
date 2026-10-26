---
name: acmg-variant-curation
description: Use when deriving ACMG codes for a variant from ClinVar, gnomAD, VEP and PubMed, or curating a VUS. Fetches all four sources live and prints the basis for every code.
license: MIT
compatibility: Requires python3 (stdlib only, no pip install) and outbound HTTPS to eutils.ncbi.nlm.nih.gov, gnomad.broadinstitute.org and rest.ensembl.org. A sandbox that blocks network access disables the whole skill.
metadata:
  version: "2.0.0"
  author: Swapnil Padmane
  repository: https://github.com/darshak26/acmg-variant-curation
  spec: https://agentskills.io/specification
---

# ACMG variant curation

Derive an evidence-backed ACMG lean for a germline SNV/indel using four public sources
and the ACMG/AMP 2015 rulebook. Built for the TP53 project but gene-agnostic.

Scope, stated honestly: the seven codes below are the ones derivable from public
databases. Curated criteria — PS4, PM6, PP1, PS2 — need case-level data and are **not**
covered. This skill does not handle CNVs, structural variants or somatic curation; the
PP3/BP4 machinery means nothing there. Out-of-scope input is rejected with exit code 2 and a
pointer to the right standard, never half-curated.

Five steps, in this order: **ClinVar → gnomAD → VEP → PubMed → ACMG**. Order matters —
population frequency can end curation early (BA1), and functional evidence outranks every
predictor.

## Hard rules

1. **Never state a score you did not fetch in this session.** No recalled REVEL values, no
   remembered frequencies. Run the script, quote its output. An identifier is the same kind
   of value: resolve the ClinVar id and *derive* the gnomAD id from ClinVar's own SPDI — a
   coordinate written out by hand is a guess that looks like data.
2. **NA is not zero.** A missing predictor is excluded from numerator *and* denominator. It
   never votes damaging, never votes tolerated.
3. **PS3/BS3 require a paper you actually read.** A title match is not a functional assay.
   If the exact protein change returns zero hits, the honest answer is "no functional data".
4. **Pin the versions.** gnomAD dataset (`gnomad_r4`), transcript accession, predictor score
   versions, and the date fetched. Unpinned frequencies are unreproducible.
5. **Declare your scheme.** PP3 has two defensible implementations (below) that disagree on
   real variants. Pick one, say which in Methods, never mix.
6. **Check `consistency.ok` before you read any code.** The three identifiers must be
   proven to describe one variant. `fetch_evidence.py` now emits a `consistency` block;
   `ok: false` means the payload must not be curated. `acmg_codes.py` refuses PM2 on it.
7. **A curation is a proposal, not a classification.** VUS work in a student project is
   evidence assembly; clinical use needs a ClinGen expert panel.

## Run it

Paths below are relative to this skill's root, per the Agent Skills spec, so the block
runs unchanged on any host. Hermes installs the skill at
`~/.hermes/profiles/athena/skills/research/acmg-variant-curation/` — `cd` there first if
your host does not make the skill root the working directory.

```bash
# 0. confirm all four endpoints answer before trusting a batch
python3 scripts/fetch_evidence.py --selftest

# 1-3. fetch ClinVar + gnomAD + VEP (same variant, three identifiers)
python3 scripts/fetch_evidence.py \
    --clinvar 3634421 \
    --gnomad 17-7669613-T-A \
    --hgvs "NM_000546.6:c.1178A>T" --transcript ENST00000269305 \
    --pubmed-terms "p.Asp393Val,D393V" \
    --out D393V.json

# 4. read PubMed by hand (see the return payload's action_required)

# 4b. trust nothing until this says ok: true
python3 -c "import json;c=json.load(open('D393V.json'))['consistency'];print(c['ok'],c['problems'])"

# 5. derive codes
python3 scripts/acmg_codes.py D393V.json --scheme count3
```

Identifiers must describe the same variant: ClinVar variation ID, GRCh38 gnomAD id
(`chr-pos-ref-alt`), HGVS c. with accession, and the MANE Select ENST.

### Resolve identifiers, never compose them

The guard catches bad identifiers; it cannot supply good ones. Get each id from a source of
record, in this order — resolving takes ~15 s and removes the whole failure class.

```bash
# 1. variation id from a tight gene + protein-change query (expect 1 hit, not 50)
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&retmode=json&term=TP53%5Bgene%5D+AND+G245S"

# 2. canonical GRCh38 SPDI for that id
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&retmode=json&id=12365" \
 | python3 -c "import json,sys; r=json.load(sys.stdin)['result']['12365'];
print(r['title']); print(r['variation_set'][0]['canonical_spdi'])"
# NM_000546.6(TP53):c.733G>A (p.Gly245Ser)
# NC_000017.11:7674229:C:T     <- 0-based position, ref, alt
# gnomAD id = 17-7674230-C-T  <- position + 1, alleles unchanged
```

- Read the ClinVar `title` and confirm the protein change is the one you meant *before*
  using its id. `esearch` on a protein change can return a neighbouring codon.
- Derive the gnomAD id from the SPDI: **0-based position + 1**, alleles as printed. Never
  count base pairs from a remembered codon number.
- Check the accession build: `NC_000017.11` is GRCh38 (gnomAD v4). A `.10` accession is
  GRCh37 and its coordinates will be shifted.
- Take the `--hgvs` accession from the same ClinVar title, and the ENST from the transcript
  the record actually annotates — do not pair a fresh id with a cached transcript.

### What the consistency guard proves, and what it does not

`check_identity()` compares three things: ClinVar's title HGVS against `--hgvs`, the
gnomAD id's chromosome and position against VEP's resolved locus, and the gnomAD alleles
against VEP's `allele_string` (complemented when `strand == -1`). That covers a wrong
variation ID, a wrong coordinate, and a wrong alt allele.

It cannot catch a wrong ClinVar id when `--hgvs` is omitted. Pass all four identifiers or
the guard is blind on that axis.

### After editing either script, prove them this way

A happy-path run exercises no guard — it only proves the endpoints answered. Run three
payloads and check each one:

| Case | Expect |
|---|---|
| all identifiers correct | `consistency.ok = true`; codes print normally |
| wrong `--clinvar`, rest correct | `ok = false`, stderr `WRONG VARIATION ID` |
| wrong `--gnomad` position | `ok = false`, stderr `DIFFERENT POSITION`, and PM2 withheld with an explanatory note |

Then `python3 -m py_compile` both scripts and re-run `--selftest`. Allow 1–2 minutes per
payload: the elapsed time is VEP retries, not a hang.

## The 7 derivable codes

| Code | Direction | Strength | Source | Fires when |
|---|---|---|---|---|
| BA1 | Benign | Standalone | gnomAD | grpmax AF > 0.05 → curation ends immediately |
| BS1 | Benign | Strong | gnomAD | AF > gene-specific prevalence (default 0.01) |
| PM2 | Pathogenic | Supporting | gnomAD | Absent (AF = 0). "Ultra-rare" needs a declared threshold |
| PP3 | Pathogenic | Supporting | VEP | ≥3 available predictors damaging, none tolerated |
| BP4 | Benign | Supporting | VEP | ≥3 available predictors tolerated, none damaging |
| PS3 | Pathogenic | Strong | literature | Functional assay shows loss of function (human-read) |
| BS3 | Benign | Strong | literature | Functional assay shows no effect (human-read) |

Predictor thresholds: REVEL ≥ 0.5 damaging / < 0.5 tolerated · CADD PHRED ≥ 20 / < 20 ·
BayesDel addAF ≥ 0.13 damaging / ≤ −0.36 tolerated (between = indeterminate, votes neither) ·
PolyPhen D or P damaging / B benign · SIFT deleterious / tolerated.

Two PP3 schemes:
- `count3` — project default, ≥3 of available tools damaging. Works without a calibrated
  meta-predictor; can fire when REVEL alone sits in the middle.
- `revel` — single calibrated predictor, PP3 at REVEL ≥ 0.75, BP4 at ≤ 0.15, middle is
  indeterminate. ClinGen SVI 2022 preference; stricter.

Decision logic — **points, never a code count** (ClinGen SVI / Tavtigian 2020 weights):
Very Strong 8, Strong 4, Moderate 2, Supporting 1. BA1 → benign, done. Pathogenic at
≥ 10 points, benign at ≥ 7, otherwise **VUS**. Direction is reported separately and is
not a classification.

The seven public-data codes max out at 6 pathogenic points (PM2 1 + PP3 1 + PS3 4), so
this skill **cannot** reach a pathogenic call alone. Do not "improve" the lean to fix
that: the gap is closed by case data (PS4) or a read assay (PS3), not by re-weighting.
When the codes and ClinVar point different ways, print `CONFLICT WITH CLINVAR` and do not
resolve it in the tool's favour.

## Pitfalls (all verified against the live APIs)

- **gnomAD v4 `joint` has no `af` field.** Compute `ac/an`. Use `dataset: gnomad_r4`;
  `popmax` lives under `faf95`, not on the block. A bare `variant(id:...)` query fails —
  the argument is `variantId` plus a required `dataset`.
- **gnomAD version drift moves the number.** The same variant reported as AF 1.37e-06 from
  one release reads 1.24e-06 in `gnomad_r4` joint (2/1,613,494). Always record dataset + date.
- **VEP flags are camelCase and opt-in.** Lowercase flags are silently ignored. BayesDel is
  not a VEP plugin — it needs `dbNSFP=BayesDel_addAF_score`, and that param suppresses
  AlphaMissense on the same request, so fetch it in a second call and merge.
- **VEP coverage is per transcript × per position.** PolyPhen returning nothing for one
  codon is normal, not a command error. Filter with `transcript_id=` to the MANE Select
  isoform; without it a variant spans 100+ transcripts.
- **PubMed has a precision cliff.** `"TP53" AND (...) OR (...)` returned 18,103 hits of noise;
  tight `[tiab]` phrase queries over the exact protein change returned 0. Zero is a result:
  it means no one has published an assay. Never widen to gene-only and automate the reading.
- **NCBI rate-limits at 3 requests/second** without an API key. Batch jobs need sleep or a key.
- **VEP flaps.** Identical URLs return 200, 500 and 503 within seconds — measured at
  roughly one success in two. `_get()` skips retries on 4xx (a real answer) but retries
  5xx six times with jitter for VEP, so a run takes 1–2 min. A 500 is not a bad parameter.
- **A typo'd gnomAD id is indistinguishable from a rare variant.** `variant(variantId:)`
  answers `Variant not found` for a wrong position *and* for a genuinely unobserved allele.
  PM2 is pathogenic evidence, so it may only fire on an id that `consistency.ok` passed.
  Worked example: `17-7676154-C-T` looks plausible, is 1.9 kb off, resolved to "absent",
  and manufactured a PM2 that the real variant (G245S, AC=5, AF 3.1e-06) does not deserve.
- **VEP `allele_string` is transcript-oriented.** On a minus-strand gene it is the reverse
  complement of the gnomAD id: `c.1178A>T` is VEP `A/T` at 7669613 but gnomAD `T-A`.
  Complement when `strand == -1`, or every TP53 variant reads as a mismatch.
- **`dbNSFP` must not ride on the AlphaMissense request.** Same URL = HTTP 500 and, worse,
  AlphaMissense silently disappears from the merged result. Two URLs, then merge.
- **A swallowed dbNSFP failure reads as "no BayesDel coverage".** A 5xx on the flaky merge
  call leaves BayesDel `None`, which looks exactly like a position VEP has no dbNSFP data
  for. Since count3 counts only tools that returned data, this silently shrinks the
  denominator and can flip PP3/BP4. The payload now carries `vep.bayesdel_merge`
  (`ok` / `not returned` / `merge failed: ...`) — check it before quoting the "n/available"
  count, and re-run the fetch if it says the merge failed.
- **A guard that compares nothing passes everything.** The identity check originally parsed
  only `c.` HGVS, so a ClinVar record filed as `NC_000017.11:g.7579817G>C` left
  `clinvar_hgvs: null` and returned `ok: true` without comparing a single value. Any "no
  problems found" result must be backed by a list of checks actually performed — the payload
  now carries `consistency.checks` and reports `no cross-check was possible` when empty.
- **A wrong MANE transcript silently empties VEP.** Filtering `--transcript` to another
  gene's accession drops every consequence: `most_severe_consequence` becomes `"?"`, scores
  are `{}`, and downstream that reads as "this position has no predictors" instead of "you
  gave me the wrong transcript". The payload flags `vep.transcript_mismatch`. A one-character
  typo in an accession otherwise produces a confident wrong answer with no error.
- **CNV and SV records live in ClinVar `obj_type`, not `variant_type`.** `variant_type`
  is null in esummary, so a check against it never fires. Match `copy number gain`,
  `copy number loss`, `inversion`, `translocation`, `complex` and abort with exit 2.
- **Compare parsed HGVS coordinates, not raw strings.** ClinVar writes
  `NM_000546.6(TP53):c.733G>A` — accession outside the parens — so a string compare marks
  identical variants as mismatched and blocks a correct run.
- **Recompute every code from the numbers printed beside it.** Hand-written evidence
  tables routinely pair a code with a value that contradicts its own threshold — a variant
  with 2 alleles in gnomAD is not "absent", so PM2 cannot fire under an AF = 0 threshold.
  Either drop the code or declare an ultra-rare cutoff; never carry both. When auditing an
  existing table, check each code against its threshold and its printed value *before*
  adding new evidence, and report the contradiction rather than silently adopting the
  table's lean.

## Deliverable

An evidence matrix, one row per variant, columns: variant · ClinVar classification and
review status · gnomAD AF, grpmax, homozygotes · predictor scores and how many called
damaging · functional evidence · codes applied · lean. Every cell traceable to a fetched
value with source and date.

Methods sentence to carry: *"Variants were annotated against ClinVar, gnomAD (dataset
gnomad_r4, accessed YYYY-MM-DD) and Ensembl VEP (transcript NM_…/ENST…). In-silico evidence
was assessed on algorithms with annotation for that position; algorithms without data were
excluded and recorded as NA. PP3/BP4 was applied on the {count3|revel} scheme."*

## Files

- `scripts/fetch_evidence.py` — all four sources, stdlib only, NA never becomes 0
- `scripts/acmg_codes.py` — code derivation, auditable basis line per code
- `scripts/validate_skill.py` — Agent Skills spec conformance + host compatibility table
- `tests/test_engine.py` — offline: does PM2 still refuse an unverified id? (`python3 tests/test_engine.py`)
- `references/thresholds.md` — full tables, scheme comparison, two worked examples
  (D393V: scheme divergence; G245S: the identity trap that manufactures PM2)
- `examples/g245s.json` — a real fetched payload, used as the test fixture

Published as a standalone repo: https://github.com/darshak26/acmg-variant-curation
(version in `metadata.version`; see CHANGELOG.md — threshold or guard changes are breaking).
