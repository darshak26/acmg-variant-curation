# Changelog

Version lives in `metadata.version` in `SKILL.md`. Thresholds, schemes and guard
behaviour are part of the contract: changing them changes the codes a given variant
receives, so pin the tool version in Methods alongside the dataset and access date.

## [2.0.0] - 2026-10-26

**Breaking: the codes a variant receives can differ from v1.0.0.** The lean is now ClinGen SVI
points (Tavtigian 2020) instead of a code count, and the identifier guard can no longer pass
without comparing anything. Re-run any v1.0.0 curation you intend to quote.

### Changed
- **Lean: points, not votes.** A single Supporting code used to read "pathogenic-leaning"
  because one pathogenic code outnumbered zero benign ones. Any variant that was merely absent
  from gnomAD and had no predictors therefore looked pathogenic-leaning — synonymous, intronic
  and UTR variants included. Reported case: CD68 `NC_000017.11:g.7579817G>C`, ClinVar Benign,
  absent from gnomAD, CADD 6.68 — produced "Pathogenic-leaning". Now: 1 point, VUS, with the
  ClinVar conflict printed.
- **Direction is no longer called a classification.** The output states the point totals, the
  direction, and that direction is not a classification. `Verdict:` is printed as a stable line
  (`PATHOGENIC` / `BENIGN` / `VUS`).
- **The 6-point ceiling is now documented and tested.** The seven public-data codes cannot reach
  the 10-point pathogenic threshold, by construction. A test pins this so it is not "fixed" into
  over-calling.

### Fixed
- **The guard could be bypassed.** `check_identity` parsed only `c.` HGVS, so a ClinVar record
  filed as genomic (`NC_000017.11:g.7579817G>C`) left `clinvar_hgvs: null` and returned
  `ok: true` having compared nothing. It now parses `c.`/`g.`/`n.`, cross-checks ClinVar
  against VEP's resolved locus, and returns **not ok** when no cross-check was possible.
  `consistency.checks` lists the comparisons actually performed.
- **Accession prefix caused false mismatches.** Comparing raw HGVS strings marked
  `NM_000546.6(TP53):c.733G>A` against `c.733G>A` as a wrong variant. Comparison now uses
  parsed coordinates.
- **A wrong MANE transcript silently emptied VEP.** No consequence returned for the requested
  accession produced blank scores and read as "no predictors". Now flagged as
  `vep.transcript_mismatch`.
- **CNVs and structural variants are rejected, not half-curated.** ClinVar files them under
  `obj_type` (`variant_type` is null, so a check against that never fires). The run stops with
  exit code 2, the variant type and size, and a pointer to Riggs 2020.

### Added
- Genomic (`g.`) and non-coding (`n.`) HGVS support.
- Offline guard tests: genomic cross-check, wrong locus, nothing-to-compare, accession prefix,
  and out-of-scope rejection with a canned ClinVar response.

## [1.0.0] - 2026-10-26

First public release. Seven codes (BA1, BS1, PM2, PP3, BP4, PS3, BS3) derived from
ClinVar, gnomAD v4, Ensembl VEP and PubMed; two declared PP3 schemes (`count3`,
`revel`); identifier consistency guard; offline test suite; spec-conformance validator.

### Fixed during v1.0.0 testing

These were found by running the tool against real variants, not by reading it.

- **A typo'd gnomAD id silently became PM2.** gnomAD answers `Variant not found` for a
  wrong position *and* for a genuinely unobserved allele, so a mistyped id read as
  "absent" and manufactured pathogenic evidence. `acmg_codes.py` now refuses PM2 when
  `consistency.ok` is false. Demonstrated on TP53 p.Gly245Ser: the real variant has
  AC=5 (AF 3.1e-06) and does not deserve PM2.
- **No cross-check between identifiers.** ClinVar could report a classification for
  variant A while gnomAD and VEP answered for variant B. Added `check_identity()`,
  which compares the ClinVar title's HGVS against the requested HGVS and VEP's locus
  and alleles against the gnomAD id.
- **`dbNSFP` rode on the AlphaMissense request.** The same URL both returned HTTP 500
  and suppressed AlphaMissense, silently thinning the predictor set. Split into two
  calls and merged; BayesDel now appears correctly (4/4 tools on G245S).
- **Retries too thin for a flapping endpoint.** Ensembl VEP returns 200/500/503 within
  seconds on identical requests. 5xx is now retried 6x with jitter; 4xx fails fast
  instead of burning retries.
- **Minus-strand allele orientation.** VEP returns `allele_string` in transcript
  orientation, so on TP53 every variant compared as a mismatch until complemented
  against `strand == -1`.

### Added

- `scripts/validate_skill.py` — spec conformance and host-compatibility check.
- `tests/test_engine.py` — offline assertions that PM2 refuses unverified ids.
- `examples/g245s.json` — real fetched payload used as a test fixture.
- `compatibility` and `metadata.version` frontmatter fields.
