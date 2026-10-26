# Changelog

Version lives in `metadata.version` in `SKILL.md`. Thresholds, schemes and guard
behaviour are part of the contract: changing them changes the codes a given variant
receives, so pin the tool version in Methods alongside the dataset and access date.

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
